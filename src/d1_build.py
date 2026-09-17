"""Exp 1 (D1) data build: SEQC2 SomaticSeq ensemble tables -> annotated per-site tables.

For each sequencing run we produce, for SNVs and indels:
  * a UNION table: every site called (PASS-like) by >=1 caller in the run
  * a TRUTH table: every SEQC2 v1.2.1 high-confidence variant, with per-caller flags
    (0 for all callers if the site is absent from the run's ensemble table)
Both are annotated with GIAB v3.3 strata (GRCh38) and SEQC2 HC regions.

Outputs: results/d1/{run}_{vtype}_{union|truth}.parquet
"""
from __future__ import annotations

import logging
import time

import numpy as np
import pandas as pd
import pysam

from common import (HC_BED, MAIN_CHROMS, RES, SEQC2, TRUTH_INDEL, TRUTH_SNV, annotate_regions,
                    ensure_dirs, giab_beds)

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
log = logging.getLogger("d1_build")

OUT = RES / "d1"
CONS = SEQC2 / "somaticseq_consensus"
RUNS = {
    "IL_T_1": CONS / "IL_T_1_vs_IL_N_1" / "WGS.bwa.dedup-IL_T_1_vs_IL_N_1-Ensemble.{t}.tsv.gz",
    "NV_T_1": CONS / "NV_T_1_vs_NV_N_1" / "WGS.bwa.dedup-NV_T_1_vs_NV_N_1-Ensemble.{t}.tsv.gz",
    "FD_T_1": CONS / "FD_T_1_vs_FD_N_1" / "WGS.bwa.dedup-FD_T_1_vs_FD_N_1-Ensemble.{t}.tsv.gz",
    "NS_T_1": CONS / "NS_T_1_vs_NS_N_1" / "WGS.bwa.dedup-NS_T_1_vs_NS_N_1-Ensemble.{t}.tsv.gz",
    "NS_380X": CONS / "Combine_9_NS_380X" / "CombineNova_380X.bwa.Ensemble.{t}.tsv.gz",
}
SNV_CALLERS = ["MuTect2", "SomaticSniper", "VarDict", "MuSE", "Strelka"]
INDEL_CALLERS = ["MuTect2", "VarDict", "Strelka"]
LENIENT = ("VarDict", "SomaticSniper", "MuSE")  # callers with a separate lenient flag
SNIPER_MIN_SSC = 40
KEY = ["CHROM", "POS", "REF", "ALT"]


def load_ensemble(path, vtype: str) -> pd.DataFrame:
    cols = KEY + ["if_MuTect", "if_VarDict", "if_Strelka", "T_DP", "N_DP", "T_ALT_FOR", "T_ALT_REV",
                  "SiteHomopolymer_Length", "MaxHomopolymer_Length", "InDel_Length"]
    if vtype == "snv":
        cols += ["if_SomaticSniper", "Sniper_Score", "MuSE_Tier"]
    d = pd.read_csv(path, sep="\t", usecols=cols, low_memory=False)
    d = d[d.CHROM.isin(MAIN_CHROMS)].copy()
    # --- primary (strict, PASS-like) caller flags
    d["MuTect2"] = (d.if_MuTect == 1).astype(np.int8)
    d["VarDict"] = (d.if_VarDict == 1).astype(np.int8)
    d["Strelka"] = (d.if_Strelka == 1).astype(np.int8)
    # --- lenient flags (sensitivity analysis)
    d["VarDict_len"] = (d.if_VarDict >= 0.5).astype(np.int8)
    if vtype == "snv":
        d["SomaticSniper"] = ((d.if_SomaticSniper == 1) & (d.Sniper_Score >= SNIPER_MIN_SSC)).astype(np.int8)
        d["MuSE"] = (d.MuSE_Tier == 1).astype(np.int8)
        d["SomaticSniper_len"] = (d.if_SomaticSniper == 1).astype(np.int8)
        d["MuSE_len"] = (d.MuSE_Tier > 0).astype(np.int8)
    d["vaf_obs"] = (d.T_ALT_FOR + d.T_ALT_REV) / d.T_DP.replace(0, np.nan)
    drop = [c for c in d.columns if c.startswith("if_")] + ["T_ALT_FOR", "T_ALT_REV", "Sniper_Score", "MuSE_Tier"]
    return d.drop(columns=[c for c in drop if c in d.columns])


def load_truth(path) -> pd.DataFrame:
    recs = []
    with pysam.VariantFile(str(path)) as vf:
        for r in vf:
            info = r.info
            filt = list(r.filter.keys())
            conf = next((f for f in filt if f.endswith("Conf")), "NA")
            recs.append((r.chrom, r.pos, r.ref, r.alts[0], float(info.get("TVAF", np.nan)),
                         float(info.get("LC", np.nan)), conf, int(info.get("nPASSES", 0))))
    t = pd.DataFrame(recs, columns=KEY + ["TVAF", "LC_score", "conf", "nPASSES"])
    return t[t.CHROM.isin(MAIN_CHROMS)].reset_index(drop=True)


def main():
    ensure_dirs()
    OUT.mkdir(parents=True, exist_ok=True)
    beds = giab_beds("GRCh38")
    beds["hc"] = HC_BED
    truth = {"snv": load_truth(TRUTH_SNV), "indel": load_truth(TRUTH_INDEL)}
    for vt, t in truth.items():
        log.info("truth %s: %d variants, conf=%s", vt, len(t), t.conf.value_counts().to_dict())

    for vt, callers in [("snv", SNV_CALLERS), ("indel", INDEL_CALLERS)]:
        tag = "sSNV" if vt == "snv" else "sINDEL"
        len_callers = [c + "_len" if c in LENIENT else c for c in callers]
        tables = {}
        for run, pat in RUNS.items():
            t0 = time.time()
            d = load_ensemble(str(pat).format(t=tag), vt)
            d["n_callers"] = d[callers].sum(axis=1)
            d["n_callers_len"] = d[len_callers].sum(axis=1)
            tables[run] = d
            log.info("%s %s: %d rows, strict-union %d, lenient-union %d (%.1fs)", run, vt, len(d),
                     (d.n_callers > 0).sum(), (d.n_callers_len > 0).sum(), time.time() - t0)

        # Annotate all unique sites once (union across runs + truth)
        keep = pd.concat([d.loc[d.n_callers_len > 0, KEY] for d in tables.values()] + [truth[vt][KEY]])
        sites = keep.drop_duplicates().reset_index(drop=True)
        t0 = time.time()
        ann = annotate_regions(sites, beds)
        log.info("annotated %d unique %s sites in %.1fs; LC frac=%.3f", len(sites), vt, time.time() - t0,
                 ann.lc.mean())

        tr = truth[vt].merge(ann, on=KEY, how="left")
        tr["is_truth"] = True
        for run, d in tables.items():
            d = d.merge(tr[KEY + ["is_truth", "TVAF", "LC_score", "conf"]], on=KEY, how="left")
            d["is_truth"] = d.is_truth.fillna(False).astype(bool)
            # union tables (strict and lenient unions share the file; filter at analysis time)
            u = d[d.n_callers_len > 0].merge(ann, on=KEY, how="left")
            u.to_parquet(OUT / f"{run}_{vt}_union.parquet", index=False)
            # truth tables: left join truth -> run flags; absent sites = called by nobody
            cols = [c for c in d.columns if c not in ("is_truth", "TVAF", "LC_score", "conf")]
            tt = tr.merge(d[cols], on=KEY, how="left", indicator=True)
            tt["in_table"] = tt.pop("_merge") == "both"
            fill = callers + [c for c in len_callers if c not in callers] + ["n_callers", "n_callers_len"]
            tt[fill] = tt[fill].fillna(0).astype(np.int8)
            tt.to_parquet(OUT / f"{run}_{vt}_truth.parquet", index=False)
            log.info("%s %s: union=%d (truth among union=%d), truth rows=%d, truth absent from table=%d",
                     run, vt, len(u), int(u.is_truth.sum()), len(tt), int((~tt.in_table).sum()))


if __name__ == "__main__":
    main()

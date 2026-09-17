"""Exp 1b: robustness on raw IL_T_1 caller VCFs (PASS-only definitions; adds TNscope).

Calls:
  MuSE, MuTect2 (GATK4), Strelka snv/indel, TNscope : FILTER == PASS
  VarDict       : FILTER == PASS and INFO STATUS contains 'Somatic' (StrongSomatic/LikelySomatic)
  SomaticSniper : SS == 2 and SSC >= 40 (VCF has no FILTER)
SNV callers (6): MuTect2, SomaticSniper, VarDict, MuSE, Strelka, TNscope
Indel callers (4): MuTect2, VarDict, Strelka, TNscope

VAF: truth -> TVAF; union -> observed VAF from the IL_T_1 SomaticSeq table when the site is present,
otherwise the median caller-reported VAF.
Also: leave-one-caller-out (LOCO) all-agree to show how much the weakest caller drives disagreement.

Outputs: results/d1b/*.parquet, results/d1b_*.csv
"""
from __future__ import annotations

import logging
import re
from itertools import combinations
from multiprocessing import Pool

import numpy as np
import pandas as pd

from common import HC_BED, MAIN_CHROMS, RES, SEQC2, agreement_summary, annotate_regions, giab_beds
from d1_build import KEY
from d2_purity import fmt, normalize, stream_pass, truth_table
from hypothesis_tests import h1_h2, prep, within_bin_contrast

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
log = logging.getLogger("d1b")

RAW = SEQC2 / "raw_caller_vcfs" / "IL_T_1"
OUT = RES / "d1b"
PFX = "WGS.bwa.dedup-IL_T_1_vs_IL_N_1-"
FILES = {  # caller -> (files, tumor sample name or index)
    "MuSE": ([PFX + "MuSE.vcf.gz"], "TUMOR"),
    "MuTect2": ([PFX + "MuTect2.vcf.gz"], "WGS_IL_T_1"),
    "SomaticSniper": ([PFX + "SomaticSniper.vcf.gz"], "TUMOR"),
    "Strelka": ([PFX + "Strelka.snv.vcf.gz", PFX + "Strelka.indel.vcf.gz"], "TUMOR"),
    "TNscope": ([PFX + "TNScope.201711.02.vcf.gz"], "WGS_IL_T_1"),
    "VarDict": ([PFX + "VarDict.vcf.gz"], 0),
}
SNV_CALLERS = ["MuTect2", "SomaticSniper", "VarDict", "MuSE", "Strelka", "TNscope"]
INDEL_CALLERS = ["MuTect2", "VarDict", "Strelka", "TNscope"]


def caller_vaf(caller, f, ref, alt):
    try:
        if caller in ("MuTect2", "TNscope", "VarDict"):
            return float(f["AF"].split(",")[0])
        if caller == "MuSE":
            r, a = (int(x) for x in f["AD"].split(",")[:2]); return a / (r + a) if r + a else np.nan
        if caller == "SomaticSniper":
            d = [int(x) for x in f["DP4"].split(",")]; return (d[2] + d[3]) / sum(d) if sum(d) else np.nan
        if caller == "Strelka":
            if len(ref) == 1 and len(alt) == 1:
                r = int(f[f"{ref}U"].split(",")[0]); a = int(f[f"{alt}U"].split(",")[0])
            else:
                r = int(f["TAR"].split(",")[0]); a = int(f["TIR"].split(",")[0])
            return a / (r + a) if r + a else np.nan
    except (KeyError, ValueError, ZeroDivisionError):
        return np.nan
    return np.nan


def parse(caller):
    files, tname = FILES[caller]
    recs = []
    for fn in files:
        # SomaticSniper has no FILTER -> stream everything ('somaticSniper' key disables PASS filter)
        mode = "somaticSniper" if caller == "SomaticSniper" else caller
        tcol = None
        for x in stream_pass(RAW / fn, mode):
            if x[0].startswith("#"):
                s = [v.strip('"') for v in x[9:]]
                tcol = 9 + (tname if isinstance(tname, int) else s.index(tname))
                continue
            if x[0] not in MAIN_CHROMS:
                continue
            f = fmt(x[8], x[tcol])
            if caller == "SomaticSniper":
                if f.get("SS") != "2" or not f.get("SSC", ".").isdigit() or int(f["SSC"]) < 40:
                    continue
            if caller == "VarDict" and not re.search(r"STATUS=\w*Somatic", x[7]):
                continue
            ref, alt = x[3].upper(), x[4].upper().split(",")[0]
            if len(ref) == len(alt) and len(ref) > 1:
                continue  # MNP
            vaf = caller_vaf(caller, f, ref, alt)
            p, r, a = normalize(x[0], int(x[1]), ref, alt)
            recs.append((x[0], p, r, a, vaf))
    d = pd.DataFrame(recs, columns=KEY + [f"vaf_{caller}"]).drop_duplicates(KEY)
    d[caller] = np.int8(1)
    log.info("%s: %d PASS calls", caller, len(d))
    return d


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    f = OUT / "calls.parquet"
    if not f.exists():
        with Pool(len(FILES)) as pool:
            parts = pool.map(parse, list(FILES))
        wide = parts[0]
        for p in parts[1:]:
            wide = wide.merge(p, on=KEY, how="outer")
        for c in FILES:
            wide[c] = wide[c].fillna(0).astype(np.int8)
        wide["vtype"] = np.where((wide.REF.str.len() == 1) & (wide.ALT.str.len() == 1), "snv", "indel")
        # caller-independent observed VAF from the SomaticSeq table where available
        obs = []
        for vt, tag in (("snv", "sSNV"), ("indel", "sINDEL")):
            o = pd.read_parquet(RES / "d1" / f"IL_T_1_{vt}_union.parquet", columns=KEY + ["vaf_obs"])
            obs.append(o)
        obs = pd.concat(obs).drop_duplicates(KEY)
        wide = wide.merge(obs, on=KEY, how="left")
        wide["vaf_med"] = wide[[f"vaf_{c}" for c in FILES]].median(axis=1)
        wide["vaf_use"] = wide.vaf_obs.fillna(wide.vaf_med)
        wide.to_parquet(f, index=False)
    calls = pd.read_parquet(f)
    log.info("union sites: %d; with table VAF: %.3f", len(calls), calls.vaf_obs.notna().mean())

    truth = truth_table()
    sites = pd.concat([truth[KEY], calls[KEY]]).drop_duplicates().reset_index(drop=True)
    beds = giab_beds("GRCh38"); beds["hc"] = HC_BED
    ann = annotate_regions(sites, beds)
    calls = calls.merge(ann, on=KEY, how="left")
    tt = truth.merge(calls[KEY + list(FILES)], on=KEY, how="left").merge(ann, on=KEY, how="left")
    tt[list(FILES)] = tt[list(FILES)].fillna(0).astype(np.int8)
    tt.to_parquet(OUT / "truth.parquet", index=False)

    summ, h12, h3, h4, loco = [], [], [], [], []
    for denom, base, vcol in (("truth", tt, "TVAF"), ("union", calls, "vaf_use"), ("union_hc", calls[calls.hc], "vaf_use")):
        for cs_name in ("native", "shared4"):
            parts = []
            for vt in ("snv", "indel"):
                cols = INDEL_CALLERS if cs_name == "shared4" else (SNV_CALLERS if vt == "snv" else INDEL_CALLERS)
                sub = base[base.vtype == vt]
                if denom != "truth":
                    sub = sub[sub[cols].sum(axis=1) > 0]
                p = prep(sub, cols, vcol, vt)
                p["all"] = "all"
                for grp in (["all"], ["vaf_bin"], ["vaf_bin", "lc"]):
                    summ.append(agreement_summary(p, cols, grp).assign(vtype=vt, grouping="+".join(grp),
                                                                       denom=denom, callerset=cs_name))
                parts.append(p)
                if cs_name == "native":
                    # leave-one-caller-out (and best subset of K-1) within the same site set
                    for drop in cols:
                        keep = [c for c in cols if c != drop]
                        s = agreement_summary(p, keep, ["vaf_bin"])
                        loco.append(s[["vaf_bin", "n", "all_agree", "majority"]].assign(
                            vtype=vt, denom=denom, dropped=drop, K=len(keep)))
            both = pd.concat(parts, ignore_index=True)
            tag = dict(denom=denom, callerset=cs_name)
            h12.append(h1_h2(both, ["vtype"]).assign(**tag))
            h3.append(within_bin_contrast(both, "vtype", "indel", "snv", []).assign(**tag))
            h4.append(within_bin_contrast(both, "lc", True, False, ["vtype"]).assign(**tag))
    pd.concat(summ).to_csv(RES / "d1b_summary.csv", index=False)
    pd.concat(h12).to_csv(RES / "d1b_h1h2.csv", index=False)
    pd.concat(h3).to_csv(RES / "d1b_h3.csv", index=False)
    pd.concat(h4).to_csv(RES / "d1b_h4.csv", index=False)
    pd.concat(loco).to_csv(RES / "d1b_leave_one_out.csv", index=False)

    # pairwise Jaccard per VAF bin on the union (which caller pairs agree?)
    pw = []
    for vt in ("snv", "indel"):
        cols = SNV_CALLERS if vt == "snv" else INDEL_CALLERS
        p = prep(calls[(calls.vtype == vt) & (calls[cols].sum(axis=1) > 0)], cols, "vaf_use", vt)
        for b, g in p.groupby("vaf_bin", observed=True):
            for a, c in combinations(cols, 2):
                x, y = g[a].astype(bool), g[c].astype(bool)
                pw.append(dict(vtype=vt, vaf_bin=b, a=a, b=c, jaccard=(x & y).sum() / max(1, (x | y).sum())))
    pd.DataFrame(pw).to_csv(RES / "d1b_pairwise_jaccard.csv", index=False)
    log.info("done")


if __name__ == "__main__":
    main()

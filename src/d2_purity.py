"""Exp 2 (D2): SEQC2 tumor-purity dilution series (100x WGS, BWA) - controlled-VAF test.

Five callers per purity (TNscope, Lancet, MuTect2, SomaticSniper[SNV only], Strelka).
Calls = FILTER==PASS (SomaticSniper files are pre-filtered to SS=2 & SSC>=40, so every row is a call).
Indels are left-normalised against hg38 before exact matching.

Denominators:
  truth : SEQC2 v1.2.1 truth variants, expected VAF = TVAF x purity (caller-independent)
  union : all PASS calls; VAF = median of caller-reported VAFs (documented caveat)

Outputs: results/d2/calls_{purity}.parquet, results/d2_*.csv
"""
from __future__ import annotations

import logging
import subprocess
from multiprocessing import Pool
from functools import lru_cache

import numpy as np
import pandas as pd
from pyfaidx import Fasta
from statsmodels.stats.contingency_tables import mcnemar

from common import (DATA, MAIN_CHROMS, RES, SEQC2, agreement_summary, annotate_regions, giab_beds, HC_BED)
from d1_build import KEY, load_truth, TRUTH_SNV, TRUTH_INDEL
from hypothesis_tests import h1_h2, prep, within_bin_contrast

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
log = logging.getLogger("d2")

PDIR = SEQC2 / "purity_series_100X"
OUT = RES / "d2"
PURITY = {"1-0": 1.0, "3-1": 0.75, "1-1": 0.5, "1-4": 0.2, "1-9": 0.1, "1-19": 0.05}
CALLERS = ["TNscope", "lancet", "muTect2", "somaticSniper", "strelka"]
SNV_CALLERS = CALLERS
INDEL_CALLERS = ["TNscope", "lancet", "muTect2", "strelka"]
_FA = {}


def _fa():
    # opened lazily so each worker process has its own file handle
    if "hg38" not in _FA:
        _FA["hg38"] = Fasta(str(DATA / "reference" / "hg38.fa"), as_raw=True, sequence_always_upper=True)
    return _FA["hg38"]


@lru_cache(maxsize=None)
def _base(chrom, pos):  # 1-based single base
    return _fa()[chrom][pos - 1]


def normalize(chrom, pos, ref, alt):
    """vt-style normalisation: right-trim / left-extend until stable, then left-trim (keep 1 anchor)."""
    if len(ref) == 1 and len(alt) == 1:
        return pos, ref, alt
    changed = True
    while changed:
        changed = False
        if ref and alt and ref[-1] == alt[-1] and not (len(ref) == 1 and len(alt) == 1):
            ref, alt = ref[:-1], alt[:-1]
            changed = True
        if len(ref) == 0 or len(alt) == 0:
            pos -= 1
            b = _base(chrom, pos)
            ref, alt = b + ref, b + alt
            changed = True
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref, alt, pos = ref[1:], alt[1:], pos + 1
    return pos, ref, alt


def stream_pass(path, caller):
    """Yield data lines that are calls (PASS; all lines for pre-filtered SomaticSniper)."""
    cmd = f"zcat {path} | grep -v '^##'"
    if caller != "somaticSniper":
        cmd += " | awk -F'\\t' '$1 ~ /^#/ || $7==\"PASS\"'"
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
    for line in p.stdout:
        yield line.rstrip("\n").split("\t")
    p.wait()


def fmt(fmt_str, sample_str):
    return dict(zip(fmt_str.split(":"), sample_str.split(":")))


def tumor_vaf(caller, f, ref, alt):
    try:
        if caller in ("TNscope", "muTect2"):
            return float(f["AF"].split(",")[0])
        if caller == "lancet":
            r, a = (int(x) for x in f["AD"].split(",")[:2])
            return a / (r + a) if r + a else np.nan
        if caller == "somaticSniper":
            d = [int(x) for x in f["DP4"].split(",")]
            s = sum(d)
            return (d[2] + d[3]) / s if s else np.nan
        if caller == "strelka":
            if len(ref) == 1 and len(alt) == 1:
                r = int(f[f"{ref}U"].split(",")[0]); a = int(f[f"{alt}U"].split(",")[0])
            else:
                r = int(f["TAR"].split(",")[0]); a = int(f["TIR"].split(",")[0])
            return a / (r + a) if r + a else np.nan
    except (KeyError, ValueError):
        return np.nan
    return np.nan


def parse_caller(purity_code, caller):
    path = PDIR / f"SPP_GT_{purity_code}_100X.bwa.{caller}.vcf.gz"
    tcol = None
    recs = []
    n_mnp = 0
    for x in stream_pass(path, caller):
        if x[0].startswith("#"):
            samples = x[9:]
            if caller in ("TNscope", "lancet"):
                tcol = 9 + next(i for i, s in enumerate(samples) if f"_{purity_code}_" in s)
            else:
                tcol = 9 + samples.index("TUMOR")
            continue
        chrom = x[0]
        if chrom not in MAIN_CHROMS:
            continue
        pos, ref = int(x[1]), x[3].upper()
        for alt in x[4].upper().split(","):
            if len(ref) == len(alt) and len(ref) > 1:
                n_mnp += 1
                continue
            vaf = tumor_vaf(caller, fmt(x[8], x[tcol]), ref, alt)
            p2, r2, a2 = normalize(chrom, pos, ref, alt)
            recs.append((chrom, p2, r2, a2, vaf))
            break  # first ALT only
    d = pd.DataFrame(recs, columns=KEY + ["vaf"]).drop_duplicates(KEY)
    log.info("  %s %s: %d calls (skipped MNPs %d)", purity_code, caller, len(d), n_mnp)
    return d


def build_one(code):
    f = OUT / f"calls_{code}.parquet"
    if f.exists():
        return
    wide = None
    for c in CALLERS:
        d = parse_caller(code, c).rename(columns={"vaf": f"vaf_{c}"})
        d[c] = np.int8(1)
        wide = d if wide is None else wide.merge(d, on=KEY, how="outer")
    for c in CALLERS:
        wide[c] = wide[c].fillna(0).astype(np.int8)
    wide["vtype"] = np.where((wide.REF.str.len() == 1) & (wide.ALT.str.len() == 1), "snv", "indel")
    wide["vaf_med"] = wide[[f"vaf_{c}" for c in CALLERS]].median(axis=1)
    wide.to_parquet(f, index=False)
    log.info("purity %s: %d union sites", code, len(wide))


def build():
    OUT.mkdir(parents=True, exist_ok=True)
    with Pool(len(PURITY)) as pool:
        pool.map(build_one, list(PURITY))


def truth_table():
    t = []
    for vt, path in (("snv", TRUTH_SNV), ("indel", TRUTH_INDEL)):
        d = load_truth(path)
        if vt == "indel":
            n = [normalize(c, p, r, a) for c, p, r, a in d[KEY].itertuples(index=False)]
            d["POS"], d["REF"], d["ALT"] = zip(*n)
        d["vtype"] = vt
        t.append(d)
    return pd.concat(t, ignore_index=True)


def main():
    build()
    beds = giab_beds("GRCh38"); beds["hc"] = HC_BED
    truth = truth_table()
    calls = {code: pd.read_parquet(OUT / f"calls_{code}.parquet") for code in PURITY}
    sites = pd.concat([truth[KEY]] + [c[KEY] for c in calls.values()]).drop_duplicates().reset_index(drop=True)
    ann = annotate_regions(sites, beds)
    log.info("annotated %d sites", len(ann))

    summ, h12, h3, h4, match = [], [], [], [], []
    truth_long = []
    for code, pur in PURITY.items():
        c = calls[code].merge(ann, on=KEY, how="left")
        c["is_truth"] = c.set_index(KEY).index.isin(truth.set_index(KEY).index)
        tt = truth.merge(calls[code][KEY + CALLERS], on=KEY, how="left").merge(ann, on=KEY, how="left")
        tt[CALLERS] = tt[CALLERS].fillna(0).astype(np.int8)
        tt["vaf_exp"] = tt.TVAF * pur
        tt["purity"] = pur
        truth_long.append(tt)
        for vt in ("snv", "indel"):
            match.append(dict(purity=pur, vtype=vt, truth_n=int((tt.vtype == vt).sum()),
                              truth_called_any=int(((tt.vtype == vt) & (tt[CALLERS].sum(axis=1) > 0)).sum())))
        for denom, base, vcol in (("truth", tt, "vaf_exp"), ("union", c, "vaf_med"), ("union_hc", c[c.hc], "vaf_med")):
            for cs_name, cs in (("native", None), ("shared4", INDEL_CALLERS)):
                parts = []
                for vt in ("snv", "indel"):
                    cols = cs or (SNV_CALLERS if vt == "snv" else INDEL_CALLERS)
                    sub = base[base.vtype == vt]
                    if denom != "truth":
                        sub = sub[sub[cols].sum(axis=1) > 0]
                    p = prep(sub, cols, vcol, vt)
                    p["all"] = "all"
                    for grp in (["all"], ["vaf_bin"], ["vaf_bin", "lc"]):
                        summ.append(agreement_summary(p, cols, grp).assign(
                            vtype=vt, grouping="+".join(grp), purity=pur, denom=denom, callerset=cs_name))
                    parts.append(p)
                both = pd.concat(parts, ignore_index=True)
                tag = dict(purity=pur, denom=denom, callerset=cs_name)
                h12.append(h1_h2(both, ["vtype"]).assign(**tag))
                h3.append(within_bin_contrast(both, "vtype", "indel", "snv", []).assign(**tag))
                h4.append(within_bin_contrast(both, "lc", True, False, ["vtype"]).assign(**tag))
        log.info("analysed purity %s", code)

    pd.concat(summ).to_csv(RES / "d2_summary.csv", index=False)
    pd.concat(h12).to_csv(RES / "d2_h1h2.csv", index=False)
    pd.concat(h3).to_csv(RES / "d2_h3.csv", index=False)
    pd.concat(h4).to_csv(RES / "d2_h4.csv", index=False)
    pd.DataFrame(match).to_csv(RES / "d2_truth_matching.csv", index=False)

    # Pooled-over-purity truth analysis by expected VAF (same variants, VAF varied by dilution)
    tl = pd.concat(truth_long, ignore_index=True)
    tl.to_parquet(OUT / "truth_long.parquet", index=False)
    pooled = []
    for vt in ("snv", "indel"):
        cols = SNV_CALLERS if vt == "snv" else INDEL_CALLERS
        p = prep(tl[tl.vtype == vt], cols, "vaf_exp", vt)
        p["all"] = "all"
        pooled.append(agreement_summary(p, cols, ["vaf_bin", "lc"]).assign(vtype=vt, callerset="native"))
        pooled.append(agreement_summary(p, cols, ["vaf_bin"]).assign(vtype=vt, callerset="native"))
        p = prep(tl[tl.vtype == vt], INDEL_CALLERS, "vaf_exp", vt)
        pooled.append(agreement_summary(p, INDEL_CALLERS, ["vaf_bin", "lc"]).assign(vtype=vt, callerset="shared4"))
        pooled.append(agreement_summary(p, INDEL_CALLERS, ["vaf_bin"]).assign(vtype=vt, callerset="shared4"))
    pd.concat(pooled).to_csv(RES / "d2_pooled_expected_vaf.csv", index=False)

    # Paired within-variant test: clonal truth variants (TVAF>=0.4): all-agree at 100% vs each purity
    paired = []
    for vt in ("snv", "indel"):
        cols = SNV_CALLERS if vt == "snv" else INDEL_CALLERS
        w = tl[(tl.vtype == vt) & (tl.TVAF >= 0.4)].copy()
        w["agree"] = (w[cols].sum(axis=1) == len(cols)).astype(int)
        piv = w.pivot_table(index=KEY, columns="purity", values="agree")
        for pur in PURITY.values():
            if pur == 1.0:
                continue
            a, b = piv[1.0], piv[pur]
            tab = pd.crosstab(a, b).reindex(index=[0, 1], columns=[0, 1], fill_value=0).values
            res = mcnemar(tab, exact=False, correction=True)
            paired.append(dict(vtype=vt, purity=pur, n=len(piv), agree_100=a.mean(), agree_p=b.mean(),
                               expected_vaf_median=float(w[w.purity == pur].vaf_exp.median()),
                               lost=int(tab[1, 0]), gained=int(tab[0, 1]), mcnemar_p=float(res.pvalue)))
    pd.DataFrame(paired).to_csv(RES / "d2_paired_clonal.csv", index=False)
    log.info("done")


if __name__ == "__main__":
    main()

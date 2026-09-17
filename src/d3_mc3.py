"""Exp 3 (D3): TCGA MC3 pan-cancer replication (10,295 exomes, GRCh37).

Callers: SNV = MUTECT, MUSE, VARSCANS, SOMATICSNIPER, RADIA ; indel = INDELOCATOR, VARSCANI, PINDEL.
CENTERS entries with a trailing '*' are counted as calls in the primary analysis (consistent with
NCALLERS); a strict variant counts only un-starred entries.

To avoid mixing different caller-set sizes, each variant type is analysed only in samples where
ALL callers for that type were active (contributed >=1 call of that type anywhere in the MAF).

CAVEAT: the public MAF only contains variants with NCALLERS>=2, so agreement here is conditional on
at least two callers already agreeing (an upper bound relative to an untruncated union).

Outputs: results/d3_*.csv
"""
from __future__ import annotations

import logging

import numpy as np
import pandas as pd

from common import DATA, MAIN_CHROMS_37, RES, SEED, agreement_summary, annotate_regions, giab_beds, vaf_bin
from hypothesis_tests import h1_h2, within_bin_contrast

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
log = logging.getLogger("d3")

MAF = DATA / "tcga_mc3" / "mc3.v0.2.8.PUBLIC.maf.gz"
SNV_CALLERS = ["MUTECT", "MUSE", "VARSCANS", "SOMATICSNIPER", "RADIA"]
INDEL_CALLERS = ["INDELOCATOR", "VARSCANI", "PINDEL"]
COLS = ["Chromosome", "Start_Position", "End_Position", "Variant_Type", "Reference_Allele",
        "Tumor_Seq_Allele2", "Tumor_Sample_Barcode", "t_depth", "t_alt_count", "FILTER", "CENTERS", "NCALLERS"]
N_BOOT = 200


def load():
    d = pd.read_csv(MAF, sep="\t", usecols=COLS, low_memory=False, dtype={"Chromosome": str})
    log.info("MAF rows: %d", len(d))
    d = d[d.Chromosome.isin(MAIN_CHROMS_37)]
    d["vtype"] = np.where(d.Variant_Type == "SNP", "snv", np.where(d.Variant_Type.isin(["INS", "DEL"]), "indel", "other"))
    d = d[d.vtype != "other"].copy()
    cent = d.CENTERS.fillna("").str.split("|")
    for c in SNV_CALLERS + INDEL_CALLERS:
        d[c] = cent.map(lambda xs, c=c: int(c in xs or c + "*" in xs)).astype(np.int8)
        d[c + "_strict"] = cent.map(lambda xs, c=c: int(c in xs)).astype(np.int8)
    d["vaf"] = d.t_alt_count / d.t_depth.replace(0, np.nan)
    return d


def active_samples(d, vt, callers):
    s = d[d.vtype == vt].groupby("Tumor_Sample_Barcode")[callers].max()
    return set(s.index[(s == 1).all(axis=1)]), len(s)


def cluster_boot(p, callers, rng):
    """Sample-level cluster bootstrap CI of all-agree per VAF bin (variants are clustered in tumours)."""
    samples = p.Tumor_Sample_Barcode.unique()
    g = p.groupby(["Tumor_Sample_Barcode", "vaf_bin"], observed=True).all_agree.agg(["sum", "size"])
    g = g.unstack("vaf_bin", fill_value=0)
    sums, sizes = g["sum"], g["size"]
    out = []
    for _ in range(N_BOOT):
        idx = rng.choice(len(sums), len(sums), replace=True)
        out.append(sums.iloc[idx].sum() / sizes.iloc[idx].sum())
    b = pd.DataFrame(out)
    return pd.DataFrame({"boot_lo": b.quantile(0.025), "boot_hi": b.quantile(0.975)})


def main():
    rng = np.random.default_rng(SEED)
    d = load()
    # region annotation on unique spans (MAF: Start..End 1-based inclusive; INS has Start/End flanking)
    span = pd.DataFrame({"chrom": d.Chromosome.values, "start": d.Start_Position.values - 1,
                         "end": np.maximum(d.End_Position.values, d.Start_Position.values)})
    u = span.drop_duplicates().reset_index(drop=True)
    ann = annotate_regions(u, giab_beds("GRCh37"), span=u)
    d = _merge_ann(d, span, ann)
    log.info("annotated; LC fraction snv=%.3f indel=%.3f", d[d.vtype == "snv"].lc.mean(), d[d.vtype == "indel"].lc.mean())

    info, summ, h12, h3, h4, boots, macro = [], [], [], [], [], [], []
    for strict in (False, True):
        sfx = "_strict" if strict else ""
        parts = []
        for vt, callers in (("snv", SNV_CALLERS), ("indel", INDEL_CALLERS)):
            act, n_samp = active_samples(d, vt, callers)
            cols = [c + sfx for c in callers]
            for flt in ("PASS", "all"):
                sub = d[(d.vtype == vt) & d.Tumor_Sample_Barcode.isin(act)]
                if flt == "PASS":
                    sub = sub[sub.FILTER == "PASS"]
                sub = sub[sub.vaf.notna()]
                n_c = sub[cols].sum(axis=1)
                p = pd.DataFrame({"all_agree": (n_c == len(cols)).astype(int).values,
                                  "vaf": sub.vaf.values, "vtype": vt, "lc": sub.lc.values.astype(bool),
                                  "lowmap": sub.lowmap.values, "segdup": sub.segdup.values,
                                  "Tumor_Sample_Barcode": sub.Tumor_Sample_Barcode.values})
                for c, cc in zip(callers, cols):
                    p[c] = sub[cc].values
                p["vaf_bin"] = vaf_bin(p.vaf)
                p["all"] = "all"
                tag = dict(filter=flt, centers="strict" if strict else "starred_counted")
                info.append(dict(**tag, vtype=vt, samples_with_calls=n_samp, samples_all_callers_active=len(act),
                                 variants=len(p), frac_zero_callers=float((n_c == 0).mean()),
                                 frac_one_caller=float((n_c == 1).mean())))
                for grp in (["all"], ["vaf_bin"], ["lc"], ["vaf_bin", "lc"]):
                    summ.append(agreement_summary(p, callers, grp).assign(vtype=vt, grouping="+".join(grp), **tag))
                if flt == "PASS":
                    parts.append(p)
                    if not strict:
                        bt = cluster_boot(p, callers, rng).reset_index().rename(columns={"index": "vaf_bin"})
                        boots.append(bt.assign(vtype=vt, **tag))
                        # macro average: per-sample all-agree (samples with >=10 variants in the bin)
                        ps = p.groupby(["Tumor_Sample_Barcode", "vaf_bin"], observed=True).all_agree.agg(["mean", "size"])
                        ps = ps[ps["size"] >= 10].reset_index()
                        m = ps.groupby("vaf_bin", observed=True)["mean"].agg(["median", "mean", "size"]).reset_index()
                        macro.append(m.assign(vtype=vt, **tag))
        both = pd.concat(parts, ignore_index=True)
        tag = dict(filter="PASS", centers="strict" if strict else "starred_counted")
        h12.append(h1_h2(both, ["vtype"]).assign(**tag))
        h3.append(within_bin_contrast(both, "vtype", "indel", "snv", []).assign(**tag))
        h4.append(within_bin_contrast(both, "lc", True, False, ["vtype"]).assign(**tag))
        log.info("done strict=%s", strict)

    pd.DataFrame(info).to_csv(RES / "d3_info.csv", index=False)
    pd.concat(summ).to_csv(RES / "d3_summary.csv", index=False)
    pd.concat(h12).to_csv(RES / "d3_h1h2.csv", index=False)
    pd.concat(h3).to_csv(RES / "d3_h3.csv", index=False)
    pd.concat(h4).to_csv(RES / "d3_h4.csv", index=False)
    pd.concat(boots).to_csv(RES / "d3_cluster_bootstrap.csv", index=False)
    pd.concat(macro).to_csv(RES / "d3_macro_per_sample.csv", index=False)
    log.info("wrote results")


def _merge_ann(d, span, ann):
    """Attach region flags (computed on unique spans) back to every MAF row."""
    key = ["chrom", "start", "end"]
    a = ann[key + [c for c in ann.columns if c not in key]]
    s = span.reset_index(drop=True).merge(a, on=key, how="left")
    d = d.reset_index(drop=True)
    for c in a.columns:
        if c not in key:
            d[c] = s[c].values
    return d


if __name__ == "__main__":
    main()

"""Exp 1 (D1) analysis: stratified agreement metrics + H1-H4 tests + GLMs on SEQC2 runs.

Configurations analysed (each gives a complete set of tables):
  denom   : truth (VAF = pooled TVAF), truth_obsvaf (VAF observed in the run; sites present in table),
            union (all sites called by >=1 caller; VAF observed), union_hc (union within SEQC2 HC regions)
  callers : native (SNV 5 callers, indel 3 callers) | shared3 (MuTect2, VarDict, Strelka for both types)
  flags   : strict (PASS-like, primary) | lenient (VarDict 0.5, any MuSE tier, any SomaticSniper SS=2)

Outputs in results/: d1_summary.csv, d1_h1h2.csv, d1_h3.csv, d1_h4.csv, d1_highest_cell.csv, d1_glm.csv
"""
from __future__ import annotations

import logging

import numpy as np
import pandas as pd

from common import RES, agreement_summary
from d1_build import INDEL_CALLERS, LENIENT, OUT, RUNS, SNV_CALLERS
from hypothesis_tests import fit_glm, h1_h2, prep, within_bin_contrast

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
log = logging.getLogger("d1_analyze")

SHARED3 = ["MuTect2", "VarDict", "Strelka"]
CONFIGS = (
    [("strict", "native", d) for d in ("truth", "truth_obsvaf", "union", "union_hc")]
    + [("strict", "shared3", d) for d in ("truth", "union", "union_hc")]
    + [("lenient", "native", d) for d in ("truth", "union")]
)
GLM_COVARS = "np.log1p(depth) + lowmap + segdup"


def caller_cols(vt, callerset, flags):
    base = SHARED3 if callerset == "shared3" else (SNV_CALLERS if vt == "snv" else INDEL_CALLERS)
    return [c + "_len" if (flags == "lenient" and c in LENIENT) else c for c in base]


def load_config(run, vt, flags, callerset, denom):
    kind = "truth" if denom.startswith("truth") else "union"
    d = pd.read_parquet(OUT / f"{run}_{vt}_{kind}.parquet")
    cols = caller_cols(vt, callerset, flags)
    if kind == "union":
        d = d[d[cols].sum(axis=1) > 0]          # union for this caller set / flag definition
        if denom == "union_hc":
            d = d[d.hc]
        vaf_col = "vaf_obs"
    else:
        if denom == "truth_obsvaf":
            d = d[d.in_table]
            vaf_col = "vaf_obs"
        else:
            vaf_col = "TVAF"
            # depth for GLM: sites absent from the table have no depth -> use median
            d = d.copy()
            d["T_DP"] = d["T_DP"].fillna(d["T_DP"].median())
    p = prep(d, cols, vaf_col, vt)
    return p, [c.replace("_len", "") for c in cols]


def main():
    summ, h12, h3, h4, hi_cell, glms = [], [], [], [], [], []
    for flags, callerset, denom in CONFIGS:
        for run in RUNS:
            tag = dict(run=run, flags=flags, callerset=callerset, denom=denom)
            parts = {}
            for vt in ("snv", "indel"):
                p, cols = load_config(run, vt, flags, callerset, denom)
                parts[vt] = p
                p["all"] = "all"
                for grp in (["all"], ["vaf_bin"], ["lc"], ["vaf_bin", "lc"]):
                    s = agreement_summary(p, cols, grp)
                    s = s.assign(vtype=vt, grouping="+".join(grp), **tag)
                    summ.append(s)
            both = pd.concat(parts.values(), ignore_index=True)
            both["is_indel"] = (both.vtype == "indel").astype(int)
            both["cell"] = both.vtype + np.where(both.lc, "_LC", "_nonLC")

            h12.append(h1_h2(both, ["vtype"]).assign(**tag))
            h12.append(h1_h2(both.assign(vtype="pooled"), ["vtype"]).assign(**tag))
            h3.append(within_bin_contrast(both, "vtype", "indel", "snv", []).assign(**tag))
            h4.append(within_bin_contrast(both, "lc", True, False, ["vtype"]).assign(**tag))
            # which (type x LC) cell has the highest disagreement in each VAF bin?
            c = (both.groupby(["vaf_bin", "cell"], observed=True)
                 .agg(n=("all_agree", "size"), disagree=("all_agree", lambda x: 1 - x.mean()))
                 .reset_index())
            c = c[c.n >= 20]
            top = c.loc[c.groupby("vaf_bin", observed=True).disagree.idxmax()]
            hi_cell.append(top.assign(**tag))

            # GLMs (strict flags only; keeps runtime bounded)
            if flags != "strict" or denom == "truth_obsvaf":
                continue
            for vt, p in parts.items():
                p = p.copy()
                p["lc"] = p.lc.astype(int); p["lowmap"] = p.lowmap.astype(int); p["segdup"] = p.segdup.astype(int)
                if p.all_agree.nunique() < 2:
                    continue
                f = f"all_agree ~ C(vaf_bin, Treatment('>=20%')) * lc + {GLM_COVARS}"
                r = f"all_agree ~ C(vaf_bin, Treatment('>=20%')) + lc + {GLM_COVARS}"
                try:
                    g = fit_glm(p, f, r)
                    g0 = fit_glm(p, r)
                except Exception as e:  # separation in tiny truth-indel cells etc.
                    log.warning("GLM failed %s %s: %s", tag, vt, e)
                    continue
                glms.append(dict(**tag, model=f"{vt}_lcXvaf", n=g["n"], lrt_interaction_p=g["lrt_p"],
                                 lrt_df=g["lrt_df"], OR_lc_main=g0["table"].loc["lc", "OR"],
                                 OR_lc_lo=g0["table"].loc["lc", "OR_lo"], OR_lc_hi=g0["table"].loc["lc", "OR_hi"],
                                 p_lc_main=g0["table"].loc["lc", "p"]))
            if callerset == "shared3":
                b = both.copy()
                for col in ("lc", "lowmap", "segdup"):
                    b[col] = b[col].astype(int)
                f = (f"all_agree ~ C(vaf_bin, Treatment('>=20%')) * is_indel + "
                     f"C(vaf_bin, Treatment('>=20%')) * lc + {GLM_COVARS}")
                r = f"all_agree ~ C(vaf_bin, Treatment('>=20%')) + is_indel + C(vaf_bin, Treatment('>=20%')) * lc + {GLM_COVARS}"
                r0 = f"all_agree ~ C(vaf_bin, Treatment('>=20%')) + is_indel + lc + {GLM_COVARS}"
                try:
                    g = fit_glm(b, f, r)
                    g0 = fit_glm(b, r0)
                    glms.append(dict(**tag, model="pooled_indelXvaf", n=g["n"], lrt_interaction_p=g["lrt_p"],
                                     lrt_df=g["lrt_df"], OR_indel_main=g0["table"].loc["is_indel", "OR"],
                                     OR_indel_lo=g0["table"].loc["is_indel", "OR_lo"],
                                     OR_indel_hi=g0["table"].loc["is_indel", "OR_hi"],
                                     OR_lc_main=g0["table"].loc["lc", "OR"],
                                     OR_lc_lo=g0["table"].loc["lc", "OR_lo"], OR_lc_hi=g0["table"].loc["lc", "OR_hi"]))
                except Exception as e:
                    log.warning("pooled GLM failed %s: %s", tag, e)
            log.info("done %s", tag)

    pd.concat(summ).to_csv(RES / "d1_summary.csv", index=False)
    pd.concat(h12).to_csv(RES / "d1_h1h2.csv", index=False)
    pd.concat(h3).to_csv(RES / "d1_h3.csv", index=False)
    pd.concat(h4).to_csv(RES / "d1_h4.csv", index=False)
    pd.concat(hi_cell).to_csv(RES / "d1_highest_cell.csv", index=False)
    pd.DataFrame(glms).to_csv(RES / "d1_glm.csv", index=False)
    log.info("wrote results")


if __name__ == "__main__":
    main()

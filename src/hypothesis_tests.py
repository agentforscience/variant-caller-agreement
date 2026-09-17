"""Reusable hypothesis tests (H1-H4) and GLMs on a per-site agreement table.

Input table conventions (one row per variant site):
  vaf_bin (categorical, VAF_LABELS), vtype ('snv'/'indel'), lc (bool),
  all_agree (0/1), plus optional covariates depth, lowmap, segdup.
"""
from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy import stats

from common import VAF_LABELS, binom_one_sided, holm, newcombe_diff, wilson


def h1_h2(df: pd.DataFrame, by: list[str]) -> pd.DataFrame:
    """H1: all-agree > 0.80 in >=20% bin; H2: all-agree < 0.50 in <5% bin (exact one-sided)."""
    rows = []
    for key, g in df.groupby(by, observed=True):
        key = key if isinstance(key, tuple) else (key,)
        for hyp, b, p0, alt in [("H1", ">=20%", 0.80, "greater"), ("H2", "<5%", 0.50, "less")]:
            s = g[g.vaf_bin == b]
            k, n = int(s.all_agree.sum()), len(s)
            p, lo, hi = wilson(k, n)
            pv = binom_one_sided(k, n, p0, alt)
            rows.append({**dict(zip(by, key)), "hypothesis": hyp, "vaf_bin": b, "threshold": p0,
                         "n": n, "all_agree": p, "lo": lo, "hi": hi, "p_one_sided": pv,
                         "supported": bool(n > 0 and pv < 0.05)})
    out = pd.DataFrame(rows)
    out["p_holm"] = holm(out.p_one_sided)
    out["supported_holm"] = out.p_holm < 0.05
    return out


def within_bin_contrast(df: pd.DataFrame, factor: str, level_hi, level_lo, by: list[str]) -> pd.DataFrame:
    """Disagreement(level_hi) - disagreement(level_lo) within each VAF bin (Newcombe CI, chi2 test).

    Disagreement = 1 - all_agree. A positive difference supports H3/H4 in that bin.
    """
    rows = []
    for key, g in df.groupby(by + ["vaf_bin"], observed=True):
        key = key if isinstance(key, tuple) else (key,)
        a = g[g[factor] == level_hi]
        b = g[g[factor] == level_lo]
        ka, na = int((1 - a.all_agree).sum()), len(a)
        kb, nb = int((1 - b.all_agree).sum()), len(b)
        if na == 0 or nb == 0:
            d = lo = hi = pv = np.nan
        else:
            d, lo, hi = newcombe_diff(ka, na, kb, nb)
            tab = np.array([[ka, na - ka], [kb, nb - kb]])
            if (tab.sum(0) == 0).any():
                pv = 1.0
            elif tab.min() < 5:
                pv = stats.fisher_exact(tab)[1]
            else:
                pv = stats.chi2_contingency(tab)[1]
        rows.append({**dict(zip(by + ["vaf_bin"], key)), "factor": factor,
                     "n_hi": na, "disc_hi": ka / na if na else np.nan,
                     "n_lo": nb, "disc_lo": kb / nb if nb else np.nan,
                     "diff": d, "diff_lo": lo, "diff_hi": hi, "p": pv})
    out = pd.DataFrame(rows)
    out["p_holm"] = holm(out.p)
    out["direction_supports"] = out["diff"] > 0
    out["sig_supports"] = (out["diff_lo"] > 0) & (out.p_holm < 0.05)
    return out


def fit_glm(df: pd.DataFrame, formula: str, reduced: str | None = None) -> dict:
    """Logistic GLM with odds ratios; optional LRT vs a reduced (no-interaction) model."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        full = smf.glm(formula, data=df, family=sm.families.Binomial()).fit()
        res = {"formula": formula, "n": int(full.nobs), "aic": full.aic}
        tab = pd.DataFrame({"coef": full.params, "se": full.bse, "p": full.pvalues})
        ci = full.conf_int()
        tab["OR"] = np.exp(tab.coef)
        tab["OR_lo"] = np.exp(ci[0])
        tab["OR_hi"] = np.exp(ci[1])
        res["table"] = tab
        if reduced:
            red = smf.glm(reduced, data=df, family=sm.families.Binomial()).fit()
            lr = 2 * (full.llf - red.llf)
            dfd = full.df_model - red.df_model
            res["lrt_stat"] = lr
            res["lrt_df"] = dfd
            res["lrt_p"] = float(stats.chi2.sf(lr, dfd)) if dfd > 0 else np.nan
    return res


def prep(df: pd.DataFrame, callers: list[str], vaf_col: str, vtype: str) -> pd.DataFrame:
    """Standardise a per-site table for the tests above."""
    from common import vaf_bin
    out = pd.DataFrame({
        "vaf": df[vaf_col].values,
        "all_agree": (df[callers].sum(axis=1).values == len(callers)).astype(int),
        "n_callers": df[callers].sum(axis=1).values,
        "vtype": vtype,
        "lc": df["lc"].astype(bool).values,
    })
    for c in callers:  # canonical caller names (lenient "_len" suffix dropped)
        out[c.replace("_len", "")] = df[c].values
    for c in ("lowmap", "segdup", "hc", "homopol_ge7", "tandem_rep", "is_truth"):
        if c in df:
            out[c] = df[c].astype(bool).values
    if "T_DP" in df:
        out["depth"] = df["T_DP"].values
    out = out[out.vaf.notna()].copy()
    out["vaf_bin"] = pd.Categorical(vaf_bin(out.vaf), categories=VAF_LABELS)
    return out

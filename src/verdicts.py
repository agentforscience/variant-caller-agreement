"""Consolidate H1-H4 evidence across experiments into results/verdicts.csv (+ markdown tables).

H1/H2 : one-sided exact binomial test of all-agree vs 0.80 (>=20% bin) / 0.50 (<5% bin).
H3/H4 : within-VAF-bin disagreement differences (indel - SNV; LC - non-LC). "all bins" means the
        difference is positive AND its 95% Newcombe CI excludes 0 (Holm-adjusted p<0.05) in every bin
        with data (both groups n>=20).
"""
from __future__ import annotations

import pandas as pd

from common import RES

MIN_N = 20


def h12_rows(df, source, keep):
    df = df[keep(df)]
    out = []
    for _, r in df.iterrows():
        out.append(dict(source=source, context=r.get("ctx", ""), hypothesis=r.hypothesis, vtype=r.vtype,
                        n=int(r.n), estimate=r.all_agree, lo=r.lo, hi=r.hi, p_holm=r.p_holm,
                        verdict="supported" if r.supported_holm else "not supported"))
    return out


def contrast_rows(df, source, hyp, group_cols):
    out = []
    df = df[(df.n_hi >= MIN_N) & (df.n_lo >= MIN_N)]
    for key, g in df.groupby(group_cols, dropna=False):
        key = key if isinstance(key, tuple) else (key,)
        n_pos = int(g.direction_supports.sum())
        n_sig = int(g.sig_supports.sum())
        n_rev = int(((g.diff_hi < 0) & (g.p_holm < 0.05)).sum())
        bins = len(g)
        verdict = ("supported (all bins)" if n_sig == bins else
                   "partial" if n_sig > 0 and n_rev == 0 else
                   "mixed (reversed in some bins)" if n_rev > 0 and n_sig > 0 else
                   "reversed" if n_rev > 0 else "not supported")
        per_bin = "; ".join(f"{b}: {d:+.2f} [{lo:+.2f},{hi:+.2f}]" for b, d, lo, hi in
                            g.sort_values("vaf_bin", key=lambda s: s.map({"<5%": 0, "5-10%": 1, "10-20%": 2, ">=20%": 3}))
                            [["vaf_bin", "diff", "diff_lo", "diff_hi"]].itertuples(index=False))
        out.append(dict(source=source, hypothesis=hyp, context=" | ".join(str(k) for k in key),
                        bins_tested=bins, bins_sig_support=n_sig, bins_sig_reversed=n_rev,
                        verdict=verdict, per_bin=per_bin))
    return out


def main():
    rows12, rows34 = [], []
    # --- D1
    d = pd.read_csv(RES / "d1_h1h2.csv")
    d["ctx"] = d.run + " " + d["flags"] + " " + d.callerset + " " + d.denom
    rows12 += h12_rows(d, "D1 SEQC2 ensemble", lambda x: x.denom.isin(["truth", "union", "union_hc"]) & (x.vtype != "pooled"))
    h3 = pd.read_csv(RES / "d1_h3.csv")
    rows34 += contrast_rows(h3, "D1 SEQC2 ensemble", "H3 indel>SNV", ["flags", "callerset", "denom", "run"])
    h4 = pd.read_csv(RES / "d1_h4.csv")
    rows34 += contrast_rows(h4, "D1 SEQC2 ensemble", "H4 LC>nonLC", ["flags", "callerset", "denom", "vtype", "run"])
    # --- D1b
    d = pd.read_csv(RES / "d1b_h1h2.csv")
    d["ctx"] = "IL_T_1 raw PASS " + d.callerset + " " + d.denom
    rows12 += h12_rows(d, "D1b raw VCFs (6/4 callers)", lambda x: x.vtype != "pooled")
    rows34 += contrast_rows(pd.read_csv(RES / "d1b_h3.csv"), "D1b raw VCFs", "H3 indel>SNV", ["callerset", "denom"])
    rows34 += contrast_rows(pd.read_csv(RES / "d1b_h4.csv"), "D1b raw VCFs", "H4 LC>nonLC", ["callerset", "denom", "vtype"])
    # --- D2 (per purity)
    d = pd.read_csv(RES / "d2_h1h2.csv")
    d["ctx"] = "purity " + d.purity.astype(str) + " " + d.callerset + " " + d.denom
    rows12 += h12_rows(d, "D2 purity series", lambda x: x.n > 0)
    rows34 += contrast_rows(pd.read_csv(RES / "d2_h3.csv"), "D2 purity series", "H3 indel>SNV", ["callerset", "denom", "purity"])
    rows34 += contrast_rows(pd.read_csv(RES / "d2_h4.csv"), "D2 purity series", "H4 LC>nonLC", ["callerset", "denom", "vtype", "purity"])
    # --- D3
    d = pd.read_csv(RES / "d3_h1h2.csv")
    d["ctx"] = "MC3 PASS " + d.centers
    rows12 += h12_rows(d, "D3 TCGA MC3", lambda x: x.n > 0)
    rows34 += contrast_rows(pd.read_csv(RES / "d3_h3.csv"), "D3 TCGA MC3", "H3 indel>SNV", ["centers"])
    rows34 += contrast_rows(pd.read_csv(RES / "d3_h4.csv"), "D3 TCGA MC3", "H4 LC>nonLC", ["centers", "vtype"])

    a = pd.DataFrame(rows12)
    b = pd.DataFrame(rows34)
    a.to_csv(RES / "verdicts_h1h2.csv", index=False)
    b.to_csv(RES / "verdicts_h3h4.csv", index=False)
    # compact tallies
    t12 = a.groupby(["source", "hypothesis", "vtype"]).verdict.apply(
        lambda s: f"{(s == 'supported').sum()}/{len(s)}").unstack("vtype")
    t34 = b.groupby(["source", "hypothesis"]).verdict.value_counts().unstack(fill_value=0)
    with open(RES / "verdicts_tally.md", "w") as f:
        f.write("## H1/H2: configurations where the hypothesis is supported (Holm-adjusted one-sided binomial)\n\n")
        f.write(t12.to_markdown() + "\n\n")
        f.write("## H3/H4: verdict counts across contexts\n\n")
        f.write(t34.to_markdown() + "\n")
    print(t12.to_string())
    print(t34.to_string())


if __name__ == "__main__":
    main()

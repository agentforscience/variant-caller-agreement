"""Figures for the report (static PNGs; every plotted number is also in results/*.csv).

Encoding: SNV = blue, indel = orange (fixed categorical order); non-LC = solid line + circle,
LC = dashed line + triangle (secondary encoding so identity is never colour-only).
"""
from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from common import FIG, RES, VAF_LABELS

CAT = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300"]
TYPE_COL = {"snv": CAT[0], "indel": CAT[1]}
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
plt.rcParams.update({"font.size": 10, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
                     "ytick.color": INK2, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
                     "axes.spines.top": False, "axes.spines.right": False, "lines.linewidth": 2,
                     "figure.facecolor": "#fcfcfb", "axes.facecolor": "#fcfcfb", "legend.frameon": False})
X = np.arange(len(VAF_LABELS))


def _order(df):
    return df.set_index("vaf_bin").reindex(VAF_LABELS)


def ref_lines(ax):
    ax.axhline(0.8, color=INK2, lw=0.8, ls=":")
    ax.axhline(0.5, color=INK2, lw=0.8, ls=":")
    ax.text(len(X) - 0.55, 0.81, "0.80 (H1)", color=INK2, fontsize=8, ha="right", va="bottom")
    ax.text(len(X) - 0.55, 0.51, "0.50 (H2)", color=INK2, fontsize=8, ha="right", va="bottom")


def fig1():
    s = pd.read_csv(RES / "d1_summary.csv")
    s = s[(s.grouping == "vaf_bin") & (s["flags"] == "strict")]
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.2), sharey=True)
    panels = [("truth", "native", "Truth variants (native callers)"),
              ("truth", "shared3", "Truth variants (3 shared callers)"),
              ("union", "native", "Union of all calls (native callers)")]
    for ax, (denom, cs, title) in zip(axes, panels):
        for vt, lab in (("snv", "SNV"), ("indel", "Indel")):
            d = s[(s.denom == denom) & (s.callerset == cs) & (s.vtype == vt)]
            il = _order(d[d.run == "IL_T_1"])
            rng = d.groupby("vaf_bin").all_agree.agg(["min", "max"]).reindex(VAF_LABELS)
            off = -0.06 if vt == "snv" else 0.06
            ax.vlines(X + off, rng["min"], rng["max"], color=TYPE_COL[vt], lw=5, alpha=0.3)
            ax.plot(X + off, il.all_agree, marker="o", ms=8, color=TYPE_COL[vt],
                    label=f"{lab} (IL_T_1; bar = range over 5 runs)")
        ref_lines(ax)
        ax.set_title(title, fontsize=10, color=INK)
        ax.set_xticks(X, VAF_LABELS)
        ax.set_xlabel("VAF bin")
    axes[0].set_ylabel("Fraction called by ALL callers")
    axes[0].set_ylim(-0.02, 1.02)
    axes[0].legend(loc="upper left", fontsize=8)
    fig.suptitle("SEQC2 HCC1395: caller unanimity by VAF (SNV: 5 callers; indel: 3 callers)", color=INK)
    fig.tight_layout()
    fig.savefig(FIG / "fig1_agreement_by_vaf.png", dpi=160)
    plt.close(fig)


def fig2():
    s = pd.read_csv(RES / "d1_summary.csv")
    s = s[(s.grouping == "vaf_bin+lc") & (s["flags"] == "strict") & (s.run == "IL_T_1") & (s.callerset == "native")]
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2), sharey=True)
    for ax, denom, title in ((axes[0], "truth", "Truth variants"), (axes[1], "union", "Union of all calls")):
        for vt, lab in (("snv", "SNV"), ("indel", "Indel")):
            for lc, ls, mk in ((False, "-", "o"), (True, "--", "^")):
                d = _order(s[(s.denom == denom) & (s.vtype == vt) & (s.lc == lc)])
                ok = d.n >= 20  # hide unstable cells
                y = d.all_agree.where(ok)
                ax.plot(X, 1 - y, ls=ls, marker=mk, ms=8, color=TYPE_COL[vt],
                        label=f"{lab}, {'low-complexity' if lc else 'other regions'}")
        ax.set_title(title, fontsize=10)
        ax.set_xticks(X, VAF_LABELS)
        ax.set_xlabel("VAF bin")
    axes[0].set_ylabel("Disagreement (1 - all-callers fraction)")
    axes[0].set_ylim(-0.02, 1.02)
    axes[1].legend(loc="lower left", fontsize=8)
    fig.suptitle("IL_T_1: disagreement by region type (GIAB tandem repeats + homopolymers); cells with n<20 hidden",
                 fontsize=10)
    fig.tight_layout()
    fig.savefig(FIG / "fig2_lowcomplexity_by_vaf.png", dpi=160)
    plt.close(fig)


def fig3():
    p = pd.read_csv(RES / "d2_pooled_expected_vaf.csv")
    pc = pd.read_csv(RES / "d2_paired_clonal.csv")
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    ax = axes[0]
    p = p[(p.callerset == "shared4") & p.lc.notna()]
    p["lc"] = p.lc.astype(str) == "True"
    for vt, lab in (("snv", "SNV"), ("indel", "Indel")):
        for lc, ls, mk in ((False, "-", "o"), (True, "--", "^")):
            d = _order(p[(p.vtype == vt) & (p.lc == lc)])
            y = d.all_agree.where(d.n >= 20)
            ax.plot(X, y, ls=ls, marker=mk, ms=8, color=TYPE_COL[vt],
                    label=f"{lab}, {'LC' if lc else 'non-LC'}")
    ref_lines(ax)
    ax.set_xticks(X, VAF_LABELS)
    ax.set_xlabel("Expected VAF (pooled TVAF x purity)")
    ax.set_ylabel("Fraction called by all 4 shared callers")
    ax.set_ylim(-0.02, 1.02)
    ax.set_title("Truth variants pooled over 6 purities", fontsize=10)
    ax.legend(fontsize=8, loc="upper left")
    ax = axes[1]
    for vt, lab in (("snv", "SNV (5 callers)"), ("indel", "Indel (4 callers)")):
        d = pc[pc.vtype == vt].sort_values("purity")
        first = d.iloc[0]
        xs = list(d.purity) + [1.0]
        ys = list(d.agree_p) + [first.agree_100]
        ax.plot(xs, ys, marker="o", ms=8, color=TYPE_COL[vt], label=f"{lab}, n={int(first.n)}")
    ax.set_xscale("log")
    ax.set_xticks([0.05, 0.1, 0.2, 0.5, 0.75, 1.0], ["5%", "10%", "20%", "50%", "75%", "100%"])
    ax.minorticks_off()
    ax.set_xlabel("Tumour purity (clonal variants with TVAF>=0.4)")
    ax.set_ylabel("Fraction called by all callers")
    ax.set_ylim(-0.02, 1.02)
    ax.set_title("Same clonal variants diluted (paired)", fontsize=10)
    ax.legend(fontsize=8, loc="upper left")
    fig.suptitle("SEQC2 purity series (100x): controlled VAF", fontsize=11)
    fig.tight_layout()
    fig.savefig(FIG / "fig3_purity_series.png", dpi=160)
    plt.close(fig)


def fig4():
    f = RES / "d3_summary.csv"
    if not f.exists():
        return
    s = pd.read_csv(f)
    s = s[(s.grouping == "vaf_bin+lc") & (s["filter"] == "PASS") & (s.centers == "starred_counted")]
    s["lc"] = s.lc.astype(str) == "True"
    fig, ax = plt.subplots(figsize=(6, 4.2))
    for vt, lab in (("snv", "SNV (5 callers)"), ("indel", "Indel (3 callers)")):
        for lc, ls, mk in ((False, "-", "o"), (True, "--", "^")):
            d = _order(s[(s.vtype == vt) & (s.lc == lc)])
            ax.errorbar(X, d.all_agree, yerr=[(d.all_agree - d.all_agree_lo).clip(lower=0), (d.all_agree_hi - d.all_agree).clip(lower=0)],
                        ls=ls, marker=mk, ms=8, color=TYPE_COL[vt], capsize=0,
                        label=f"{lab}, {'LC' if lc else 'non-LC'}")
    ref_lines(ax)
    ax.set_xticks(X, VAF_LABELS)
    ax.set_xlabel("VAF (t_alt_count / t_depth)")
    ax.set_ylabel("Fraction called by all callers")
    ax.set_ylim(-0.02, 1.02)
    ax.set_title("TCGA MC3 (PASS; variants with >=2 callers only)", fontsize=10)
    ax.legend(fontsize=8, loc="lower right")
    fig.tight_layout()
    fig.savefig(FIG / "fig4_mc3.png", dpi=160)
    plt.close(fig)


def fig5():
    s = pd.read_csv(RES / "d1b_summary.csv")
    s = s[(s.grouping == "vaf_bin") & (s.denom == "truth") & (s.callerset == "native")]
    callers = ["MuTect2", "SomaticSniper", "VarDict", "MuSE", "Strelka", "TNscope"]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), sharey=True)
    for ax, vt, title in ((axes[0], "snv", "Truth SNVs"), (axes[1], "indel", "Truth indels")):
        d = _order(s[s.vtype == vt])
        for i, c in enumerate(callers):
            col = f"rate_{c}"
            if col not in d or d[col].isna().all():
                continue
            ax.plot(X, d[col], marker="os^Dvx"[i], ms=8, color=CAT[i], label=c)
        ax.set_xticks(X, VAF_LABELS)
        ax.set_xlabel("VAF bin (pooled TVAF)")
        ax.set_title(title, fontsize=10)
    axes[0].set_ylabel("Per-caller recall (PASS calls)")
    axes[0].legend(fontsize=8, loc="lower right")
    fig.suptitle("IL_T_1 raw VCFs: per-caller sensitivity drives disagreement at low VAF", fontsize=11)
    fig.tight_layout()
    fig.savefig(FIG / "fig5_caller_recall.png", dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    FIG.mkdir(exist_ok=True)
    for f in (fig1, fig2, fig3, fig4, fig5):
        f()
        print("done", f.__name__)

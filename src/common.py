"""Shared helpers: paths, VAF binning, region annotation (bedtools), agreement metrics, CIs.

Every experiment script imports from here so that definitions (VAF bins, "agree",
low-complexity strata) are identical across experiments.
"""
from __future__ import annotations

import math
import os
import subprocess
import tempfile
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "datasets"
RES = ROOT / "results"
FIG = ROOT / "figures"
BEDTOOLS = str(ROOT / "tools" / "bedtools")
SEED = 42

SEQC2 = DATA / "seqc2_hcc1395"
TRUTH_SNV = SEQC2 / "release_v1.2.1" / "high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz"
TRUTH_INDEL = SEQC2 / "release_v1.2.1" / "high-confidence_sINDEL_in_HC_regions_v1.2.1.vcf.gz"
HC_BED = SEQC2 / "release_v1.2.1" / "High-Confidence_Regions_v1.2.bed"

# Main contigs only (both builds: with and without "chr")
MAIN_CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX"]
MAIN_CHROMS_37 = [str(i) for i in range(1, 23)] + ["X"]

# VAF bins used everywhere. Hypothesis thresholds are <5% and >=20%.
VAF_EDGES = [0.0, 0.05, 0.10, 0.20, 1.0000001]
VAF_LABELS = ["<5%", "5-10%", "10-20%", ">=20%"]


def giab_beds(build: str = "GRCh38") -> dict[str, Path]:
    """GIAB v3.3 stratification BEDs used as region annotations."""
    g = DATA / "giab_stratifications" / build
    return {
        "lc": g / "LowComplexity" / f"{build}_AllTandemRepeatsandHomopolymers_slop5.bed.gz",
        "homopol_ge7": g / "LowComplexity" / f"{build}_AllHomopolymers_ge7bp_imperfectge11bp_slop5.bed.gz",
        "tandem_rep": g / "LowComplexity" / f"{build}_AllTandemRepeats.bed.gz",
        "lowmap": g / "Mappability" / f"{build}_lowmappabilityall.bed.gz",
        "segdup": g / "SegmentalDuplications" / f"{build}_segdups.bed.gz",
    }


def vaf_bin(v: pd.Series) -> pd.Categorical:
    return pd.cut(v.clip(0, 1), VAF_EDGES, labels=VAF_LABELS, right=False)


def variant_span(df: pd.DataFrame, chrom="CHROM", pos="POS", ref="REF") -> pd.DataFrame:
    """0-based half-open interval covering the REF allele (1bp for SNVs)."""
    start = df[pos].astype(np.int64) - 1
    end = start + df[ref].str.len().astype(np.int64)
    return pd.DataFrame({"chrom": df[chrom].astype(str).values, "start": start.values, "end": end.values})


def annotate_regions(df: pd.DataFrame, beds: dict[str, Path], chrom="CHROM", pos="POS", ref="REF",
                     span: pd.DataFrame | None = None) -> pd.DataFrame:
    """Add one boolean column per BED: does the variant's REF span overlap the region?

    Uses `bedtools intersect -c` (no sorting needed); rows are tracked by an index column.
    `span` (chrom/start/end, 0-based half-open) may be given to override the REF-based span.
    """
    span = variant_span(df, chrom, pos, ref) if span is None else span[["chrom", "start", "end"]].copy()
    span["idx"] = np.arange(len(df))
    out = df.copy()
    with tempfile.TemporaryDirectory(dir=str(RES)) as td:
        a = os.path.join(td, "a.bed")
        span.to_csv(a, sep="\t", header=False, index=False)
        for name, bed in beds.items():
            res = subprocess.run([BEDTOOLS, "intersect", "-c", "-a", a, "-b", str(bed)],
                                 check=True, capture_output=True, text=True).stdout
            cnt = pd.read_csv(pd.io.common.StringIO(res), sep="\t", header=None,
                              names=["c", "s", "e", "idx", "n"])
            hit = np.zeros(len(df), dtype=bool)
            hit[cnt["idx"].values] = cnt["n"].values > 0
            out[name] = hit
    return out


# ----------------------------------------------------------------------------- statistics
def wilson(k, n, alpha=0.05):
    """Wilson score interval; returns (p, lo, hi). Works on scalars."""
    if n == 0:
        return (np.nan, np.nan, np.nan)
    z = stats.norm.ppf(1 - alpha / 2)
    p = k / n
    den = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / den
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return (p, max(0.0, centre - half), min(1.0, centre + half))


def newcombe_diff(k1, n1, k2, n2, alpha=0.05):
    """Newcombe (hybrid Wilson) CI for p1 - p2."""
    p1, l1, u1 = wilson(k1, n1, alpha)
    p2, l2, u2 = wilson(k2, n2, alpha)
    d = p1 - p2
    lo = d - math.sqrt((p1 - l1) ** 2 + (u2 - p2) ** 2)
    hi = d + math.sqrt((u1 - p1) ** 2 + (p2 - l2) ** 2)
    return d, lo, hi


def fleiss_kappa_binary(M: np.ndarray) -> float:
    """Fleiss' kappa for items x raters binary matrix (each item rated by all K raters)."""
    if len(M) == 0:
        return np.nan
    K = M.shape[1]
    yes = M.sum(1)
    counts = np.stack([K - yes, yes], 1).astype(float)
    P_i = ((counts ** 2).sum(1) - K) / (K * (K - 1))
    p_j = counts.sum(0) / (len(M) * K)
    Pe = (p_j ** 2).sum()
    if Pe >= 1:
        return np.nan
    return float((P_i.mean() - Pe) / (1 - Pe))


def mean_pairwise_jaccard(M: np.ndarray) -> float:
    """Mean Jaccard over caller pairs, computed within the rows given (e.g. a stratum)."""
    vals = []
    for i, j in combinations(range(M.shape[1]), 2):
        a, b = M[:, i].astype(bool), M[:, j].astype(bool)
        u = (a | b).sum()
        if u:
            vals.append((a & b).sum() / u)
    return float(np.mean(vals)) if vals else np.nan


def agreement_summary(df: pd.DataFrame, callers: list[str], group_cols: list[str]) -> pd.DataFrame:
    """Per stratum agreement metrics. `callers` columns must be 0/1.

    all_agree  : called by all K callers
    majority   : called by >= ceil(K/2)+ (strict majority, i.e. > K/2)
    detected_discord : among rows with >=1 caller, fraction NOT called by all
    """
    K = len(callers)
    maj = K // 2 + 1
    rows = []
    for key, g in df.groupby(group_cols, observed=True):
        M = g[callers].to_numpy(dtype=np.int8)
        n_c = M.sum(1)
        N = len(g)
        k_all = int((n_c == K).sum())
        p, lo, hi = wilson(k_all, N)
        det = n_c > 0
        rec = dict(zip(group_cols, key if isinstance(key, tuple) else (key,)))
        rec.update(
            n=N, K=K, n_all_agree=k_all, all_agree=p, all_agree_lo=lo, all_agree_hi=hi,
            majority=float((n_c >= maj).mean()), mean_frac_callers=float(n_c.mean() / K),
            detected=int(det.sum()),
            detected_discord=float((n_c[det] < K).mean()) if det.any() else np.nan,
            singleton=float((n_c == 1).mean()),
            jaccard=mean_pairwise_jaccard(M), fleiss_kappa=fleiss_kappa_binary(M),
        )
        for c in callers:
            rec[f"rate_{c}"] = float(g[c].mean())
        rows.append(rec)
    return pd.DataFrame(rows)


def binom_one_sided(k, n, p0, alternative):
    if n == 0:
        return np.nan
    return float(stats.binomtest(int(k), int(n), p0, alternative=alternative).pvalue)


def holm(pvals):
    p = np.asarray(pvals, dtype=float)
    out = np.full_like(p, np.nan)
    ok = ~np.isnan(p)
    idx = np.argsort(p[ok])
    m = ok.sum()
    adj = np.empty(m)
    running = 0.0
    for r, i in enumerate(idx):
        running = max(running, (m - r) * p[ok][i])
        adj[i] = min(1.0, running)
    out[ok] = adj
    return out


def ensure_dirs():
    for d in [RES, FIG, ROOT / "logs"]:
        d.mkdir(exist_ok=True)

# Do Somatic Variant Callers Agree on the Same Tumor Sample?

This project tests one hypothesis: somatic callers agree on more than 80% of variants with VAF >20%, on fewer than 50% of variants with VAF <5%, and disagree most for indels and in low-complexity (LC) regions regardless of VAF.

We used released per-caller outputs from three sources:
- SEQC2 HCC1395: 5 sequencing runs, a truth set, and a purity-dilution series.
- TCGA MC3: 10k exomes.

Low-complexity regions come from the GIAB v3.3 strata.

Full write-up: **[REPORT.md](REPORT.md)**.

## Key findings

**1. VAF (H1/H2) depends on the denominator.**

| Variant set | Agreement at VAF ≥20% | Agreement at VAF <5% |
|---|---|---|
| True (truth-set) variants, standard-depth runs | 88–96% SNVs, 86–91% indels | ≤0.3% SNVs, 19–37% indels |
| Union of everything the callers call | 29–56% | ≈0 |

- On true variants, H1 and H2 hold.
- On the call union, H1 is refuted: 63–70% of calls are single-caller calls.

**2. VAF is causal.** In the purity-dilution series, the same 16,636 clonal SNVs went from 93% five-caller unanimity at 100% purity to 0.6% at 5% purity.

**3. Low complexity (H4) is supported.** In GIAB tandem-repeat and homopolymer regions:
- disagreement is higher in every VAF bin where agreement is measurable, in all datasets, with no reversed bin;
- the adjusted odds ratio for unanimity is 0.02–0.72;
- the gap varies with VAF, so the effect is not constant.

**4. Indels (H3) are only conditionally worse.**
- Indels disagree more than SNVs at ≥20% VAF everywhere.
- With matched caller sets, indels disagree more in most bins.
- With native caller sets, SNVs look worse at low VAF, because SNV-only SomaticSniper makes almost no true calls below 10% VAF.

**5. Disagreement at low VAF is a sensitivity problem.** It is driven by the weakest caller:
- SNVs: SomaticSniper.
- Indels: TNscope.

## Reproduce

Requirements:
- the `.venv` from `uv`;
- the datasets in `datasets/` (see `datasets/README.md`);
- `tools/bedtools`.

Run:

```bash
uv venv && source .venv/bin/activate && uv sync   # or use the existing .venv
bash run_all.sh                                    # ~45 min, CPU only
```

## Layout

```
src/common.py            shared definitions (VAF bins, bedtools annotation, metrics, Wilson/Newcombe CIs)
src/hypothesis_tests.py  H1-H4 tests, GLMs
src/d1_build.py          Exp 1 build: SEQC2 ensemble tables -> annotated per-site tables
src/d1_analyze.py        Exp 1 metrics/tests/GLMs (5 runs; truth vs union; native vs shared callers; strict vs lenient)
src/d1b_raw.py           Exp 1b: raw IL_T_1 VCFs, PASS-only, TNscope added, leave-one-caller-out
src/d2_purity.py         Exp 2: purity dilution series (controlled VAF)
src/d3_mc3.py            Exp 3: TCGA MC3 replication
src/verdicts.py          consolidated verdict tables
src/make_figures.py      figures/fig1-5
results/                 all CSV outputs (+ per-site parquet files)
figures/                 PNG figures used in REPORT.md
planning.md              pre-registered plan and direction ranking
literature_review.md, resources.md   pre-gathered background
```

# Research State

- Current phase: `None`
- Pipeline completed: `True`

## Previous phases

resource_finder (succeeded), experiment_runner (succeeded)

## Current phase context

- Phase: `experiment_runner`
- Status: `completed`
- Started: `2026-09-17T09:25:53.082936Z`
- Next steps:
  - Validate the report and experimental artifacts before finalizing.

## Workspace check

- Root: `/workspaces/do_variant_callers_agree_on_so_20260917_084425_f6df3f98`
- Directory usable: `True`

## Output validation

- Valid: `True`
- Expected: `REPORT.md`
- Missing: None
- Outside workspace: None

## Agent notes

<!-- NEURICO_AGENT_NOTES_START -->
### resource_finder
<!-- NEURICO_AGENT_NOTES_START:resource_finder -->
**Phase:** resource_finder. COMPLETE (2026-09-17).

**Done**
- 13 PDFs plus 10 Europe PMC full texts.
- Notes and synthesis: `literature_review.md`, with per-paper details in `notes_lit_A.md` and `notes_lit_B.md`.
- 4 datasets, about 10 GB, all passing gzip checks: SEQC2 HCC1395, TCGA MC3, GIAB v3.3 strata, hg38/hg19.
- 6 repositories cloned into `code/`.
- `tools/bedtools` (static v2.31.0).
- The venv has pandas, pysam, pyranges, pyfaidx, scipy, statsmodels, seaborn and pydantic.

**Direction ranking** (details in `planning.md`)
- Kept:
  - D1: SEQC2 same-sample multi-caller agreement, stratified by VAF × type × low complexity.
  - D2: SEQC2 purity-dilution series, a controlled test of VAF.
  - D3: TCGA MC3 pan-cancer replication.
- Rejected:
  - BAMSurgeon simulation, and re-running callers: no toolchain or Docker available.
  - PCAWG: the portal has been retired.
  - DREAM: data are controlled access.
  - WES-aligner analysis and meta-analysis: pruned.

**Key evidence and caveats**
- Literature: Chen 2020 reports Strelka2/Mutect2 concordance above 90% at VAF ≥10%, 75–89% at 5% and 20–41% at 1%.
  - Goode 2013 finds only 46% three-caller agreement even at VAF ≥20% when all calls are counted. The denominator matters, so report both the union of calls and truth-only.
  - Indels are consistently worse.
  - Low-complexity effects have not been tested directly in somatic data. This is the gap.
- SEQC2 ensemble TSVs have usable flags for 5 SNV callers (MuTect2, SomaticSniper, VarDict, MuSE, Strelka) and 3 indel callers (MuTect2, VarDict, Strelka).
  - TNscope, VarScan2 and LoFreq columns are all NaN.
  - A VarDict value of 0.5 is ambiguous.
  - The tables include unfiltered calls.
- SomaticSniper VCFs have no FILTER (use SS=2 plus an SSC threshold).
- Sample column order differs between caller VCFs.
- The MC3 public MAF contains only variants with NCALLERS ≥2. SNV calls have 5 possible callers and indel calls have 3, so normalise.
- Pilot results (a sanity check only, not final):
  - SEQC2 IL_T_1 truth SNVs, mean callers out of 5: 2.22 at VAF ≤5%, 3.65 at 5–10%, 4.68 at 10–20%, 4.96 at >20%.
  - MC3 PASS (first 300k rows), fraction of SNVs with ≥4 callers: 0.79 outside low complexity vs 0.72 inside.

**Next phase: experiment_runner**
1. Implement D1 using `datasets/README.md`: load the code, join the truth set and GIAB BEDs, compute agreement metrics with Wilson CIs, and fit the GLM interaction tests.
2. D2: parse the purity VCFs with caller-specific PASS/VAF rules.
3. D3: MC3 normalised agreement by VAF and strata.
4. Use the hypothesis operationalisation in `planning.md`.

**Unresolved**
- The SEQC2 Nat Biotech papers (Fang/Xiao 2021), MC3 (Ellrott 2018) and Strelka2 PDFs could not be downloaded; we have abstract-level knowledge only.
- The hg19 vs GRCh37 chrM difference is irrelevant for this work.
<!-- NEURICO_AGENT_NOTES_END:resource_finder -->

### experiment_runner
<!-- NEURICO_AGENT_NOTES_START:experiment_runner -->
**Phase:** experiment_runner. COMPLETE (2026-09-17). All phases are done. REPORT.md, README.md, results/ and figures/ exist.

**What was done**
- `planning.md` was rewritten. It adds a Motivation & Novelty section and a hypothesis operationalisation, and keeps the direction ranking. The original is saved at `logs/planning_resource_finder_orig.md`.
- Top 3 directions implemented:
  - D1: SEQC2 ensemble tables, 5 runs (`src/d1_*.py`), plus robustness check D1b on the raw IL_T_1 VCFs with TNscope added and leave-one-caller-out (`src/d1b_raw.py`).
  - D2: purity series (`src/d2_purity.py`).
  - D3: MC3 (`src/d3_mc3.py`).
  - Verdicts: `src/verdicts.py`. Figures: `src/make_figures.py`. Full pipeline: `run_all.sh`.

**Key findings** (details and evidence in REPORT.md; tables in `results/verdicts_*`)
- **H1/H2:**
  - Supported on the truth denominator in the standard-depth runs: ≥20% VAF gives 0.88–0.96 for SNVs and 0.87–0.91 for indels; <5% gives ≤0.003 for SNVs and 0.19–0.37 for indels.
  - H1 is refuted on the call union (0.29–0.56).
  - The purity series shows VAF is causal.
  - NS_380X fails H1, because MuTect2 and MuSE lose high-VAF calls at 380×.
- **H4:** LC is worse in every measurable bin in all datasets, with no reversals. The effect size varies with VAF (interaction significant). The "partial" verdicts are floor effects.
- **H3:** Indels are worse at ≥20% VAF everywhere, and in most bins with matched caller sets. With native caller sets the ranking reverses at low VAF, because SomaticSniper and TNscope are weak there.

**Caveats**
- The truth set is biased toward agreement, and truth LC cells are small.
- The MC3 public MAF is truncated to variants with ≥2 callers, and the meaning of the `*` suffix is uncertain.
- Indels were matched exactly after vt-style normalisation.

**Validation**
- The D3 rerun was byte-identical (`logs/d3_mc3_rerun.log`).
- The report numbers were cross-checked against the CSVs.

**Next steps (optional):** paper-writer; finer repeat strata; haplotype-aware matching; modern callers.
<!-- NEURICO_AGENT_NOTES_END:experiment_runner -->

<!-- NEURICO_AGENT_NOTES_END -->

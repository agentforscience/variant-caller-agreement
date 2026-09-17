# Research Plan: Do somatic variant callers agree on the same tumor sample?

## Motivation & Novelty Assessment

### Why This Research Matters
Clinical and research pipelines report somatic mutations from one caller, or from a consensus of several. Two things matter to those users:
- **Where** callers disagree, which sets the reliability of a reported mutation.
- **Whether a consensus loses real variants.** Subclonal (low-VAF) variants and indels in repetitive DNA are the drivers of therapy resistance and of microsatellite instability (MSI) phenotypes, and they are the variants most at risk.

Quantifying agreement jointly by VAF, variant type and sequence context tells pipeline builders when single-caller output can be trusted and when orthogonal evidence is needed.

### Gap in Existing Work
From `literature_review.md`:
- **VAF.** Only Chen 2020 reports VAF-stratified concordance, and it covers 2 callers and SNVs only.
- **Indels.** Many papers show indels are harder (Alioto, Krøigård, Guille), but none tests whether the indel penalty persists **within** each VAF stratum.
- **Low complexity.** The effect is inferred from germline GIAB work and from the exclusion of such regions in somatic truth sets (Jones 2021). It has not been measured for multi-caller somatic agreement.
- **Denominator.** Published numbers mix denominators: "all calls" in Goode (46% at VAF ≥20%) versus "validated variants" in Wang 2013 (82%). This makes the 80% / 50% thresholds ambiguous.

### Our Novel Contribution
We run a joint VAF × variant-type × low-complexity analysis of agreement among 5 SNV callers and 3 indel callers on the same tumor (SEQC2 HCC1395).
- **Two denominators:** the union of calls, and the truth-set variants.
- **Replication:** across 5 sequencing runs.
- **Controlled-VAF test:** the purity-dilution series. The same variants are diluted, which gives a within-variant test of the VAF effect.
- **Pan-cancer replication:** 10k TCGA exomes (MC3).

We test the "regardless of VAF" claims explicitly, using within-bin contrasts and interaction terms.

### Experiment Justification
- **Exp 1 (D1, SEQC2 ensemble tables, 5 runs).** The core test of H1–H4 on one tumor with a truth set. It is needed because only here do we have ≥5 callers on identical reads, plus truth labels that separate sensitivity-driven from false-positive-driven disagreement.
- **Exp 1b (robustness, raw IL_T_1 VCFs).** A PASS-only 6-caller set that adds TNscope. It checks that conclusions do not depend on SomaticSeq's call-flag definitions.
- **Exp 2 (D2, purity series).** Observed VAF is confounded with site properties such as copy number, mappability and context. Diluting the same variants changes VAF alone, so it is a causal test of H1/H2 and of whether the indel and low-complexity (LC) penalties hold at every VAF.
- **Exp 3 (D3, TCGA MC3).** Checks generalisation to other tumors, callers (MuTect1, RADIA, VarScan2, Pindel, Indelocator) and exome data. It has no truth set, and only calls made by ≥2 callers are public.

## Research Question
For the same tumor sample, what fraction of somatic variants do different callers agree on, and how does this depend on:
- VAF (<5%, 5–10%, 10–20%, ≥20%);
- variant type (SNV vs indel);
- low-complexity context (GIAB tandem-repeat/homopolymer regions)?

## Hypothesis Decomposition
"Agree" (primary) means the variant is called by **all** callers in the caller set (5 for SNVs, 3 for indels). Secondary definitions: called by a majority of callers, and normalised mean callers (n/K).
- **H1:** all-agree fraction > 0.80 for VAF ≥ 0.20.
- **H2:** all-agree fraction < 0.50 for VAF < 0.05.
- **H3:** disagreement(indel) > disagreement(SNV) in **every** VAF bin. Tested two ways: with the native caller sets, and like-for-like with the 3 shared callers (MuTect2, VarDict, Strelka).
- **H4:** disagreement(LC) > disagreement(non-LC) in every VAF bin, for SNVs and for indels separately.
- **"Highest" claim:** the indel-in-LC cell has the highest disagreement within each VAF bin.

Variables:
- Independent: VAF bin, variant type, LC status.
- Dependent: all-agree indicator.
- Covariates: depth, low mappability, segmental duplications.

Each hypothesis is evaluated on two denominators:
- **U (union of calls):** every site called by ≥1 caller.
- **T (truth):** SEQC2 v1.2.1 high-confidence variants. The VAF used is the pooled TVAF, which is caller-independent and not conditioned on detection.

## Proposed Methodology

### Approach
We use released per-caller outputs, because the tools needed to re-run callers are not available (see the rejected directions below).

"Called" definitions (primary = PASS-like):
- **MuTect2:** `if_MuTect==1`.
- **SomaticSniper:** `SS==2` and `SSC≥40`.
- **VarDict:** `==1` (PASS & Somatic).
- **MuSE:** `MuSE_Tier==1` (PASS).
- **Strelka:** `==1` (PASS).

Lenient variant, used as a sensitivity analysis: VarDict ≥0.5, any MuSE tier, and any SomaticSniper SS=2.

Low complexity (LC):
- **Primary:** GIAB v3.3 `AllTandemRepeatsandHomopolymers_slop5`.
- **Secondary:** homopolymer-length strata, and SEQC2 LC score.

VAF bins: [0, 0.05), [0.05, 0.10), [0.10, 0.20), [0.20, 1].

Scope: chromosomes 1–22 and X only.

### Experimental Steps
1. **Annotate the ensemble tables.** For each of the 5 SEQC2 runs (IL, NV, FD, NS WGS, and the 380× combined sample), load the SNV and indel tables and derive caller flags and observed VAF. Annotate with truth membership, HC region, and GIAB strata (LC, homopolymer, low mappability, segdup) using bedtools.
2. **Build the truth table.** Add truth variants absent from the table (0 callers). Use TVAF, and match on CHROM/POS/REF/ALT.
3. **Stratified metrics:**
   - all-agree, majority and n/K;
   - pairwise Jaccard and Fleiss' κ;
   - per-caller recall on T and positive fraction on U;
   - Wilson 95% CIs.
4. **Hypothesis tests:**
   - One-sided exact binomial tests against 0.80 (H1) and 0.50 (H2).
   - Within-bin risk differences for H3/H4 with Newcombe CIs, and χ² tests with Holm correction.
   - Logistic GLMs: `all_agree ~ C(vaf_bin)*is_indel + C(vaf_bin)*lc + log(depth) + lowmap + segdup`, fitted separately for T and U. Interaction likelihood-ratio tests (LRTs) are used to test "regardless of VAF".
5. **Exp 1b.** Parse the raw IL_T_1 VCFs (PASS only, SomaticSniper SSC≥40) and repeat the key tables with 6 SNV callers and 4 indel callers (MuTect2, VarDict, Strelka, TNscope).
6. **Exp 2.** For each of 6 purities:
   - Parse the PASS calls from 5 callers: TNscope, Lancet, MuTect2, SomaticSniper (SNVs only), and Strelka. The indel caller set is therefore 4.
   - Match calls to the truth variants and set expected VAF = TVAF × purity.
   - Compute agreement by expected-VAF bin × type × LC. For the union analysis, use the median caller-reported VAF.
7. **Exp 3.** From MC3, keep PASS rows.
   - Callers: SNV = {MUTECT, MUSE, VARSCANS, SOMATICSNIPER, RADIA}; indel = {INDELOCATOR, VARSCANI, PINDEL}.
   - Per sample, K = the callers active in that sample, i.e. that contributed ≥1 call of that type.
   - All-agree = n == K.
   - Annotate with GRCh37 GIAB LC; VAF = t_alt/t_depth.
   - Starred callers are counted as calls, consistent with NCALLERS. A strict version that excludes them is a sensitivity analysis.

### Baselines / reference points
- **Literature values:** Chen 2020 concordance (>90% at ≥10%, 75–89% at 5%); Goode 2013 (46% at ≥20%, all calls); Wang 2013 (82% of validated SNVs).
- **Internal baselines:**
  - the lenient vs strict call definitions;
  - SNV vs indel, and non-LC vs LC, as the within-study controls;
  - the like-for-like 3-caller comparison.

### Evaluation Metrics
- **All-agree fraction.** The primary metric, and the direct operationalisation of the hypothesis.
- **Majority fraction and n/K.** Less sensitive to a single outlier caller.
- **Pairwise Jaccard.** Standard for set overlap.
- **Fleiss' κ.** Chance-corrected agreement.
- **Per-caller recall.** Explains the disagreement.

### Statistical Analysis Plan
- α = 0.05, with Holm correction within each hypothesis family.
- Wilson CIs for proportions; Newcombe CIs for differences.
- GLM odds ratios with 95% CIs.
- Replicate consistency: a claim of "in every bin" requires the same sign in all 5 runs.
- Large n makes tiny effects significant, so effect sizes (risk differences) are emphasised.

## Expected Outcomes
- **Truth denominator:** H1 is likely supported for SNVs and possibly not for indels. H2 is likely supported.
- **Union denominator:** H1 is likely *refuted*, because private false positives dominate.
- **H3:** likely supported.
- **H4:** possibly supported for indels but weaker for SNVs; the effect may vary with VAF.

Any of these can be refuted, and all outcomes will be reported.

## Timeline and Milestones
| Step | Time |
|---|---|
| Exp 1 | 1 h |
| Exp 1b | 30 min |
| Exp 2 | 45 min |
| Exp 3 | 30 min |
| Stats and figures | 45 min |
| Report | 45 min |

The plan includes a debugging buffer of about 25%.

## Potential Challenges
- **Indel representation mismatches** between callers or truth. SomaticSeq normalised the tables, but the purity-series VCFs were not normalised. We apply left-normalisation with pysam and the reference, then exact matching.
- **Truth-set bias toward agreement.** The truth set was built from callers and excludes some difficult regions. The union denominator reduces this bias; it is noted in the report as a caveat.
- **Few truth indels in LC regions / at low VAF.** We report CIs and n per cell, and do not over-interpret small cells.
- **Large files** (Lancet 7M rows, VarDict 550 MB). We stream them and pre-filter PASS rows with zcat/awk.

## Success Criteria
- Each of H1–H4 receives a supported, refuted or partially-supported verdict, with CIs, on both denominators.
- The verdicts are replicated across runs, the purity series and MC3.
- All results are saved in `results/` and figures in `figures/`.

---

## Direction Ranking (from resource_finder phase; retained, top 3 kept)

| # | Direction | Lit | Rel | Gain | Feas | Total | Decision |
|---|-----------|-----|-----|------|------|-------|----------|
| D1 | SEQC2 same-sample multi-caller agreement, stratified | 5 | 5 | 5 | 5 | 20 | **KEEP** (Exp 1, 1b) |
| D2 | SEQC2 purity-dilution series | 4 | 5 | 5 | 4 | 18 | **KEEP** (Exp 2) |
| D3 | TCGA MC3 pan-cancer replication | 4 | 4 | 4 | 5 | 17 | **KEEP** (Exp 3) |
| D4 | BAMSurgeon spike-ins and re-calling | 4 | 5 | 3 | 1 | 13 | Rejected: no aligner, caller or Docker toolchain |
| D5 | Re-run modern callers on SEQC2 BAMs | 3 | 4 | 3 | 1 | 11 | Rejected: same toolchain barrier |
| D6 | PCAWG consensus | 3 | 4 | 4 | 1 | 12 | Rejected: portal retired / controlled access |
| D7 | DREAM SMC | 4 | 3 | 3 | 1 | 11 | Rejected: controlled access, no per-caller calls |
| D8 | SEQC2 WES aligner × caller | 3 | 3 | 3 | 5 | 14 | Pruned: outside the top 3 |
| D9 | Meta-analysis of published numbers | 3 | 3 | 2 | 4 | 12 | Pruned: used as context only |

The original resource-finder version of this file is kept at `logs/planning_resource_finder_orig.md`.

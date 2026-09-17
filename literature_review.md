# Literature Review: Do somatic variant callers agree on the same tumor sample?

**Hypothesis under test.**
- H1: Callers agree on more than 80% of variants with VAF above 20%.
- H2: Callers agree on fewer than 50% of variants with VAF below 5%.
- H3: Indels disagree more than SNVs, regardless of VAF.
- H4: Variants in low-complexity regions disagree more, regardless of VAF.

Detailed per-paper notes with page references are in `notes_lit_A.md` (8 papers) and `notes_lit_B.md` (16 sources). This file is the synthesis.

## 1. Research area overview

Somatic small-variant calling compares tumor and matched-normal reads to find mutations that are present only in the tumor. The callers use different statistical models:
- Bayesian / likelihood: MuTect, Mutect2, MuSE, Strelka2.
- Heuristic / threshold: VarScan2, VarDict, SomaticSniper.
- Local assembly: Mutect2, Lancet, TNscope, Scalpel.
- Deep learning: NeuSomatic, DeepSomatic.

Benchmarks since 2013 have repeatedly found poor raw overlap between callers. That finding motivated consensus and ensemble methods (Goode 2013; SomaticSeq, Fang 2015; SMuRF, Huang 2019; SomaticCombiner, Wang 2020) and community truth sets:
- ICGC benchmark (Alioto 2015).
- DREAM SMC synthetic tumors (Ewing 2015).
- TCGA MC3 multi-caller calls (Ellrott 2018).
- SEQC2 HCC1395 reference (Fang 2021; Xiao 2021; Zhao 2021).

Most papers report sensitivity and precision against a truth set. Few report **inter-caller agreement stratified by VAF, variant type and genomic context at the same time**. That joint stratification is the gap this project targets.

## 2. Key papers (condensed)

| Paper | Data / callers | Key hypothesis-relevant result |
|---|---|---|
| **Chen et al. 2020**, Sci Rep | WES in-silico mixtures (1–40% VAF, 100–800×). Strelka2 vs Mutect2. SNVs only. | Concordance **>90% except at low VAF**: **75–89% at 5% VAF**, **20–41% at 1% VAF**. Recall 92–97% at ≥20% VAF vs 2.7–35% at 1%. This is the only direct VAF-stratified concordance. |
| **Alioto et al. 2015**, Nat Commun | ICGC medulloblastoma WGS, 18 SNV and 16 indel pipelines, curated Gold Set. | Only **205 SNVs and 1 indel** were called by all pipelines. Best F1 was 0.79 for SNVs and 0.65 for indels. **83% of true indels lie in tandem repeats (71% in homopolymers)** and are under-called. False positives are enriched at low AF. |
| **Krøigård et al. 2016**, PLoS One | 5 breast tumors, WES 80× plus deep panel 362×, 9 callers. | 472 SNVs were called by all 8 SNV callers vs 22,032 by only one. For indels, 27 were called by all 5 callers vs 10,541 by one. Indel sensitivity was 45–74% vs up to 96% for SNVs. Agreement rises with depth. |
| **Goode et al. 2013**, Genome Med | 3 callers on exomes. | Three-caller agreement was **11.5%** of calls. It rose to only **46.3%** after restricting to tumor VAF ≥20%. Calls made by all 3 callers validated 98.9% of the time, vs 5.7% for singletons. |
| **Wang et al. 2013**, Genome Med | 5 callers, WGS melanoma + WES lung. | 82% (WGS) and 64% (WES) of *validated* SNVs were called by all 5. |
| **Cai et al. 2016**, Sci Rep | 57 gallbladder tumors, WES + 370× panel, 4 callers. | Only ~21% of SNVs were called by ≥2 callers. Consensus calls cluster at mid VAF. |
| **Fang et al. 2015** (SomaticSeq), Genome Biol | DREAM synthetic data. | Five-caller agreement was **3–18%** of raw calls. |
| **de Schaetzen van Brienen 2020**, BMC Med Genomics | FF vs FFPE WGS, 4 callers. | Only 8 calls were shared by all 4 callers. A ~0.25 VAF clone was reproducible; a ~0.05 VAF subclone was mostly lost. F1 was 0.82–0.95 for clonal variants vs 0–0.72 for subclonal. |
| **Cibulskis et al. 2013** (MuTect), Nat Biotechnol | Virtual tumors. | At 5% VAF and 30×, sensitivity was MuTect 16%, Strelka 4.6%, others ≤2%. At AF 0.4, all callers exceeded 91%. |
| **Sergi et al. 2024**, BMC Bioinf | Artificial low-tumor-fraction samples. | Below 5% VAF, TPR ranged 0.22–0.70. **Only Mutect2 called any low-VAF indels** at default settings (sensitivity 0.198). |
| **Guille et al. 2025**, Brief Bioinform | 20 callers on WES, plus voting ensembles. | 11/18 callers reached F1 >0.8 for SNVs, vs **2/15** for indels. |
| **Huang et al. 2019** (SMuRF), Bioinformatics | Ensemble of callers. | Individual callers had indel precision <8%. Three-caller indel consensus reached F1 0.46. |
| **Wang et al. 2020** (SomaticCombiner), Sci Rep | 8 callers, mixture data. | VarScan and SomaticSniper collapse at ≤10% VAF. Indel false positives are "more severe". |
| **Jones et al. 2021**, Genome Biol | SEQC2 pooled cell-line reference, 22 pipelines. | Low-complexity regions were **excluded** as irreproducible. Outside the high-confidence regions, spurious calls were >20% vs <1% inside. SNVs are easier than indels "regardless of the VAF". |
| **Ewing et al. 2015** (DREAM), Nat Methods | BAMSurgeon synthetic WGS, 248 submissions. | Only 31% of TCGA lung SNVs were called by all 4 centers (cited result). False-positive artifacts include contexts that create homopolymers. Mapping quality drives false negatives. |
| **Bian et al. 2018**, BMC Bioinf | DREAM sets, several callers, 10 Mb segmentation. | Sensitivity spread widened from 93–98% (clonal) to 54–77% (subclonal). No GC or gene-density effect at 10 Mb scale. Repeat strata were not tested. |
| **Fan et al. 2016** (MuSE), Genome Biol | TCGA WES. | Shared calls validated at 92% vs ~25% for private calls. Validation rate drops in CpG islands. |
| **Benjamin et al. 2019** (Mutect2), bioRxiv | DREAM, mixtures, MC3. | Indel sensitivity is 10–30 points below SNV sensitivity at matched AF. Mutect2 has a dedicated STR/repeat filter. |
| **Dwarshuis et al. 2024** (GIAB strata), Nat Commun | Germline benchmarks. | Tandem repeats and homopolymers show below-average precision and recall. Errors increase more than 10-fold in tandem repeats with more than 10 variants. Defines the BEDs used here. |
| **Zhao 2021 / Fang 2021 / Xiao 2021** (SEQC2) | HCC1395 WGS/WES from 6 centers, 3 aligners, 5–7 callers, purity/depth series. | The truth set reaches down to ~5% VAF. Biological replicates matter more than bioinformatic ones. Depth improves low-VAF reproducibility. (Fang and Xiao were read from abstracts and secondary sources only.) |
| **Ellrott et al. 2018** (MC3), Cell Systems | 10,295 TCGA exomes, 7 callers. | Public MAF of calls made by ≥2 callers, with per-variant `CENTERS`. (Read from abstract and secondary sources only.) |

## 3. Synthesis against the hypothesis

- **H1 (>80% at VAF >20%): conditionally supported. It depends on the denominator.**
  - *True* high-VAF SNVs: callers have 91–99% sensitivity (Cibulskis, Chen), so agreement should be at least 80%.
  - *All calls*, including each caller's private false positives: agreement is far lower. Goode found 46% even at VAF ≥20%.
  - The experiment must report both denominators.
- **H2 (<50% at VAF <5%): indirectly supported.**
  - Sensitivities at ≤5% VAF vary several-fold across callers (16% vs 4.6%; 0.22–0.70), which forces agreement below 50%.
  - Chen shows 75–89% agreement at exactly 5% VAF (2 similar callers, high depth), so the threshold depends on depth and on the caller set.
- **H3 (indels worse regardless of VAF): well supported overall.**
  - Alioto, Krøigård, Guille, SMuRF and Mutect2 all show it.
  - No paper tests the indel × VAF interaction directly.
- **H4 (low-complexity worse regardless of VAF): plausible, but not directly tested in somatic data.**
  - Evidence comes from germline work (GIAB), from indel-in-repeat under-calling (Alioto), and from the exclusion of these regions in Jones.
  - Confound: most indels sit in repeats (83%). Variant type and region effects must therefore be separated, e.g. SNVs inside vs outside low-complexity regions, and indels inside vs outside.
  - Truth sets usually exclude difficult regions, so this stratum is under-sampled.

## 4. Common methodologies

- Intersection counts, Venn/UpSet plots, and "called by k of n callers" histograms (Krøigård, Alioto, Cai, Goode).
- Precision, recall and F1 against a truth set, computed with som.py/hap.py or vcfeval (Chen, Alioto, Ewing).
- Mixture or dilution designs to control VAF: BAM mixing (Chen, Xu 2014), DNA dilution (SEQC2 purity series), and spike-ins (BAMSurgeon; Ewing, Sergi).
- Region stratification with GIAB BEDs (Dwarshuis), or with feature-importance models of FP/FN (Ewing, Alioto).

## 5. Standard baselines (callers)

- **SNV callers:** MuTect/Mutect2, Strelka2, SomaticSniper, VarScan2, VarDict, MuSE, RADIA, LoFreq, TNscope, Lancet, NeuSomatic.
- **Indel callers:** Mutect2, Strelka2, VarDict, Pindel, Indelocator, VarScan2-indel, Scalpel, Lancet, TNscope.
- **Consensus baselines:** majority vote (≥k of n), SomaticSeq, SMuRF.

## 6. Evaluation metrics recommended for this project

1. **Agreement on the union of calls.** Fraction of the call union supported by all n callers, and by at least ⌈n/2⌉ callers. Also mean pairwise Jaccard.
2. **Agreement on true variants.** Fraction of truth variants detected by all callers, and mean number of detecting callers divided by n. This separates sensitivity-driven disagreement from false-positive-driven disagreement.
3. **Chance-corrected agreement.** Fleiss' κ or Krippendorff's α across callers, per stratum.
4. **Per-caller recall and precision** per stratum. These give the context needed to interpret the agreement numbers.
5. **Statistics.**
   - Wilson 95% CIs per bin.
   - Logistic regression or GLM: `all_agree ~ VAF_bin * variant_type + low_complexity * VAF_bin + depth`. A non-significant interaction term supports "regardless of VAF".
   - χ² / Cochran–Armitage trend test across VAF bins.
6. **VAF assignment.** Use one caller-independent VAF per variant: SEQC2 `TVAF` from the pooled replicates, or tumor read counts from the SomaticSeq TSV, or MAF `t_alt_count/t_depth`. **Do not use each caller's own AF field.**
7. **VAF bins.**
   - Minimum set: <5%, 5–10%, 10–20%, ≥20%, with ≥20% optionally split into 20–50% and >50%.
   - Finer set: <1%, 1–5%, 5–10%, 10–20%, 20–50%, >50%.

## 7. Datasets in the literature, and what we use

- **SEQC2 HCC1395** (used; primary).
  - The same tumor called by up to 7 callers.
  - Truth set v1.2.1 with 39,447 SNVs and 1,625 indels.
  - Pooled VAF and linguistic-complexity (LC) annotations in the superSet.
  - Purity series from 100% to 5% tumor.
- **TCGA MC3** (used; secondary).
  - Pan-cancer, per-variant list of calling callers.
  - The public version is restricted to calls made by ≥2 callers.
- **Not used:** DREAM SMC and ICGC benchmark data (controlled/large); PCAWG consensus (portal retired); Jones 2021 Sample A (exome only, low-complexity regions excluded, mostly germline variants).

## 8. Gaps and opportunities

- No published analysis jointly stratifies multi-caller (≥5) agreement by VAF × variant type × low-complexity context on the same tumor.
- Consensus-derived truth sets are biased toward agreement. Using the superSet or ensemble *union*, rather than only the truth set, reduces that bias.
- The purity dilution series allows a within-variant causal test: the same variants observed at lower VAF. Published caller comparisons have not used SEQC2 this way.

## 9. Recommendations for the experiment

- **Primary dataset.** SEQC2 SomaticSeq ensemble TSVs: 4 WGS replicates plus the 380× combined sample.
  - Caller sets: 5 SNV callers (MuTect2, SomaticSniper, VarDict, MuSE, Strelka) and 3 indel callers (MuTect2, VarDict, Strelka).
  - Also run a PASS-only sensitivity analysis on the IL_T_1 raw VCFs, adding TNscope.
- **Controlled VAF test.** SEQC2 purity series at 100×, with 5 callers at 6 purities. Track truth variants whose expected VAF is `TVAF × purity`.
- **Generalisation.** TCGA MC3, restricted to PASS calls. Normalise agreement by the number of callers applicable to each variant type, and state the ≥2-caller truncation.
- **Low-complexity definition.**
  - Primary: GIAB `AllTandemRepeatsandHomopolymers_slop5`.
  - Secondary: the homopolymer-length strata, SEQC2 `LC` (low value = low complexity), and `SiteHomopolymer_Length` from the TSVs.
  - Control for mappability and segmental duplications using the GIAB BEDs.
- **Pitfalls.**
  - VarDict's 0.5 flag is ambiguous. Run the analysis both ways.
  - SomaticSniper has no PASS filter. Use SS=2 plus an SSC threshold.
  - Indel representation differs between callers. Normalise, or match within a ±5 bp window.
  - The number of callers differs between SNVs and indels, so compare normalised agreement, or compare only the callers that call both (MuTect2, VarDict, Strelka).

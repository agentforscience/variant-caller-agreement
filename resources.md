# Resources Catalog

## Summary
Resources gathered to test whether somatic variant callers agree on the same tumor sample, broken down by VAF, variant type (SNV vs indel), and low-complexity context:
- 13 PDFs plus 10 open-access full texts
- 4 datasets (about 10 GB)
- 6 cloned repositories
- a static bedtools binary

Direction ranking and pruning are in `planning.md`.

## Papers
In total: 13 PDFs in `papers/` and 10 full texts in `papers/fulltext_xml/`. The table lists the most important ones; see `papers/README.md` for the complete list.

| Title (short) | Authors | Year | File | Key info |
|---|---|---|---|---|
| Somatic calling vs depth and mutation frequency | Chen et al. | 2020 | papers/Chen2020_depth_mutation_frequency.pdf | Only paper that reports caller concordance by VAF: >90% at ≥10%, 75–89% at 5%, 20–41% at 1% |
| ICGC WGS somatic benchmark | Alioto et al. | 2015 | papers/Alioto2015_ICGC_somatic_benchmark.pdf | 205 SNVs vs 1 indel called by all pipelines; 83% of true indels are in tandem repeats |
| Nine somatic variant callers | Krøigård et al. | 2016 | papers/Kroigard2016_nine_somatic_callers.pdf | Most calls come from a single caller; indels agree less |
| Callers vs sequencing depth | Cai et al. | 2016 | papers/Cai2016_callers_sequencing_depth.pdf | About 21% of calls made by ≥2 callers |
| DREAM SMC challenge | Ewing et al. | 2015 | papers/Ewing2015_DREAM_SMC.pdf | Synthetic truth; ensembles |
| FF vs FFPE somatic calling | de Schaetzen van Brienen et al. | 2020 | papers/deSchaetzenvanBrienen2020_FF_vs_FFPE_somatic.pdf | 8 calls shared by all 4 callers; subclones lost |
| Low-VAF verified reference | Jones et al. | 2021 | papers/Jones2021_verified_reference_low_VAF_panels.pdf | Low-complexity regions excluded because they were irreproducible |
| SEQC2 WGS/WES data | Zhao et al. | 2021 | papers/Zhao2021_SEQC2_WGS_WES_data.pdf | Describes the primary dataset |
| SomaticCombiner | Wang et al. | 2020 | papers/Wang2020_SomaticCombiner.pdf | Consensus; low-VAF collapse |
| NeuSomatic | Sahraeian et al. | 2019 | papers/Sahraeian2019_NeuSomatic.pdf | Indels are harder than SNVs |
| Mutect2 | Benjamin et al. | 2019 | papers/Benjamin2019_Mutect2_somatic_SNVs_indels.pdf | Sensitivity by AF and depth; STR filter |
| MuTect | Cibulskis et al. | 2013 | papers/Cibulskis2013_MuTect.pdf | At 5% VAF, sensitivities of 16% vs 4.6% |
| GIAB stratifications | Dwarshuis et al. | 2024 | papers/Dwarshuis2024_GIAB_stratifications.pdf | Definitions of the region strata |
| Consensus approach | Goode et al. | 2013 | papers/fulltext_xml/Goode2013_consensus_somatic.txt | 3-caller agreement is 11.5%, and 46.3% at VAF ≥20% |
| Comparison of mutation callers | Wang et al. | 2013 | papers/fulltext_xml/Wang2013_comparison_mutation_callers.txt | 82% / 64% of validated SNVs called by all 5 callers |
| SomaticSeq | Fang et al. | 2015 | papers/fulltext_xml/Fang2015_SomaticSeq.txt | 5-way agreement is 3–18% |
| 20 callers + ensembles on WES | Guille et al. | 2025 | papers/fulltext_xml/Guille2025_callers_voting_ensembles_WES.txt | F1 >0.8 for 11/18 SNV callers vs 2/15 indel callers |
| Low-tumor-fraction benchmark | Sergi et al. | 2024 | papers/fulltext_xml/Sergi2024_artificial_low_tumor_fraction.txt | True-positive rate (TPR) 0.22–0.70 below 5% VAF |

## Datasets
Total: 4. See `datasets/README.md` for details, download commands, and caveats.

| Name | Source | Size | Task | Location | Notes |
|---|---|---|---|---|---|
| SEQC2 HCC1395 multi-caller calls | NCBI FTP (SEQC2) | 4.7 GB | Agreement between callers on the same sample; truth set | datasets/seqc2_hcc1395/ | Truth v1.2.1 (39,447 SNVs, 1,625 indels). SomaticSeq ensemble TSVs for 5 replicates. Raw VCFs from 7 callers. Purity series with 30 VCFs. 108 WES VCFs. |
| TCGA MC3 public MAF | GDC | 753 MB (gz) | Pan-cancer agreement between callers | datasets/tcga_mc3/ | 3.6M variants from 10,295 tumors. Caller lists in `CENTERS` and `NCALLERS`. Only variants called by ≥2 callers are included. |
| GIAB stratifications v3.3 | NIST FTP | 867 MB | Low-complexity, mappability and segmental-duplication strata | datasets/giab_stratifications/ | GRCh38 + GRCh37 |
| UCSC hg38 / hg19 FASTA | UCSC | 6.3 GB | Sequence context | datasets/reference/ | Indexed (.fai) |

## Code repositories
Total: 6. See `code/README.md`.

| Name | URL | Purpose | Location | Notes |
|---|---|---|---|---|
| somaticseq | github.com/bioinform/somaticseq | Generated the SEQC2 ensemble files; VCF parsers for each caller; linguistic-complexity score | code/somaticseq/ | LC function tested and working (needs pydantic, which is installed) |
| hap.py / som.py | github.com/Illumina/hap.py | Stratified somatic benchmarking | code/hap.py/ | Not built (needs C++ or Docker); reference only |
| bamsurgeon | github.com/adamewing/bamsurgeon | Spike-in simulation | code/bamsurgeon/ | Not runnable here (needs samtools and bwa) |
| mc3 | github.com/OpenGenomics/mc3 | MC3 pipeline definitions | code/mc3/ | Documentation of CENTERS and FILTER |
| neusomatic | github.com/bioinform/neusomatic | Deep-learning caller and ensemble wrappers | code/neusomatic/ | Reference only |
| giab-stratifications | github.com/usnistgov/giab-stratifications | How the strata were built | code/giab-stratifications/ | Reference only |

## Resource gathering notes

### Search strategy
The paper-finder service was offline (it returned "not running"). The Semantic Scholar API returned HTTP 429. We therefore:
- searched Europe PMC with title and keyword queries (logs in `logs/epmc_search*.txt`);
- found PDFs through Unpaywall open-access locations;
- used Europe PMC full-text XML when publisher sites blocked automated downloads (Cloudflare or reCAPTCHA on Europe PMC, PMC, Springer, OUP, and bioRxiv).

For datasets, we browsed the SEQC2 and GIAB FTP trees directly and used the GDC API for MC3.

### Selection criteria
- Papers that report agreement between callers or sensitivity by VAF, variant type, or genomic context.
- Datasets where **the same tumor** has been called by several callers and each variant has a VAF. A truth set was preferred.

### Challenges
- Publisher bot protection blocked several PDFs. Fang 2021, Xiao 2021 (SEQC2, Nat Biotechnol), Ellrott 2018 (MC3), and Strelka2 are available as abstracts or secondary information only.
- The first MC3 download was corrupted because two download processes wrote to the file at once. It was downloaded again and checked (md5 `639ad8f8…`, and `gzip -t` passes).
- The SEQC2 ensemble TSVs have no TNscope, VarScan2, or LoFreq flags, even though the columns exist. TNscope is only available from the raw IL_T_1 VCF and the purity-series VCFs.
- The SomaticSniper VCFs have no FILTER values.

### Gaps and workarounds
- We could not re-run callers because no bioinformatics toolchain or Docker is available. We use the released per-caller outputs instead.
- For controlled VAF, the SEQC2 purity series replaces spike-in simulation.

## Recommendations for experiment design
1. **Primary dataset:** SEQC2 ensemble TSVs (D1) plus the purity series (D2). MC3 serves as replication (D3).
2. **Callers compared:**
   - SNVs: MuTect2, SomaticSniper, VarDict, MuSE, Strelka (plus TNscope and Lancet in the raw and purity VCFs).
   - Indels: MuTect2, VarDict, Strelka.
   - MC3: MuTect, MuSE, VarScan2, SomaticSniper, RADIA / Pindel, Indelocator, VarScan-indel.
3. **Metrics:**
   - fraction of variants called by all callers, and by a majority of callers;
   - normalised mean number of callers per variant;
   - pairwise Jaccard index and Fleiss' κ;
   - per-caller recall and precision.
   Report each metric on two denominators: all calls (union) and truth-set variants only. Stratify by VAF bin (<5, 5–10, 10–20, ≥20%) × SNV/indel × low-complexity status, with Wilson confidence intervals and logistic regression interaction tests.
4. **Code to reuse:** `tools/bedtools` for intersecting regions. The somaticseq parsers show how to extract VAF from each caller's VCF. `somaticseq.utilities.linguistic_sequence_complexity.LC` computes a complexity score from the hg38 sequence.

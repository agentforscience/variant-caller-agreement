#!/usr/bin/env bash
# Reproduce every result and figure (about 45 min on a 32-core machine; no GPU needed).
set -euo pipefail
cd "$(dirname "$0")"
source .venv/bin/activate
cd src
python d1_build.py      # SEQC2 ensemble tables -> results/d1/*.parquet   (~2 min)
python d1_analyze.py    # metrics, H1-H4 tests, GLMs -> results/d1_*.csv  (~25 min)
python d2_purity.py     # purity series                -> results/d2_*.csv (~3 min)
python d1b_raw.py       # raw IL_T_1 VCFs + TNscope     -> results/d1b_*.csv (~4 min; needs d1 + d2 module)
python d3_mc3.py        # TCGA MC3                      -> results/d3_*.csv (~7 min)
python verdicts.py      # consolidated verdicts         -> results/verdicts_*
python make_figures.py  # figures/*.png

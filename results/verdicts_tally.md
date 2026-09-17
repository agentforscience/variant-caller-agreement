## H1/H2: configurations where the hypothesis is supported (Holm-adjusted one-sided binomial)

|                                      | indel   | snv   |
|:-------------------------------------|:--------|:------|
| ('D1 SEQC2 ensemble', 'H1')          | 12/40   | 12/40 |
| ('D1 SEQC2 ensemble', 'H2')          | 37/40   | 40/40 |
| ('D1b raw VCFs (6/4 callers)', 'H1') | 2/6     | 2/6   |
| ('D1b raw VCFs (6/4 callers)', 'H2') | 6/6     | 6/6   |
| ('D2 purity series', 'H1')           | 0/30    | 11/32 |
| ('D2 purity series', 'H2')           | 36/36   | 36/36 |
| ('D3 TCGA MC3', 'H1')                | 0/2     | 2/2   |
| ('D3 TCGA MC3', 'H2')                | nan     | 2/2   |

## H3/H4: verdict counts across contexts

|                                       |   mixed (reversed in some bins) |   not supported |   partial |   reversed |   supported (all bins) |
|:--------------------------------------|--------------------------------:|----------------:|----------:|-----------:|-----------------------:|
| ('D1 SEQC2 ensemble', 'H3 indel>SNV') |                              33 |               0 |         5 |          6 |                      1 |
| ('D1 SEQC2 ensemble', 'H4 LC>nonLC')  |                               1 |               0 |        50 |          0 |                     39 |
| ('D1b raw VCFs', 'H3 indel>SNV')      |                               6 |               0 |         0 |          0 |                      0 |
| ('D1b raw VCFs', 'H4 LC>nonLC')       |                               0 |               0 |         8 |          0 |                      4 |
| ('D2 purity series', 'H3 indel>SNV')  |                              17 |               0 |         2 |          4 |                     13 |
| ('D2 purity series', 'H4 LC>nonLC')   |                               0 |               4 |        44 |          0 |                     24 |
| ('D3 TCGA MC3', 'H3 indel>SNV')       |                               1 |               0 |         1 |          0 |                      0 |
| ('D3 TCGA MC3', 'H4 LC>nonLC')        |                               0 |               1 |         3 |          0 |                      0 |

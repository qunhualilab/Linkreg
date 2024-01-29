# Linkreg
Linking candidate cis-regulatory elements to target genes

# Installation

Open terminal and type following:

```pip install Linkreg```


# Input & Output
**Input**:

1. Gene expression data file, separated by '\t'. Format: gene by biosample. Example:

|      | chromosome | gene_start | gene_end | strand | biosample1 | biosample2 | ... |
|------|------------|------------|----------|--------|------------|------------|-----|
| gene1 | chr1       | 50         | 60       | +      | 1.5        | 2.7        | ... |
| gene2 | chr1       | 70         | 80       | -      | 3.1        | 1.8        | ... |
| ...  | ...        | ...        | ...      | ...    | ...        | ...        | ... |

2. Track data file, separated by '\t'. Format: cCRE by biosample. Example:

| chromosome | start | end | biosample1 | biosample2 | ... |
|------------|-------|-----|------------|------------|-----|
| chr1       | 100   | 150 | 0.3        | 0.4        | ... |
| chr1       | 300   | 320 | 0.6        | 0.2        | ... |
| ...        | ...   | ... | ...        | ...        | ... |


**Output**:

cCRE-gene link score files. Format: pair by biosample. Example:

|      | chromosome | start | end | biosample1 | biosample2 | ... |
|------|------------|-------|-----|------------|------------|-----|
| cCRE1 | chr1       | 100   | 150 | 0.3        | 0.4        | ... |
| cCRE2 | chr1       | 300   | 320 | 0.6        | 0.2        | ... |
| ...  |            |       |     | ...        | ...        | ... |

# Example run

In your terminal, type

```Linkreg --expression_input ./data/expression.tpm --tracks_input ./data/DNase.bed ./data/H3K4me1.bed --output ./results/```

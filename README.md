# RNA-Seq-Differential-Expression-Analysis-Stroke-Diabetes-RNASeq-Analysis

"Differential Expression Analysis of Stroke and Diabetes RNA-Seq Data with Observational Training in a Stroke Animal Model".
using animal like mouse and rat sample and i used for stroke and diabtes type 1 mouse and his brain sample tissue  


End-to-end RNA-Seq differential gene expression analysis pipeline using HISAT2, HTSeq, DESeq2, FastQC, MultiQC, and Machine Learning for Stroke vs Diabetes transcriptomic analysis.
# RNA-Seq Differential Expression Analysis

## Overview

This project performs an end-to-end RNA-Seq differential expression analysis to compare transcriptomic profiles of Stroke, Type-2 Diabetes, and Healthy Control samples.

## Workflow

SRA Download

↓

FASTQ Conversion

↓

FastQC

↓

MultiQC

↓

HISAT2 Alignment

↓

HTSeq Count

↓

DESeq2

↓

Visualization

↓

Machine Learning Validation

## Dataset

- Organism: Homo sapiens
- Samples: 9
- Platform: Illumina
- GEO: GSE186844

## Tools

- FastQC
- MultiQC
- HISAT2
- Samtools
- HTSeq
- DESeq2
- R
- Python

## Results

✔ PCA Plot

✔ Volcano Plot

✔ Heatmap

✔ Differentially Expressed Genes

✔ Pathway Enrichment

✔ Machine Learning Classification

1. 413 vs 7,812 DEGs. The abstract says 413 DEGs are between stroke and diabetes, but your results show 413 for Stroke vs Control and 7,812 for Stroke vs Diabetes. Fix the abstract.

2.Type-1 vs Type-2 diabetes. The text says Type-2, but your code and plots use DIABETES_TYPE1. Choose one and make it consistent.
Dataset source. Your stroke and control samples (SRR166381xx/2xx) are from GSE186844, but the diabetes samples (SRR62166xx) look like they come from a different study. Say so in the README and check that the metadata matches the real GEO records.
3. Healthy control label. The code mixes Control and Healthy_control, so check it runs cleanly from start to finish.
Sample size. With n=3 per group, say clearly that results are exploratory.

## Author

Vajid Sayyad

M.Sc. Systems and Computational Biology

University of Hyderabad

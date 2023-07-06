# TFM

This GitHub repository houses the code and associated resources for my Master's thesis titled **"Studying fecal microbiota transplantation in HIV patients by novel functional profiling and enrichment method"**.

## Content

The following folders are included in this repository: 

- lablogs: Bash scripts to execute the performeed standard shotgun metagenomics pipeline using SLURM as job scheduler.
- emapper-profiler: Source code and documentation of the developed functional profiling tool
- DataAnalysisR: R scripts to analyze the functional potential of the profiled samples

## Abstract

The human gut microbiome has been postulated to play an important role in human immunodeficiency virus (HIV) patients. Recent studies based on 16S rRNA amplicon sequencing have shown that the gut microbiota can be significantly modified by fecal microbiome transplantation (FMT) from healthy donors. These changes have been only analyzed from a community structure perspective, identifying changes in the taxonomic composition of the samples. However, little is known about the functional changes produced by FMT and whether those changes might be considered positive for HIV patients. Reasons for such a knowledge gap is two-folds: First, there is a need for shotgun metagenomics sequencing data that enables functional analyses. Second, comprehensive analysis of functional patterns would require from novel and consistent bioinformatics methods that enable the comparison of fine-grained functional abundances across shotgun metagenomics samples.

Here, we address these challenges by developing a novel functional profiling tool, and applying to the study of a novel set of metagenomics data obtained by the REFRESH project, a longitudinal study on the effect of FMT on 13 HIV patients. The new tool developed in this project, named emapper-profiler, allows for the calculation of functional abundances across multiple samples, taking into account functional annotation from different sources such as eggNOG or KEGG databases.

By applying emapper-profiler to the REFRESH project metagenomics data, we uncovered valuable insights into the functional dynamics of the gut microbiota in HIV patients receiving FMT. Our analyses demonstrate significant functional differences between treated and control subjects, with a notable increase in the abundance of butyrate-producing bacteria compared to the placebo group. This suggests the positive impact of FMT in enhancing the functional potential of the gut microbiota. Additionally, the study identified several unknown gene families with differential abundance profiles, opening new avenues for the identification of novel biomarkers and future investigations in the treatment and diagnosis of inflammatory conditions associated to HIV.

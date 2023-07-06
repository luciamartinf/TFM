
# Shotgun Metagenomics Analysis Pipeline

Bash scripts for the execution of all the steps of the standard Shotgun Metagenomics Analysis Pipeline. 
   
![shotgun_pipeline2](https://github.com/luciamartinf/TFM/assets/56353778/c358993c-c438-4474-a018-a95d2c094c23)

## Detailed execution order: 

0. **Raw reads**, `00_script_raw_reads`: Bash script to create a text file containing all samples identifiers. 
1. **Pre-FastQC**, `01_lablog_pre_fastqc`: Quality analysis of the FastQ files with FastQC. Report is generated with multiqc. 
2. **Trimmomatic**, `02_lablog_trimmomatic`: Trimming of bad quality samples with Trimmomatic. 
3. **Post-FastQC**, `03_lablog_post_fastqc`: Quality analysis of the trimmed reads with FastQC. Report is generated with multiqc. 
4. **Assembly**, `04_lablog_assembly`: Assembly of curated reads into contigs using MetaSPAdes. 
5. **CoverM**, `05_lablog_coverm`: Coverage calculation of the contigs using CoverM. 
6. **eggNOG-mapper**, `06_lablog_eggnogmapper`: Prediction of open reading frames and subsequent functional annotation with eggNOG-mapper.
7. **mOTUs**, `07_lablog_motus`: Taxonomic profiling using mOTUs2. 
8. **HUMAnN3**, `08_lablog_humann3`: Complete functional and taxonomic profiling from curated reads using HUMAnN3.
9. **SqueezeMeta**, `09_lablog_squeezemeta`: Complete functional and taxonomic profiling from curated reads using SqueezeMeta.

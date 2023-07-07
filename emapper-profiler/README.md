
<img width="522" alt="image" src="https://github.com/luciamartinf/TFM/assets/56353778/17243875-e6f0-4da8-aa6e-51dfeb90c7bc">

### Lucía Martín Fernández
#

Functional profiling algorithm for eggNOG-mapper functional annotations. 


# Requirements:

## Software requirements

* Python 3.7 (or greater)

# BASIC USAGE:

`./emapper_profiler.py --input_dir data`

## Data preparation

* Execution of CoverM, for coverage calculations; 

```
coverm contig -1 reads/sampleid_R1.clean_qc_pair.fastq -2 reads/sampleid_R2.clean_qc_pair.fastq -r assembly/sampleid/contigs.fasta -o data/sample_coverage_values -m length count reads_per_base mean trimmed_mean tpm rpkm
```

* Execution of eggNOG-mapper, for functional annotations, and default settings from the web server.

```
emapper.py -m diamond --itype metagenome --genepred prodigal -i assembly/sampleid/contigs.fasta --output_dir data --output sampleid --pident 40.0 --evalue 0.001 --score 60.0 --query_cover 20.0 --subject_cover 20.0
```

* Execution of eggNOG-mappaer with novel gene families as reference database, and default settings from the web server.

```
emapper.py -m novel_fams --itype metagenome --genepred prodigal -i assembly/sampleid/contigs.fasta --output_dir data/novel_families/ --output sample_id --pident 40.0 --evalue 0.001 --score 60.0 --query_cover 20.0 --subject_cover 20.0
```


### Input Files

Execution of CoverM, for coverage calculations; and eggNOG-mapper, for functional annotations. 

A sample  the files needed to run *emapper-profiler* are  need to be 
and the ones generated are contained in this github repository and detailed below. 

| **Files**                           | **Description**                                                                                                 |                                                   
|:----------------------------------------|:----------------------------------------------------------------------------------------------------------------|
|`ArabidopsisSubNetwork_GeneList.txt`                                |  txt file stored. Contains a list of Arabidopsis Thaliana co-expressed genes                          |                  
|`GFF_genes.gff3`                                  | New GFF3 file. GFF3 report that contains all the CTTCTT regions found in the genes from the list. The coordinates are relative to the gene.                 |                                                               
|`GFF_chr.gff3`                                  | New GFF3 file. GFF3 report that contains all the CTTCTT regions found in the genes from the list. The coordinates are relative to the chromosome.                |   
|`no_CTTCTT_genes_report.txt`                                  | New txt file. Report that contains all the genes that don't contain any CTTCTT region in any of their exons                 |     

### How to install it? 


## Source Files

### main.rb

The following tasks are performed:

-   Reads the genelist file with Arabidopsis Thaliana genes and stores them in an Array
-   Creates a Seq object and gets all the useful information for each gene in the genelist file 
-   Genereates the files previously mentioned

### FileMaster.rb

Class for File Management. Methods defined:

- `get_genelist_from_file`: Class method to read genelist file and create an array withall the genes
- `fetch`: Class method fetch to access an URL via code, donated by Mark Wilkinson.
- `generate_no_report`: Class method to write report with all the genes that don't contain any CTTCTT region in any of their exons
- `generate_gff3_report`: Class method to write GFF3 report with all the CTTCTT regions found in the genes from the list. The coordinates are relative to the gene.
- `generate_gff3_report_chr`: Class method to write GFF3 report with all the CTTCTT regions found in the genes from the list. The coordinates are relative to the chromosome.

### Seq.rb 

Class to represent all of the information associated with a gene. Methods defined:

- `initialize`: Definition of initialize method.
- `get_no_genes`: Class method to get the the class variable @@no_cttctt_genes.
- `get_bio`: Instance method to obtain information for the current genes.
- `get_coords_chr`: Instance method to obtain the chromose number and the chromosome relative coordinates of the gene
- `add_new_feature`: Instance method to create new CTTCTT_region feature and add it to the Bio::Sequence object @sequence.
- `match_exons`: Instance method to find every CTTCTT region in the exons of the current gene on both the + and - strands, adds a new feature to its @sequence and updates @@no_cttctt_genes array


## Other files in this repository

- Ruby scripts are commented with YARD and respective documentation can be found in this repository. 
- `Image1`: Screenshot of the EnsemblPlants webpage for the AT2G46340 gene zoomed in in Chromosome 2 19027247-19027302 region showing the CTTCTT regions can also be found in this repository. 
- `Image2`: Screenshot of the EnsemblPlants webpage for the AT2G46340 gene zoomed in in Chromosome 2 19027016-19027468 region showing the CTTCTT regions can also be found in this repository. 


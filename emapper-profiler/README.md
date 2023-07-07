
<img width="522" alt="image" src="https://github.com/luciamartinf/TFM/assets/56353778/17243875-e6f0-4da8-aa6e-51dfeb90c7bc">

### Lucía Martín Fernández
#

Functional profiling algorithm for eggNOG-mapper functional annotations. 


## Basic usage:

```
./emapper_profiler.py --input_dir data
```
`data` is a folder containing both eggNOG-mapper and CoverM output results from all samples. In the next section, how to generate this files is described.  

### Data preparation

* Execution of CoverM, for coverage calculations; 

```
coverm contig -1 reads/sampleid_R1.clean_qc_pair.fastq -2 reads/sampleid_R2.clean_qc_pair.fastq \
-r assembly/sampleid/contigs.fasta -o data/sample_coverage_values \
-m length count reads_per_base mean trimmed_mean tpm rpkm
```

* Execution of eggNOG-mapper, for functional annotations, and default settings from the web server.

```
emapper.py -m diamond --itype metagenome --genepred prodigal -i assembly/sampleid/contigs.fasta \
--output_dir data --output sampleid \ 
--pident 40.0 --evalue 0.001 --score 60.0 --query_cover 20.0 --subject_cover 20.0
```

* Execution of eggNOG-mappaer with novel gene families as reference database, and default settings from the web server.

```
emapper.py -m novel_fams --itype metagenome --genepred prodigal -i assembly/sampleid/contigs.fasta \
--output_dir data/novel_families/ --output sample_id \
--pident 40.0 --evalue 0.001 --score 60.0 --query_cover 20.0 --subject_cover 20.0
```

## Output files

| **Files**                           | **Description**                                                                                                 |                                                   
|:----------------------------------------|:----------------------------------------------------------------------------------------------------------------|
|`ko_relabundance.tsv`                                |  Relative abundance per sample of KEGG orthologs acompanied by its description and symbol. Unmapped proportion of functions is included                      |                  
|`og_relabundance.tsv`                                  | Relative abundance per sample of orthologous groups from the eggNOG database annotated at the kingdom taxonomic level acompanied by its description. Unmapped proportion of functions is included          |                                                               
|`ko_totalabundance.tsv`                                | Abundance per sample of KEGG orthologs acompanied by its description and symbol                      |                  
|`og_totalabundance.tsv`                                  | Abundance per sample of orthologous groups from the eggNOG database annotated at the kingdom taxonomic level acompanied by its description          |  
|`pathway_coverage.tsv`                                  |  KEGG pathway's completness percentage per sample according to the annotated KEGG ortholog. Pathway's description is included.                      |                  
                                                       
|`nf_relabundance.tsv` (if `--novel_fam` flag is on)                               |  Relative abundance per sample of novel gene families. Unmapped proportion of functions is included          |                                              |                  
|`nf_totalabundance.tsv` (if `--novel_fam` flag is on)                                 | Abundance per sample of novel gene families          |                                                               

## USAGE Example

The [data_resume](data_resume) folder contains an example files required for *emapper-profiler* execution. The [results_resume](results_resume) folder contains the outputs generated when executing *emapper-profiler* with the following command line:

```
./emapper_profiler.py --inputdir data_resume --outputdir results_tpm --filer_euk --unit tpm --filter_virus --novel_fam

```

Complete results generated for all dataset are included in the [results](results) folder. 

## Software requirements

* Python 3.7 (or greater)

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


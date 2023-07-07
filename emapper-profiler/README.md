
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

## Arguments

`--inputdir DIR`,`-i DIR`
    Input directory containing CoverM and eggNOG-mapper results. Required. 

`--outputdir DIR`,`-o DIR`, default ='results',
   Output directory to store the generated tsv files. By default, output directory is called results.

    add('--coverm_suffix', '-c', default='_coverage_values',
        help='The suffix of your CoverM files')

    add('--version', action='version', version='%(prog)s 0.0.1')

    add('--sample_file', '-s',
        help = 'txt file with a list of the samples that need to be processed. If not provided, the sample list will be generated given the input files')

    add('--kegg_dict', '-k',
        help = 'Directory containing KEGG_pathway and KEGG_kos json dictionaries')

    add('--unit', '-u', type=str, choices=['tpm', 'rpkm', 'tm'], default='rpkm',
        help = 'Specify the output units for the relative abundance results. Default is rpkm')

    # on/off flag
    add('--filter_euk', '-e', action='store_true',
        help='Remove eukaryotes')

    add('--filter_virus', '-v', action='store_true',
        help='Remove viruses')

    add('--novel_fam', '-f', action='store_true',
        help='Add novel families')

    add('--nf_dir',
        help='Input directory that contains the novel families annotations. By default considers a directory called novel_families inside the input directory')


## Output files

| **Files**                           | **Description**                                                                                                 |                                                   
|:----------------------------------------|:----------------------------------------------------------------------------------------------------------------|
|`ko_relabundance.tsv`                                |  Relative abundance per sample of KEGG orthologs acompanied by its description and symbol. Unmapped proportion of functions is included                      |                  
|`og_relabundance.tsv`                                  | Relative abundance per sample of orthologous groups from the eggNOG database annotated at the kingdom taxonomic level acompanied by its description. Unmapped proportion of functions is included          |                                                               
|`ko_totalabundance.tsv`                                | Abundance per sample of KEGG orthologs acompanied by its description and symbol                      |                  
|`og_totalabundance.tsv`                                  | Abundance per sample of orthologous groups from the eggNOG database annotated at the kingdom taxonomic level acompanied by its description          |  
|`pathway_coverage.tsv`                                  |  KEGG pathway's completness percentage per sample according to the annotated KEGG ortholog. Pathway's description is included.                      |     
|                |         **if `--novel_fam` flag is on**                                                                                                 |      
|`nf_relabundance.tsv`                                  |  Relative abundance per sample of novel gene families. Unmapped proportion of functions is included          |                                                       
|`nf_totalabundance.tsv`                                  | Abundance per sample of novel gene families          |                                                               

## Usage Example

The [data_resume](data_resume) folder contains an example files required for *emapper-profiler* execution. The [results_resume](results_resume) folder contains the outputs generated when executing *emapper-profiler* with the following command line:

```
./emapper_profiler.py --inputdir data_resume --outputdir results_tpm --filer_euk --unit tpm --filter_virus --novel_fam

```

Complete results generated for all dataset are included in the [results](results) folder. 

## Software requirements

* Python 3.7 (or greater)




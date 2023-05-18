#!/usr/bin/python
# -*- coding: utf-8 -*-

from eggnog_classes import Eggnog_sample
from coverm_classes import Coverm_sample
from novelfam_classes import NovelFam_sample
from ko_functions import ko_abundance, write_tsv, calculate_KEGG_pathway_completeness, og_abundance, write_json
from novelfam_fun import check_novelfam, nf_abundance, check_all_novelfam
from arg_parse import check_arg
import json
import os
import sys
import re
from total_ko_functions import total_ko_abundance

#######################
### PARSE ARGUMENTS ###
#######################

# Example: python main.py -i data_resume/ -s sample_file.txt -k /Users/lucia/Desktop/TFM/scripts/parse_KEGGpathway_db/KEGG_KOs_dict.txt

arguments = check_arg(sys.argv[1:])

## Get Kegg pathways dictionary

if arguments.kegg_dict:
    KEGG_dict_file = arguments.kegg_dict
else:
    # generar el diccionario a partir del archivo KOs.kegs_cleaned ? No tiene sentido porque tengo que proporcionar tb ese archivo
    KEGG_dict_file = "/Users/lucia/Desktop/TFM/scripts/parse_KEGGpathway_db/KEGG_pathway_dict.txt"

kos_dict_file = "/Users/lucia/Desktop/TFM/scripts/parse_KEGGpathway_db/KEGG_kos_dict.txt"

# Load the dictionary from the json file
with open(KEGG_dict_file, "r") as file:
    KEGG_dict = json.load(file)

with open(kos_dict_file, "r") as file:
    kos_dict = json.load(file)

## Get input directory and its files

# if arguments.inputdir:
#     input_eggnog = arguments.inputdir
#     input_cov = arguments.inputdir
    
# elif arguments.inputdir_eggnogmapper and arguments.inputdir_coverage:
#     input_eggnog = arguments.inputdir_eggnogmapper
#     input_cov = arguments.inputdir_coverage
# else:
#     print("error: either --inputdir -i or --input_eggnogmapper -e together with --input_coverage -c arguments required")
    
#     #help
#     sys.exit()
    
inputdir = arguments.inputdir
files = os.listdir(inputdir)

## Get output directory 

if arguments.outputdir:
    outputdir = arguments.outputdir
else: 
    outputdir = 'results'

if not os.path.isdir(outputdir):
    os.mkdir(outputdir)

ko_abun_file = os.path.join(outputdir, 'ko_abundance.tsv')
og_abun_file = os.path.join(outputdir, 'og_abundance.tsv')
nf_abun_file = os.path.join(outputdir, 'nf_abundance.tsv')
path_cov_file = os.path.join(outputdir, 'pathway_coverage.tsv')


## Get all sample names

if arguments.sample_file :
    with open(arguments.sample_file) as f:
        sample_list = f.read().splitlines()
else:
    sample_list = []
    for file_name in files:
        basename = os.path.basename(file_name)
        if '.emapper.annotations' in basename:
            samplename = re.sub(r'.emapper.annotations', '', basename)
            sample_list.append(samplename)


## Option arguments

remove_euk = arguments.filter_euk
novel_fam = arguments.novel_fam

#########################
### PROGRAM EXECUTION ###
#########################

ko_abundance_all = {}
path_coverage = {}
og_abundance_all = {}
nf_abundance_all = {}


for sample in sample_list:

    print('Starting processing sample {}'.format(sample))

    # Define eggnog and coverm filenames
    eggnog_file = os.path.join(inputdir, sample + '.emapper.annotations')
    coverm_file = os.path.join(inputdir, sample + '_coverage_values')
    

    # Load eggnog and coverm samples    
    eggnog_sample = Eggnog_sample(eggnog_file, sample, remove_euk)
    eggnog_sample.load_sample()
    coverm_sample = Coverm_sample(coverm_file, sample)
    coverm_sample.load_sample()

    # Add sample abundance and pathways coverage to complete dictionary
    #ko_abundance_all, all_ko_list = ko_abundance(eggnog_sample, coverm_sample, sample, ko_abundance_all, sample_list, kos_dict) 
    ko_abundance_all, all_ko_list = total_ko_abundance(eggnog_sample, coverm_sample, sample, ko_abundance_all, sample_list, kos_dict) 
    path_coverage = calculate_KEGG_pathway_completeness(KEGG_dict, all_ko_list, path_coverage, sample, sample_list)
    og_abundance_all = og_abundance(eggnog_sample, coverm_sample, sample, og_abundance_all, sample_list)

    if novel_fam :
        novelfam_file = os.path.join(inputdir+'/novel_families/', sample + '.emapper.annotations')
        novelfam_sample = NovelFam_sample(novelfam_file, sample)
        novelfam_sample.load_sample()

        #repeated_queries = check_all_novelfam(eggnog_sample, novelfam_sample, sample)
        nf_abundance_all = nf_abundance(novelfam_sample, coverm_sample, sample, nf_abundance_all, sample_list, eggnog_sample)

    print('Finished processing sample {}'.format(sample))


##########################
### WRITE OUTPUT FILES ###
##########################

print('Writing files. Almost done...')

ko_abundance_all = dict(sorted(ko_abundance_all.items()))
path_coverage = dict(sorted(path_coverage.items()))
og_abundance_all = dict(sorted(og_abundance_all.items()))
nf_abundance_all = dict(sorted(nf_abundance_all.items()))

header='KEGG_ko\tDescription\t'+ '\t'.join(sample_list)
write_tsv(ko_abundance_all, ko_abun_file, header, sample_list, des = True)
write_json(os.path.join(outputdir, 'ko_abundance.json'), ko_abundance_all)

header='KEGG_pathway\tDescription\t'+ '\t'.join(sample_list)
write_tsv(path_coverage, path_cov_file, header, sample_list, des = True)

header='OG\tKingdom\tDescription\t'+ '\t'.join(sample_list)
write_tsv(og_abundance_all, og_abun_file, header, sample_list, des = True, king = True)

header='Novel_Fam\t'+ '\t'.join(sample_list) + '\tCOG_match'
write_tsv(nf_abundance_all, nf_abun_file, header, sample_list, cog=True)



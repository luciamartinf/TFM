#!/usr/bin/python
# -*- coding: utf-8 -*-

from eggnog_classes2 import Eggnog_sample
from novelfam_classes2 import NovelFam_sample
from ko_functions2 import write_tsv, write_json, read_coverm_as_nested_dict, extract_orf_lengths
#from novelfam_fun import check_novelfam, nf_abundance, check_all_novelfam
from arg_parse import check_arg
import json
import os
import sys
import re

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
    
inputdir = arguments.inputdir
if not os.path.isdir(inputdir):
    print('Error: Input directory does not exist')

files = os.listdir(inputdir)

if arguments.nf_dir:
    novelfam_dir = arguments.nf_dir
else:
    novelfam_dir = os.path.join(inputdir+'/novel_families/')

## Get output directory 

outputdir = arguments.outputdir

## Option arguments

remove_euk = arguments.filter_euk
novel_fam = arguments.novel_fam


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

## select units

units_option = arguments.unit

## coverm suffix

coverm_suffix = arguments.coverm_suffix

# Create output directory if it does not exist
if not os.path.isdir(outputdir):
    os.mkdir(outputdir)

ko_abun_file = os.path.join(outputdir, 'ko_abundance.tsv')
og_abun_file = os.path.join(outputdir, 'og_abundance.tsv')
path_cov_file = os.path.join(outputdir, 'pathway_coverage.tsv')

#########################
### PROGRAM EXECUTION ###
#########################

ko_abundance_all = {}
#ko_abundance_all['UNMAPPED'] = {}
path_coverage = {}
og_abundance_all = {}
nf_abundance_all = {}


Eggnog_sample.init_unit(units_option)
Eggnog_sample.init_sample_list(sample_list)

if novel_fam :
    nf_abun_file = os.path.join(outputdir, 'nf_abundance.tsv')
    NovelFam_sample.init_unit(units_option)
    NovelFam_sample.init_sample_list(sample_list)

for sample in sample_list:

    print('Starting processing sample {}'.format(sample))

    # Define eggnog and coverm filenames
    eggnog_file = os.path.join(inputdir, sample + '.emapper.annotations')
    coverm_file = os.path.join(inputdir, sample + coverm_suffix)
    genepred_file = os.path.join(inputdir, sample + '.emapper.genepred.fasta')
    

    # Load eggnog and coverm samples    
    coverm_dict = read_coverm_as_nested_dict(coverm_file, Eggnog_sample.calc_unit) 
    orf_dict, total = extract_orf_lengths(genepred_file, coverm_dict, Eggnog_sample.calc_unit)

    eggnog_sample = Eggnog_sample(eggnog_file, total, sample, remove_euk)
    og_abundance_all = eggnog_sample.load_sample(orf_dict, og_abundance_all)
    
    # Add sample abundance and pathways coverage to complete dictionary
    og_abundance_all = eggnog_sample.calculate_og_abundance(og_abundance_all)
    ko_abundance_all = eggnog_sample.calculate_ko_abundance(ko_abundance_all, kos_dict)
    path_coverage = eggnog_sample.calculate_KEGG_pathway_completeness(path_coverage, KEGG_dict)
    
    # Repeat process for novel families
    if novel_fam :
        novelfam_file = os.path.join(novelfam_dir, sample + '.emapper.annotations')
        novelfam_sample = NovelFam_sample(novelfam_file, total, sample)
        novelfam_sample.load_sample(orf_dict)
        nf_abundance_all = novelfam_sample.calculate_nf_abundance(nf_abundance_all)

    #     #repeated_queries = check_all_novelfam(eggnog_sample, novelfam_sample, sample)
          
    print('Finished processing sample {}'.format(sample))


##########################
### WRITE OUTPUT FILES ###
##########################

print('Writing files. Almost done...')

# ko_abundance_all = dict(sorted(ko_abundance_all.items()))
# path_coverage = dict(sorted(path_coverage.items()))
# og_abundance_all = dict(sorted(og_abundance_all.items()))
# nf_abundance_all = dict(sorted(nf_abundance_all.items()))

header='KEGG_ko\tDescription\t'+ '\t'.join(sample_list)
write_tsv(ko_abundance_all, ko_abun_file, header, sample_list, des = True)
write_json(os.path.join(outputdir, 'ko_abundance.json'), ko_abundance_all)

header='KEGG_pathway\tDescription\t'+ '\t'.join(sample_list)
write_tsv(path_coverage, path_cov_file, header, sample_list, des = True)

header='OG\tKingdom\tDescription\t'+ '\t'.join(sample_list)
write_tsv(og_abundance_all, og_abun_file, header, sample_list, des = True, king = True)

if novel_fam :
    header='Novel_Fam\t'+ '\t'.join(sample_list) #+ '\tCOG_match'
    write_tsv(nf_abundance_all, nf_abun_file, header, sample_list) #, cog=True)



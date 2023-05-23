#!/usr/bin/python
# -*- coding: utf-8 -*-

import json
import csv 
import re


def unique(list1):

    """
    Function to transform a list so its elements are not repeated
    """
 
    # insert the list to the set
    list_set = set(list1)
    # convert the set to the list
    unique_list = (list(list_set))

    return unique_list

def clean_header_names(name):

    """
    Function to clean the header names of the CoverM file
    """

    if ' ' in name:
        name = name.split(' ', 1)[1]
        clean_name = re.sub(' ', '_', name)
    return clean_name

def read_coverm_as_nested_dict(file_path, unit):

    """
    Function to read the CoverM file as nested dictionary
    """

    nested_dict = {}
    #total = 0

    with open(file_path, 'r') as tsv_file:
        reader = csv.DictReader(tsv_file, delimiter='\t') # read tsv file

        for row in reader: # row is a dictionary 
            contig = row['Contig'] # save contig
            nested_dict[contig] = {}

            for key, value in row.items():
                if key != 'Contig': 
                    clean_key = clean_header_names(key)
                    if clean_key == unit: # get desired unit
                        
                        nested_dict[contig][unit] = value # add value to nested dictionary
                        #total += float(value) # sum value to total

    return nested_dict

def get_ko_list(raw_ko):

    """
    Function to obtain a clean ko list
    """

    kos_list = []
    if ',' in raw_ko:
        kos = raw_ko.split(',')
    else: 
        kos = [raw_ko]

    for ko in kos:
        ko_id = ko.split(':')[1]
        
        kos_list.append(ko_id)

    return kos_list

def find_basal(eggnog_ogs):

    """ 
    Function to obtain the basal orthologous group and its kingdom level 
    """

    # Example of raw og line : COG0705@1|root,COG0705@2|Bacteria,1TSBP@1239|Firmicutes,249D8@186801|Clostridia
    ogs = eggnog_ogs.split(',')
    for og in ogs:
        group, kingdom = og.split('|')
        og_id, tax_code = group.split('@')
        if tax_code in ['2', '2759', '2157']: # if tax_code corresponds to Bacteria, Eukaryota or Arkea 
            # self.og = og_id
            # self.kingdom = kingdom
            return og_id, kingdom

def check_unmapped(dict1:dict, sample_list:list, des=True):
    
    if 'UNMAPPED' not in dict1.keys():
        dict1['UNMAPPED'] = {}
        for sample in sample_list:
            dict1['UNMAPPED'][sample] = 0
        if des == True:
            dict1['UNMAPPED']['description'] = '@'

    return dict1

def write_tsv (dictionnary, out_file, header, sample_list, des = False, king = False, ko= False, cog= False):

    '''
    Description:
        Function to extract the relevant part of result.txt file
    Input:
        result.txt file
    Return:
        dictionary
    '''
	
    with open (out_file,'w') as fo:
        fo.write(header+ '\n')
        for key in dictionnary.keys():
            fo.write(key)
            if king:
                fo.write('\t'+str(dictionnary[key]['kingdom']))
            if des:
                fo.write('\t'+str(dictionnary[key]['description']))
            if ko:
                fo.write('\t'+str(dictionnary[key]['ko']))
            for sample in sample_list:
                fo.write('\t'+ str(dictionnary[key][sample]))
            if cog:
                fo.write('\t'+str(dictionnary[key]['cog']).strip(', '))
            fo.write('\n')
    return


def write_json(file, dictionary):
    with open(file, 'w') as fp:
        json.dump(dictionary, fp)


def extract_orf_lengths(fasta_file, coverm_dict:dict, unit):
    
    orf_dict = {}  # Dictionary to store sequence names, lengths and abundances
    total = 0
    orf_per_contig = {}

    with open(fasta_file, 'r') as file:
        lines = file.readlines()

        for line in lines:
            line = line.strip()

            if line.startswith('>'):
                header = line[1:]  # Remove the leading '>'
                orf_name = header.split(' ')[0]  # Extract the sequence name
                info = header.split('#')[1:]  # Extract the sequence info after '#'

                # Extract the start and end positions from the sequence info
                start_pos = int(info[0].strip())
                end_pos = int(info[1].strip())

                # Calculate the sequence length
                length = end_pos - start_pos + 1

                orf_dict[orf_name] = {}
                orf_dict[orf_name]['length'] = length

                contig_name = re.sub(r'_[0-9]*$', '', orf_name)
                abundance = float(coverm_dict[contig_name][unit])
                orf_dict[orf_name]['abundance'] = abundance
                total += abundance

    #             if unit == 'Trimmed Mean':
    #                 total += abundance
    #             else:
    #                 try:
    #                     orf_per_contig[contig_name] += 1
    #                 except:
    #                     orf_per_contig[contig_name] = 1
                
    # if unit == 'RPKM':
    #     for orf in orf_dict.keys():
    #         contig_name = re.sub(r'_[0-9]*$', '', orf)
    #         abun_1 = float(orf_dict[orf]['abundance'])
    #         abundance = abun_1 / int(orf_per_contig[contig_name])
    #         orf_dict[orf]['abundance'] = abundance
    #         total += abundance

    return orf_dict, total
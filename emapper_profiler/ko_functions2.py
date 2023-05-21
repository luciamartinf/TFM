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
    if ' ' in name:
        name = name.split(' ', 1)[1]
        clean_name = re.sub(' ', '_', name)
    return clean_name

def read_coverm_as_nested_dict(file_path, unit):
    nested_dict = {}
    total = {}
    total = 0

    with open(file_path, 'r') as tsv_file:
        reader = csv.DictReader(tsv_file, delimiter='\t')

        for row in reader:
            contig = row['Contig']
            nested_dict[contig] = {}

            for key, value in row.items():
                if key != 'Contig':
                    clean_key = clean_header_names(key)
                    if clean_key == unit:
                        nested_dict[contig][clean_key] = value
                        total += float(value)

    return nested_dict, total

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

def check_unmapped(dict1, sample_list):
 
    if 'UNMAPPED' in dict1.keys():
        dict1['UNMAPPED']['description'] = '@'
    
    else:
        dict1['UNMAPPED'] = {}
        for sample in sample_list:
            dict1['UNMAPPED'][sample] = 0
    
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
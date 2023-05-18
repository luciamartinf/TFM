#!/usr/bin/python
# -*- coding: utf-8 -*-

import json

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


def unique(list1):

    """
    Function to transform a list so its elements are not repeated
    """
 
    # insert the list to the set
    list_set = set(list1)
    # convert the set to the list
    unique_list = (list(list_set))

    return unique_list



def find_basal(raw_og):
    og_org = {}
    ogs = raw_og.split(',')
    for og in ogs:
        group, kingdom = og.split('|')
        og_id, tax_code = group.split('@')
        if tax_code in ['2', '2759', '2157']:
            og_org[og_id] = kingdom
    return og_org

def add_og_abundance(og_abundance_sample:dict, og_description:dict, og_king:dict, raw_ogs, abundance, total_abun, description):

    """
    calculate the OG abundance
    """
    og_king_small = find_basal(raw_ogs)

    for og, king in og_king_small.items():

        try: 
            og_abundance_sample[og] += float(abundance)
            
        except:
            og_abundance_sample[og] = float(abundance)
            og_description[og] = description
            og_king[og] = king
        
        total_abun += float(abundance)

  
    return og_abundance_sample, total_abun

def og_abundance(eggnog_sample, coverm_sample, samplename, rel_og_abun:dict, sample_list):

    """
    Calculate every ko total abundance in a sample. 
    This function relies on contigs and orf being sorted.
    Takes rpkm value as abundance
    """
    og_description = {}
    og_king = {}
    i = 0

    total_og_abun = 0

    og_abundance_sample = {}

    for eggnog_row in eggnog_sample.rows:

        while coverm_sample.rows[i].contig != eggnog_row.contig:
            i += 1
        
        eggnog_row.abundance = coverm_sample.rows[i].rpkm
        og_abundance_sample, total_og_abun = add_og_abundance(og_abundance_sample, og_description, og_king, eggnog_row.eggnog_ogs, eggnog_row.abundance, total_og_abun, eggnog_row.description)
    
    for og in og_abundance_sample.keys():
        try:
            rel_og_abun[og]['description'] = og_description[og]
            
        except:
            rel_og_abun[og] = {}
            rel_og_abun[og]['kingdom'] = og_king[og]  
            rel_og_abun[og]['description'] = og_description[og]
            for sample_id in sample_list:
                rel_og_abun[og][sample_id] = 0

        rel_og_abun[og][samplename] = og_abundance_sample[og]/total_og_abun      

    return rel_og_abun

def add_ko_abundance(ko_abundance_sample:dict, kos, abundance, total_abun):

    """
    add abundance to each ko on the list
    """
    for ko_id in kos:
        # initialize dictionary for each ko
        if not ko_id in ko_abundance_sample.keys():
            ko_abundance_sample[ko_id] = 0
        
        
        # add ko abundance 
        ko_abundance_sample[ko_id] += float(abundance)
        total_abun += float(abundance)

  
    return ko_abundance_sample, total_abun

def ko_abundance(eggnog_sample, coverm_sample, samplename, rel_ko_abun:dict, sample_list, kos_dict):

    """
    Calculate every ko total abundance in a sample. 
    This function relies on contigs and orf being sorted.
    Takes rpkm value as abundance
    """
    # CHECK WHEN A SEED CORRESPOND TO SEVERAL KOS
    # more_kos = []
    # total_kos = 0 

    i = 0
    
    global_ko_list = []

    total_ko_abun = 0

    ko_abundance_sample = {}

    for eggnog_row in eggnog_sample.rows:

        
        if eggnog_row.kegg_ko != '-':
            # append kos from sample to global ko list
            kos = get_ko_list(eggnog_row.kegg_ko)
                
            # CHECK WHEN A SEED CORRESPOND TO SEVERAL KOS
            # total_kos +=1 
            # if len(kos) > 1:
            #     more_kos.append(len(kos))
            #     #print(len(kos), eggnog_row.query)


            for ko in kos:
                global_ko_list.append(ko)

            while coverm_sample.rows[i].contig != eggnog_row.contig:
                i += 1
                
            eggnog_row.abundance = coverm_sample.rows[i].rpkm
            ko_abundance_sample, total_ko_abun = add_ko_abundance(ko_abundance_sample, kos, eggnog_row.abundance, total_ko_abun)

    for ko_id in ko_abundance_sample.keys():
        if not ko_id in rel_ko_abun.keys():
            rel_ko_abun[ko_id] = {}
            try: 
                rel_ko_abun[ko_id]['description'] = kos_dict[ko_id]['description']
            except:
                rel_ko_abun[ko_id]['description'] = '@'
            
            for sample_id in sample_list:
                rel_ko_abun[ko_id][sample_id]=0
    
        rel_ko_abun[ko_id][samplename] = ko_abundance_sample[ko_id]/total_ko_abun
    
    # CHECK WHEN A SEED CORRESPOND TO SEVERAL KOS
    # total = len(more_kos)
    # two = more_kos.count(2)
    # three = more_kos.count(3)
    # rest = total - two - three 
    # print(total_kos, round((total/total_kos)*100,2), total, two, three, rest)


    return rel_ko_abun, unique(global_ko_list)




def calculate_KEGG_pathway_completeness(KEGG_dict, global_KOs_list, path_coverage, sample, sample_list): #output_file):

    """
    Function to calculate KEGG pathway completeness of each sample 
    """

    kegg_cov_dict = {}

    #f = open(output_file, "w")

    for kegg_p, annotation in KEGG_dict.items():
        pathway_description = annotation[0]
        kegg_number = annotation[1]
        kegg_cov_dict[kegg_p] = 0   
        for ko in annotation[2]:
            ko_id = ko['KO']
            # symbol = ko['symbol']
            # description = ko['description']
            if ko_id in global_KOs_list:
                kegg_cov_dict[kegg_p] +=1
        
        if kegg_cov_dict[kegg_p] != 0:
            coverage = kegg_cov_dict[kegg_p]/kegg_number
            if coverage > 0.1: # cut-off proporcion
                if not kegg_p in path_coverage.keys():
                    path_coverage[kegg_p] = {}
                    path_coverage[kegg_p]['description'] = pathway_description
                    for sample_id in sample_list:
                        path_coverage[kegg_p][sample_id] = 0
                
                path_coverage[kegg_p][sample] = coverage
                #f.write("K{}\t{}\t{}\t{}\t{}\n".format(kegg_p,pathway_description,str(kegg_number),str(kegg_cov_dict[kegg_p]),str(coverage))) #, kos_list)			
    return path_coverage



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
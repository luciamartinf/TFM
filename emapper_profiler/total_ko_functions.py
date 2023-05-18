#!/usr/bin/python
# -*- coding: utf-8 -*-


def unique(list1):

    """
    Function to transform a list so its elements are not repeated
    """
 
    # insert the list to the set
    list_set = set(list1)
    # convert the set to the list
    unique_list = (list(list_set))

    return unique_list 

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

def total_ko_abundance(eggnog_sample, coverm_sample, samplename, rel_ko_abun:dict, sample_list, kos_dict):

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
                
            eggnog_row.abundance = coverm_sample.rows[i].trimmed_mean
                
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
    
        rel_ko_abun[ko_id][samplename] = ko_abundance_sample[ko_id]
    
    # CHECK WHEN A SEED CORRESPOND TO SEVERAL KOS
    # total = len(more_kos)
    # two = more_kos.count(2)
    # three = more_kos.count(3)
    # rest = total - two - three 
    # print(total_kos, round((total/total_kos)*100,2), total, two, three, rest)


    return rel_ko_abun, unique(global_ko_list)



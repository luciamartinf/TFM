#!/usr/bin/python
# -*- coding: utf-8 -*-

from dataclasses import dataclass
import re
import os
from ko_function2 import get_ko_list, find_basal, check_unmapped

@dataclass
class Eggnog_orf(object):

    """
    Dataclass to store eggnog's attributes corresponding to one orf
    """

    query: str
    eggnog_ogs: str
    og: str 
    kingdom: str 
    description: str # = field(repr=False)
    preferred_name:	str
    #GOs: str	
    kegg_ko: list
    kegg_pathway: str
    #KEGG_Module: str	
    contig: str #= field(init=False) 
    abundance: float

    def __str__(self):

        """
        Instance method to define the print format of the instance
        """

        row_list = [self.query, self.seed_ortholog, self.eggnog_ogs, self.max_annot_lvl, self.kegg_ko, self.kegg_pathway, self.contig]
        s = '\t'.join([str(item) for item in row_list])
        
        return s
    
            
class Eggnog_sample(object):

    """
    Class to store all rows from an eggnog-mapper file. 
    Each row is a eggnog_orf instance
    """

    option_unit = None # the option unit as given in the arguments
    calc_unit = None # the unit we use for calculations
    sample_list = None
    

    def __init__(self, filename: str, total_sample:dict, samplename = None, remove_euk = False) -> None:
        self.rows = []
        self.filename = filename
        self.remove_euk = remove_euk
        if samplename == None :
            self.samplename = self.define_samplename() # function to extract samplename from filename
        else: 
            self.samplename = samplename
        
        self.og_abundance = {}
        self.ko_abundance = {}
        self.total = float(total_sample)
        self.mapped = 0
        self.global_ko_list = []

        # self.query_list = []
        # self.query_dict = {}
        # self.all_dict = {}

    @classmethod
    def init_unit(given_unit):

        # dictionary to transform argument options into coverm units names 
        unit_dict ={'rpkm': 'RPKM', 'tpm': 'TPM', 'tm':'Trimmed Mean'} 
        units = unit_dict[given_unit] 

        Eggnog_sample.option_unit = units
        
        # if the parsed unit is TPM, we will use RPKM for calculations and then transform it to TPM
        if given_unit == 'tpm': 
            Eggnog_sample.calc_unit = 'RPKM'
        else:
            Eggnog_sample.calc_unit = units
    
    @classmethod
    def init_sample_list(sample_list):
        Eggnog_sample.sample_list = sample_list


    def define_samplename(self):

        """
        Instance method to extract samplename from the filename
        """

        basename = os.path.basename(self.filename)
        self.samplename = re.sub(r'.emapper.*', '', basename) # filename needs to be of format: sample.emmapper.annotations
        
        return self.samplename
    
    def load_sample(self, coverm_dict, og_dict):

        """
        Instance method to load all rows from a file as eggnog_orf instances and store them as a 
        eggnog_sample instance
        """
        
        self.rows=[]
        with open(self.filename,"r") as file:
            for line in file:
                if not line.startswith("##"): # skip headers

                    if line.startswith("#"):
                        self.header = line.strip() #.split('\t') # column names as string, eliminar #??
                    else:  # read lines
                        items = line.strip().split('\t')
                        eggnog_ogs = items[4]
                        kegg_ko = get_ko_list(items[11])

                        if self.remove_euk and '@2759' in eggnog_ogs: # remove euk. Esto puede ser un argumento 
                            continue

                        else:

                            query = items[0]
                            og, kingdom = find_basal(eggnog_ogs)
                            description = items[7]
                            preferred_name = items[8]
                            kegg_pathway = items[12]
                            contig = re.sub(r'_[0-9]*$', '', query)
                            abundance = coverm_dict[contig][Eggnog_sample.calc_unit]


                            eggnog_orf = Eggnog_orf(query, og, kingdom, description, preferred_name, kegg_ko, kegg_pathway, contig, abundance)
                            self.rows.append(eggnog_orf)

                            self.add_og_abundance(og, abundance, description, kingdom, og_dict)
                            self.add_ko_abundance(kegg_ko, abundance)

                            self.mapped += abundance # estoy sumando por query !!
                            self.total += abundance # para añadirlo al total y poder calcular bien la abundancia relativa no?
                            
                            
                            # if we want to add up abundance for each ko :
                            # self.mapped_ko += abundance * len(kegg_ko)
                            # self.total_ko += abundance * len(kegg_ko)

                            # self.query_list.append(query)
                            # self.query_dict[query] = kegg_ko
                            # self.all_dict[query] = {}
                            # self.all_dict[query]['ko'] = kegg_ko
                            # self.all_dict[query]['cog'] = find_basal(eggnog_ogs)
                            # self.all_dict[query]['description'] = description
        
    
    def add_og_abundance(self, og, abundance, des, king, og_dict):

        """
        calculate the total abundance for each og in a sample
        """
        
        # try adding the abundance of the og to the og_abundance_sample dictionary
        try: 
            self.og_abundance[og] += float(abundance)

        # if the og is not already in the dictionary :    
        except:
            self.og_abundance[og] = float(abundance)
            og_dict[og] = {}
            for sample in Eggnog_sample.sample_list:
                og_dict[og][sample] = 0
            og_dict[og]['description'] = des # save the og description 
            og_dict[og]['king'] = king # save the og kingdom
    
        #return og_dict #, self.og_abundance
    
    def add_ko_abundance(self, kos, abundance):
        
        """
        add abundance to each ko on the list
        """
        # divide abundance by number of kos annotated in the contig
        abun = float(abundance)/len(kos)

        for ko_id in kos:
            # initialize dictionary for each ko
            if not ko_id in self.ko_abundance.keys():
                self.ko_abundance[ko_id] = 0
            
            # add ko abundance 
            self.ko_abundance[ko_id] += abun # dividing by ko number
            # adding up
            # ko_abundance_sample[ko_id] += float(abundance) 
    
        #return self.ko_abundance

    def calculate_og_abundance(self, og_dict):

        """
        
        """
       
        if Eggnog_sample.option_unit == 'TPM':

            for og in self.og_abundance.keys():
                    
                og_dict[og][self.samplename] = float(self.og_abundance[og])/(self.total*10**6)
                
            og_dict = check_unmapped(og_dict, Eggnog_sample.sample_list)
            og_dict['UNMAPPED'][self.samplename] = 1000000 - float(self.mapped)/(self.total*10**6)

        else: 

            for og in self.og_abundance.keys():
                    
                og_dict[og][self.samplename] = float(self.og_abundance[og])/self.total
            
            og_dict = check_unmapped(og_dict, Eggnog_sample.sample_list)
            og_dict['UNMAPPED'][self.samplename] = 1 - float(self.mapped)/self.total
        
        # CHECK UP
        if og_dict[og]['UNMAPPED'] < 0 :
            print('ERROR CALCULATING ABUNDANCES')
    
            # except:
            #     og_dict[og] = {}
            #     og_dict[og]['kingdom'] = og_dict[og]  
            #     og_dict[og]['description'] = og_description[og]
            #     for sample_id in sample_list:
            #         rel_og_abun[og][sample_id] = 0

            #rel_og_abun[og][samplename] = og_abundance_sample[og]/total_og_abun*10**6     

        return og_dict
    
    def calculate_ko_abundance(self, ko_dict, kos_legend):

        """
        Calculate every ko total abundance in a sample. 
        This function relies on contigs and orf being sorted.
        Takes rpkm value as abundance
        """
        # CHECK WHEN A SEED CORRESPOND TO SEVERAL KOS
        # more_kos = []
        # total_kos = 0 

        # when adding up
        #ko_dict['UNMAPPED'][self.samplename] = (self.total_ko) - float(self.mapped_ko)
        
        if Eggnog_sample.option_unit == 'TPM':

            for ko_id in self.ko_abundance.keys():

                self.global_ko_list.append(ko_id)

                if not ko_id in ko_dict.keys():
                    ko_dict[ko_id] = {}
                    try: 
                        ko_dict[ko_id]['description'] = kos_legend[ko_id]['description']
                    except:
                        ko_dict[ko_id]['description'] = '@'
                    
                    for sample_id in Eggnog_sample.sample_list:
                        ko_dict[ko_id][sample_id]=0
            
                ko_dict[ko_id][self.samplename] =self.ko_abundance[ko_id]/(self.total*10**6)
            
            ko_dict = check_unmapped(ko_dict, Eggnog_sample.sample_list)
            ko_dict['UNMAPPED'][self.samplename] = 1000000 - float(self.mapped)/(float(self.total)*10**6)

        else:

            for ko_id in self.ko_abundance.keys():

                self.global_ko_list.append(ko_id)

                if not ko_id in ko_dict.keys():
                    ko_dict[ko_id] = {}
                    try: 
                        ko_dict[ko_id]['description'] = kos_legend[ko_id]['description']
                    except:
                        ko_dict[ko_id]['description'] = '@'
                    
                    for sample_id in Eggnog_sample.sample_list:
                        ko_dict[ko_id][sample_id]=0
            
                ko_dict[ko_id][self.samplename] = self.ko_abundance[ko_id]/self.total
            
            ko_dict = check_unmapped(ko_dict, Eggnog_sample.sample_list)
            ko_dict['UNMAPPED'][self.samplename] = 1 - float(self.mapped)/float(self.total)


        # CHECK WHEN A SEED CORRESPOND TO SEVERAL KOS
        # total = len(more_kos)
        # two = more_kos.count(2)
        # three = more_kos.count(3)
        # rest = total - two - three 
        # print(total_kos, round((total/total_kos)*100,2), total, two, three, rest)


        return ko_dict

            

    def calculate_KEGG_pathway_completeness(self, path_coverage, KEGG_dict):

        """
        Function to calculate KEGG pathway completeness of each sample 
        """

        kegg_cov_dict = {}

        for kegg_p, annotation in KEGG_dict.items():
            pathway_description = annotation[0]
            kegg_number = annotation[1]
            kegg_cov_dict[kegg_p] = 0   
            for ko in annotation[2]:
                ko_id = ko['KO']
                # symbol = ko['symbol']
                # description = ko['description']
                if ko_id in self.global_KOs_list:
                    kegg_cov_dict[kegg_p] +=1
            
            if kegg_cov_dict[kegg_p] != 0:
                coverage = kegg_cov_dict[kegg_p]/kegg_number
                if coverage > 0.1: # cut-off proporcion
                    if not kegg_p in path_coverage.keys():
                        path_coverage[kegg_p] = {}
                        path_coverage[kegg_p]['description'] = pathway_description
                        for sample_id in Eggnog_sample.sample_list:
                            path_coverage[kegg_p][sample_id] = 0
                    
                    path_coverage[kegg_p][self.samplename] = coverage
                    
        return path_coverage


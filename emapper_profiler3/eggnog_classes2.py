#!/usr/bin/python
# -*- coding: utf-8 -*-

import re
import os
from ko_functions2 import get_ko_list, find_basal, check_unmapped

            
class Eggnog_sample(object):

    """
    Class to perform all functional profiling functions for each sample 
    annotated with eggNOG-mapper. 
    """

    option_unit = None  # the option unit as given in the arguments
    calc_unit = None    # the unit we use for calculations
    sample_list = None
    

    def __init__(
            self, filename: str, total_sample:float, 
            samplename = None, remove_euk = False, remove_virus = False) -> None:

        """
        Init instance method to initialize the class
        """

        self.filename = filename
        self.remove_euk = remove_euk
        self.remove_virus = remove_virus
        self.total = total_sample
        if samplename == None :
            self.samplename = self.define_samplename()
        else: 
            self.samplename = samplename
        
        # Initializing instance variables 
        self.og_abundance = {}
        self.ko_abundance = {}
        self.mapped_og = 0
        self.mapped_ko = 0
        self.total_ko = total_sample
        self.global_ko_list = []
        self.contig_ko = {}
        self.contig_og = {}


    @classmethod
    def init_unit(
            cls, given_unit):

        """
        Class method to initialize the units to work with
        """

        # Transform argument options into coverm units names 
        unit_dict ={'rpkm': 'RPKM', 'tpm': 'TPM', 'tm':'Trimmed_Mean'} 
        units = unit_dict[given_unit] 

        Eggnog_sample.option_unit = units
        
        # If the parsed unit is TPM, we will use RPKM for calculations and then transform it to TPM
        if given_unit == 'tpm': 
            Eggnog_sample.calc_unit = 'RPKM'
        else:
            Eggnog_sample.calc_unit = units
    

    @classmethod
    def init_sample_list(
            cls, sample_list):
        
        """
        Class method to initialize the sample list class variable
        """

        Eggnog_sample.sample_list = sample_list


    def define_samplename(
            self):

        """
        Instance method to extract samplename from the filename
        """

        # Filename needs to be of format: 'sample.emmapper.annotations'
        basename = os.path.basename(self.filename)
        self.samplename = re.sub(r'.emapper.*', '', basename) 
        
        return self.samplename
    

    def load_sample(
            self, orf_dict, coverm_dict, 
            og_dict, og_dict2, ko_dict2, kos_legend):

        """
        Instance method to add the abundances of each annotation 
        as the sample file is being read. 
        """
        
        with open(self.filename,"r") as file:
            
            for line in file:

                if not line.startswith("##"):  # Skip headers

                    if line.startswith("#"):
                        self.header = line.strip()  # Column names as string
                    
                    else:  # Read lines
                        items = line.strip().split('\t')
                        
                        eggnog_ogs = items[4] 
                        raw_ko = items[11]
                        
                        # Skip eukaryotes lines when remove_euk flag is on
                        if self.remove_euk == True and '@2759' in eggnog_ogs: 
                            continue
                        
                        # Skip viruses lines when remove_virus flag is on 
                        elif self.remove_virus == True and '@10239' in eggnog_ogs:
                            continue

                        else:
                            query = items[0]
                            description = items[7]
                            contig = re.sub(r'_[0-9]*$', '', query)
                            abundance = orf_dict[query]['abundance']
                            tpm = coverm_dict[contig]['TPM']
                            
                            try:
                                og, kingdom = find_basal(eggnog_ogs)
                                og_dict = self.add_og_abundance(og, abundance, description, kingdom, og_dict)
                                og_dict2 = self.calc_contig_og_abundance(og, contig, tpm, description, kingdom, og_dict2)
                                self.mapped_og += float(abundance)  # Add abundance to total mapped og
                            
                            # Exception happens because the OG belongs to a virus. 
                            except:
                                og = eggnog_ogs.split('@')[0]
                                if 'viri' in eggnog_ogs:
                                    if self.remove_virus == False:
                                        kingdom = 'Viruses'
                                        og_dict = self.add_og_abundance(og, abundance, description, kingdom, og_dict)
                                        og_dict2 = self.calc_contig_og_abundance(og, contig, tpm, description, kingdom, og_dict2)
                                        self.mapped_og += float(abundance)  # Add abundance to total mapped OG
                                    else:  # remove_virus == True, don't annotate the KO
                                        raw_ko = '-' 
                                
                                else:  # This else usually does not happen 
                                    kingdom = eggnog_ogs.split('|')[1]  # = 'Other'
                                    og_dict = self.add_og_abundance(og, abundance, description, kingdom, og_dict)
                                    og_dict2 = self.calc_contig_og_abundance(og, contig, tpm, description, kingdom, og_dict2)
                                    self.mapped_og += float(abundance)  # Add abundance to total mapped OG

                            if raw_ko != '-':
                                kegg_ko = get_ko_list(items[11])
                                self.add_ko_abundance(kegg_ko, abundance)
                                ko_dict2 = self.calc_contig_ko_abundance(contig, kegg_ko, tpm, ko_dict2, kos_legend)
                                self.mapped_ko += (float(abundance)*len(kegg_ko))  # Add abundance to total mapped KOs
                                self.total_ko += (float(abundance)*(len(kegg_ko)-1))  # Add abundance to total KO functions minus the one that is already counted


        return og_dict, og_dict2, ko_dict2

    
    def add_og_abundance(
            self, og, abundance, des, king, og_dict):

        """
        Instance method to calculate the total abundance for each Orthologous Group (OG) in a sample
        It takes as arguments the OG, its abundance in the ORF, the OG's description and kingdom 
        and the global dictionary for all OGs per sample
        """

        if og not in self.og_abundance.keys():
            self.og_abundance[og] = 0
            
        # For each OG, the abundance found in each Open Reading Frame is accumulated   
        self.og_abundance[og] += float(abundance)

        # The global dictionary is initialize in this step with the OG's descriptions and kingdoms
        if og not in og_dict.keys():
            og_dict[og] = {}
            og_dict[og]['description'] = des  
            og_dict[og]['kingdom'] = king  

            # Initialize OG to 0 in all samples    
            for sample in Eggnog_sample.sample_list:
                og_dict[og][sample] = 0
    
        return og_dict
    

    def add_ko_abundance(
            self, kos, abundance):
        
        """
        Instance method to calculate the total abundance for each KEGG Ortholog (KO) in a sample.
        It takes as arguments a list of KOs annotated for an ORF and the abundance of said ORF. 
        """

        # For each KO, the abundance found in the Open Reading Frame is accumulated   
        for ko_id in kos:
            
            if not ko_id in self.ko_abundance.keys():
                self.ko_abundance[ko_id] = 0
            
            self.ko_abundance[ko_id] += abundance 
    

    def calculate_og_abundance(
            self, og_dict:dict):

        """
        Instance method to calculate the relative abundance of each OG in each sample. 
        For each OG, the total abundance acumulated is divided by the total ORFs abundance. 
        """

        # If the relative abundance is calculated as TPM, all OG's abundances of a sample sums 1 million
        if Eggnog_sample.option_unit == 'TPM':

            for og in self.og_abundance.keys():
                    
                og_dict[og][self.samplename] = (self.og_abundance[og]/self.total)*10**6
                
            og_dict = check_unmapped(og_dict, Eggnog_sample.sample_list)
            
            # The unmapped proportion is calculated as 1 million minus the mapped proportion in TPM
            og_dict['UNMAPPED'][self.samplename] = float(1000000) - (self.mapped_og/self.total)*10**6
            og_dict['UNMAPPED']['kingdom'] = '@'

        
        # In any other case, all OG's abundances of a sample sums 1
        else: 

            for og in self.og_abundance.keys():
                    
                og_dict[og][self.samplename] = self.og_abundance[og]/self.total
            
            og_dict = check_unmapped(og_dict, Eggnog_sample.sample_list)
            
            # The unmapped proportion is calculated as 1 minus the mapped proportion
            og_dict['UNMAPPED'][self.samplename] = 1 - self.mapped_og/self.total
            og_dict['UNMAPPED']['kingdom'] = '@'
        

        return og_dict
    
    def calculate_ko_abundance(
            self, ko_dict, kos_legend):

        """
        Instance method to calculate the relative abundance of each KO in each sample. 
        For each KO, the total abundance acumulated is divided by the total functions (KO and unknowns) abundances. 
        """
        
        # If the relative abundance is calculated as TPM, all KO's abundances of a sample sums 1 million
        if Eggnog_sample.option_unit == 'TPM':

            for ko_id in self.ko_abundance.keys():

                self.global_ko_list.append(ko_id)

                # Each KO description and symbol is annotated considering the given kos_legend dictionary
                if not ko_id in ko_dict.keys():
                    ko_dict[ko_id] = {}
                    try: 
                        ko_dict[ko_id]['description'] = kos_legend[ko_id]['description']
                    except:
                        ko_dict[ko_id]['description'] = '@'
                    try: 
                        ko_dict[ko_id]['symbol'] = str(kos_legend[ko_id]['symbol'])
                    except:
                        ko_dict[ko_id]['symbol'] = '@'
                    
                    # Initialize KO to 0 in all samples    
                    for sample_id in Eggnog_sample.sample_list:
                        ko_dict[ko_id][sample_id]=0
            
                ko_dict[ko_id][self.samplename] = (self.ko_abundance[ko_id]/self.total_ko)*10**6

            ko_dict = check_unmapped(ko_dict, Eggnog_sample.sample_list)
            
            # The unmapped proportion is calculated as 1 million minus the mapped proportion in TPM
            ko_dict['UNMAPPED'][self.samplename] = float(1000000) - (self.mapped_ko/self.total_ko)*10**6 
            ko_dict['UNMAPPED']['symbol'] = '@'

        
        # In any other case, all OG's abundances of a sample sums 1
        else:

            for ko_id in self.ko_abundance.keys():

                self.global_ko_list.append(ko_id)

                # Each KO description and symbol is annotated considering the given kos_legend
                if not ko_id in ko_dict.keys():
                    ko_dict[ko_id] = {}
                    try: 
                        ko_dict[ko_id]['description'] = kos_legend[ko_id]['description']
                    except:
                        ko_dict[ko_id]['description'] = '@'
                    try: 
                        ko_dict[ko_id]['symbol'] = str(kos_legend[ko_id]['symbol'])
                    except:
                        ko_dict[ko_id]['symbol'] = '@'
                    
                    for sample_id in Eggnog_sample.sample_list:
                        ko_dict[ko_id][sample_id]=0
            
                ko_dict[ko_id][self.samplename] = self.ko_abundance[ko_id]/self.total
            
            ko_dict = check_unmapped(ko_dict, Eggnog_sample.sample_list)

            # The unmapped proportion is calculated as 1 minus the mapped proportion
            ko_dict['UNMAPPED'][self.samplename] = 1 - self.mapped_ko/self.total
            ko_dict['UNMAPPED']['symbol'] = '@'


        return ko_dict


    def calc_contig_og_abundance(
            self, og, contig, TPM, des, king, og_dict2):

        """
        Instance method to calculate the abundance (in TPM) per contig of each OG in each sample. 
        If a OG is present in a contig, the abundance of the contig is accumulated (only once).
        """
   
        # Each KO description and kingdom is annotated in the global dictionary
        if og not in og_dict2.keys():
            og_dict2[og] = {}
            og_dict2[og]['description'] = des  
            og_dict2[og]['kingdom'] = king 

            # Initialize to 0 in all samples
            for sample in Eggnog_sample.sample_list:
                og_dict2[og][sample] = 0 

        # Each KO description and kingdom is annotated in the sample's dictionary    
        if og not in self.contig_og.keys():
            self.contig_og[og] = {}
            self.contig_og[og]['description'] = des 
            self.contig_og[og]['kingdom'] = king 

            self.contig_og[og]['total'] = 0 

        # If the OG has not been found in the contig yet
        if contig not in self.contig_og[og].keys():
            og_dict2[og][self.samplename] += float(TPM)  # Accumulate to the total abundance in the global dictionary   
            
            self.contig_og[og][contig] = float(TPM)      # Accumulate to the total abundance in the sample's dictionary 
            self.contig_og[og]['total'] += float(TPM)    # Add the contig key with the respective abundance in the sample's dictionary

        return og_dict2
    

    def calc_contig_ko_abundance(
            self, contig, kos, TPM, ko_dict2, kos_legend):

        """
        Instance method to calculate the abundance (in TPM) per contig of each KO in each sample. 
        If a KO is present in a contig, the abundance of the contig is accumulated (only once).
        """
        
        for ko in kos:

            # Each KO description and symbol is annotated considering the given kos_legend in the global dictionary
            if not ko in ko_dict2.keys():
                ko_dict2[ko] = {}
                try: 
                    ko_dict2[ko]['description'] = kos_legend[ko]['description'] 
                except:
                    ko_dict2[ko]['description'] = '@'
                try: 
                    ko_dict2[ko]['symbol'] = str(kos_legend[ko]['symbol'])
                except:
                    ko_dict2[ko]['symbol'] = '@'
                
                # Initialize KO to 0 in all samples    
                for sample_id in Eggnog_sample.sample_list:
                    ko_dict2[ko][sample_id]=0
            
            # Each KO description and symbol is annotated considering the given kos_legend in the sample's dictionary
            if ko not in self.contig_ko.keys():
                self.contig_ko[ko] = {}
                try: 
                    self.contig_ko[ko]['description'] = kos_legend[ko]['description']
                except:
                    self.contig_ko[ko]['description'] = '@'
                try: 
                    self.contig_ko[ko]['symbol'] = kos_legend[ko]['symbol']
                except:
                    self.contig_ko[ko]['symbol'] = '@'

                self.contig_ko[ko]['total'] = 0 

            # If the KO has not been found in the contig yet
            if contig not in self.contig_ko[ko].keys():
                ko_dict2[ko][self.samplename] += float(TPM)  # Accumulate to the total abundance in the global dictionary   
                
                self.contig_ko[ko]['total'] += float(TPM)  # Accumulate to the total abundance in the sample's dictionary 
                self.contig_ko[ko][contig] = float(TPM)    # Add the contig key with the respective abundance in the sample's dictionary
  
        return ko_dict2
            

    def calculate_KEGG_pathway_completeness(self, path_coverage, KEGG_dict):

        """
        Instance method to calculate the KEGG pathway completeness of each sample. 
        Each KO found in the sample is taken into account in this step. 
        """

        kegg_cov_dict = {}

        for kegg_p, annotation in KEGG_dict.items():
            kegg_id = str('ko'+kegg_p)
            pathway_description = annotation[0]
            kegg_number = annotation[1]
            kegg_cov_dict[kegg_id] = 0   
            for ko in annotation[2]:
                ko_id = ko['KO']
                if ko_id in self.global_ko_list:
                    kegg_cov_dict[kegg_id] +=1
            
            if kegg_cov_dict[kegg_id] != 0:
                coverage = kegg_cov_dict[kegg_id]/kegg_number
                if coverage > 0.1:  # cut-off proportion
                    if not kegg_id in path_coverage.keys():
                        path_coverage[kegg_id] = {}
                        path_coverage[kegg_id]['description'] = pathway_description
                        for sample_id in Eggnog_sample.sample_list:
                            path_coverage[kegg_id][sample_id] = 0
                    
                    path_coverage[kegg_id][self.samplename] = coverage
                    
        return path_coverage


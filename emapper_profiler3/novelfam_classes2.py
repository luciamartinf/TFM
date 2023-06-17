#!/usr/bin/python
# -*- coding: utf-8 -*-

from dataclasses import dataclass, field
from ko_functions2 import check_unmapped, check_key
import numpy as np
import re
import os

# @dataclass
# class NovelFam_orf(object):

#     """
#     Dataclass to store novel families' attributes corresponding to one orf
#     """

#     query: str
#     #target: str
#     #evalue: str # = field(repr=False) #1.25e-07?? 
#     #score: float # = field(repr=False)
#     novel_fam: str
#     contig: str #= field(init=False) 

#     def __str__(self):

#         """
#         Instance method to define the print format of the instance
#         """

#         row_list = [self.query, self.novel_fam, self.contig]
#         s = '\t'.join([str(item) for item in row_list])
        
#         return s



class NovelFam_sample(object):

    """
    Class to store all rows from an eggnog-mapper file. 
    Each row is a eggnog_orf instance
    """

    option_unit = None # the option unit as given in the arguments
    calc_unit = None # the unit we use for calculations
    sample_list = None

    def __init__(self, filename: str, total_sample:float, samplename = None) -> None:
        self.filename = filename
        if samplename == None :
            self.samplename = self.define_samplename() # function to extract samplename from filename
        else: 
            self.samplename = samplename
        # self.query_list = []
        # self.query_dict = {}
        self.total = total_sample
        self.mapped = 0
        self.nf_abundance = {}
        self.contig_nf = {}

    @classmethod
    def init_unit(cls, given_unit):

        # dictionary to transform argument options into coverm units names 
        unit_dict ={'rpkm': 'RPKM', 'tpm': 'TPM', 'tm':'Trimmed Mean'} 
        units = unit_dict[given_unit] 

        NovelFam_sample.option_unit = units
        
        # if the parsed unit is TPM, we will use RPKM for calculations and then transform it to TPM
        if given_unit == 'tpm': 
            NovelFam_sample.calc_unit = 'RPKM'
        else:
            NovelFam_sample.calc_unit = units
    
    @classmethod
    def init_sample_list(cls, sample_list):
        NovelFam_sample.sample_list = sample_list

    def define_samplename(self):

        """
        Instance method to extract samplename from the filename
        """

        basename = os.path.basename(self.filename)
        self.samplename = re.sub(r'.emapper.*', '', basename) # filename needs to be of format: sample.emmapper.annotations
        
        return self.samplename
    
    def load_sample(self, orf_dict, coverm_dict, nf_dict2):

        """
        Instance method to load all rows from a file as eggnog_orf instances and store them as a 
        eggnog_sample instance
        """
        
        with open(self.filename,"r") as file:
            for line in file:
                if not line.startswith("##"): # skip headers

                    if line.startswith("#"):
                        self.header = line.strip() #.split('\t') # column names as string, eliminar #??
                    else:  # read lines
                        items = line.strip().split('\t')
                        novel_fam = items[4]
                        query = items[0]
                        abundance = orf_dict[query]['abundance']
                        contig = re.sub(r'_[0-9]*$', '', query)
                        tpm = coverm_dict[contig]['TPM']

                        self.add_nf_abundance(novel_fam, abundance)
                        self.mapped += abundance

                        nf_dict2 = self.calc_contig_nf_abundance(novel_fam, contig, tpm, nf_dict2)

                        # self.query_list.append(query)
                        # self.query_dict[novel_fam] = query
        
        return nf_dict2

    def add_nf_abundance(self, nf, abundance):
        
        """
        add abundance to each novel family id
        """

        if not nf in self.nf_abundance.keys():
            self.nf_abundance[nf] = abundance
        else:
            self.nf_abundance[nf] += abundance
    

    def calculate_nf_abundance(self, nf_dict:dict):
       
        if NovelFam_sample.option_unit == 'TPM':

            for nf in self.nf_abundance.keys():
                nf_dict = check_key(nf_dict, nf, NovelFam_sample.sample_list) 
                nf_dict[nf][self.samplename] = (self.nf_abundance[nf]/self.total)*10**6
                
            nf_dict = check_unmapped(nf_dict, NovelFam_sample.sample_list, des=False)
            #tpm_mapped_og = (self.mapped_og/self.total)*10**6 # MAPPED IN TPM
            #og_dict['UNMAPPED'][self.samplename] = 1000000 * (self.total - self.mapped_og)/self.total # es lo mismo
            nf_dict['UNMAPPED'][self.samplename] = float(1000000) - (self.mapped/self.total)*10**6

        else: 

            for nf in self.nf_abundance.keys():

                nf_dict = check_key(nf_dict, nf, NovelFam_sample.sample_list)    
                nf_dict[nf][self.samplename] = self.nf_abundance[nf]/self.total
            
            nf_dict = check_unmapped(nf_dict, NovelFam_sample.sample_list, des=False)
            nf_dict['UNMAPPED'][self.samplename] = 1 - self.mapped/self.total

        
        return nf_dict
    

    def calc_contig_nf_abundance(self, nf, contig, TPM, nf_dict2):
   
        # if the og is not in the global og_dict initialize it in the dictionary
        if nf not in nf_dict2.keys():
            nf_dict2[nf] = {}
            for sample in NovelFam_sample.sample_list:
                nf_dict2[nf][sample] = 0 # initialize to 0 in all samples

        # check if the nf is initialize in the sample's dictionary    
        if nf not in self.contig_nf.keys():
            self.contig_nf[nf] = {}
            self.contig_nf[nf]['total'] = 0 # initialize total to 0

        # if the contig hasn't been sum
        if contig not in self.contig_nf[nf].keys():
            nf_dict2[nf][self.samplename] += float(TPM)  # add to global dict   
            
            self.contig_nf[nf][contig] = float(TPM)    # add contig abundance in sample dict
            self.contig_nf[nf]['total'] += float(TPM)  # add to sample dict

        return nf_dict2
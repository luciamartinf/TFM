#!/usr/bin/python
# -*- coding: utf-8 -*-

from dataclasses import dataclass, field
from functions import check_unmapped, check_key
import numpy as np
import re
import os


class NovelFam_sample(object):

    """
    Class to perform all functional profiling functions for each sample 
    annotated with eggNOG-mapper to the Novel Gene Families database. 
    """

    option_unit = None # the option unit as given in the arguments
    calc_unit = None # the unit we use for calculations
    sample_list = None

    def __init__(
            self, filename: str, total_sample:float, samplename = None) -> None:
        
        """
        Init instance method to initialize the class
        """

        self.filename = filename
        self.total = total_sample
        if samplename == None :
            self.samplename = self.define_samplename() # function to extract samplename from filename
        else: 
            self.samplename = samplename

        # Initializing instance variables 
        self.mapped = 0
        self.nf_abundance = {}

    @classmethod
    def init_unit(
        cls, given_unit):

        """
        Class method to initialize the units to work with
        """

        # Transform argument options into coverm units names 
        unit_dict ={'rpkm': 'RPKM', 'tpm': 'TPM', 'tm':'Trimmed Mean'} 
        units = unit_dict[given_unit] 

        NovelFam_sample.option_unit = units
        
        # If the parsed unit is TPM, we will use RPKM for calculations and then transform it to TPM
        if given_unit == 'tpm': 
            NovelFam_sample.calc_unit = 'RPKM'
        else:
            NovelFam_sample.calc_unit = units
    
    @classmethod
    def init_sample_list(
        cls, sample_list):
        
        """
        Class method to initialize the sample list class variable
        """

        NovelFam_sample.sample_list = sample_list

    def define_samplename(self):

        """
        Instance method to extract samplename from the filename
        """
        
        # Filename needs to be of format: 'sample.emmapper.annotations'
        basename = os.path.basename(self.filename)
        self.samplename = re.sub(r'.emapper.*', '', basename) # filename needs to be of format: sample.emmapper.annotations
        
        return self.samplename
    
    def load_sample(
            self, orf_dict, nf_dict3):

        """
        Instance method to add the abundances of each annotation 
        as the sample file is being read. 
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

                        self.add_nf_abundance(novel_fam, abundance)
                        self.mapped += abundance

                        nf_dict3 = self.calc3_nf_abundance(novel_fam, abundance, nf_dict3)
        
        return nf_dict3

    def add_nf_abundance(
            self, nf, abundance):
        
        """
        Instance method to calculate the total abundance for each Novel Family (NF) in a sample
        """

        if not nf in self.nf_abundance.keys():
            self.nf_abundance[nf] = abundance
        else:
            self.nf_abundance[nf] += abundance
    

    def calculate_nf_abundance(self, nf_dict:dict):
        
        """
        Instance method to calculate the relative abundance of each NF in each sample. 
        For each NF, the total abundance acumulated is divided by the total ORFs abundance. 
        """

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
    
    
    def calc3_nf_abundance(self, nf, abun, nf_dict3):
   
        # if the og is not in the global og_dict initialize it in the dictionary
        if nf not in nf_dict3.keys():
            nf_dict3[nf] = {}
            for sample in NovelFam_sample.sample_list:
                nf_dict3[nf][sample] = 0 # initialize to 0 in all samples

        nf_dict3[nf][self.samplename] += float(abun)  # add to global dict   
          

        return nf_dict3
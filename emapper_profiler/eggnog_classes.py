#!/usr/bin/python
# -*- coding: utf-8 -*-

from dataclasses import dataclass
import re
import os
from ko_functions import find_basal

@dataclass
class Eggnog_orf(object):

    """
    Dataclass to store eggnog's attributes corresponding to one orf
    """

    query: str
    seed_ortholog: str
    evalue: str # = field(repr=False) #1.25e-07?? 
    score: float # = field(repr=False)
    eggnog_ogs: str
    max_annot_lvl: str
    #COG_category: str
    description: str # = field(repr=False)
    preferred_name:	str
    #GOs: str	
    #EC: str
    kegg_ko: str
    kegg_pathway: str
    #KEGG_Module: str	
    #KEGG_Reaction: str
    #KEGG_rclass: str
    #BRITE: str
    #KEGG_TC: str 	
    #CAZy: str 
    #BiGG_Reaction: str	
    #PFAMs: str
    contig: str #= field(init=False) 

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

    def __init__(self, filename: str, samplename = None, remove_euk = False) -> None:
        self.filename = filename
        self.remove_euk = remove_euk
        if samplename == None :
            self.samplename = self.define_samplename() # function to extract samplename from filename
        else: 
            self.samplename = samplename
        self.query_list = []
        self.query_dict = {}
        self.all_dict = {}

    def define_samplename(self):

        """
        Instance method to extract samplename from the filename
        """

        basename = os.path.basename(self.filename)
        self.samplename = re.sub(r'.emapper.*', '', basename) # filename needs to be of format: sample.emmapper.annotations
        
        return self.samplename
    
    def load_sample(self):

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
                        kegg_ko = items[11]

                        if self.remove_euk and '@2759' in eggnog_ogs: # remove euk. Esto puede ser un argumento 
                            continue

                        else:
                            query = items[0]
                            seed_ortholog = items[1]
                            evalue = items[2]
                            score = float(items[3])
                            max_annot_lvl = items[5]
                            description = items[7]
                            preferred_name = items[8]
                            kegg_pathway = items[12]
                            contig = re.sub(r'_[0-9]*$', '', query)

                                # line_list = line.split('\t')
                                # args_list = line_list[0:6] + line_list[7:9] + line_list[11:13]
                                # args = ' ,'.join([str(item) for item in args_list])
                                # eggnog_orf = Eggnog_orf(args)

                            eggnog_orf = Eggnog_orf(query, seed_ortholog, evalue, score, eggnog_ogs, max_annot_lvl, description, preferred_name, kegg_ko, kegg_pathway, contig)
                            self.rows.append(eggnog_orf)
                            self.query_list.append(query)
                            self.query_dict[query] = kegg_ko
                            self.all_dict[query] = {}
                            self.all_dict[query]['ko'] = kegg_ko
                            self.all_dict[query]['cog'] = find_basal(eggnog_ogs)
                            self.all_dict[query]['description'] = description

    def calculate_contig_abundance(coverm_sample, unit):

        i = 0

        for eggnog_row in self.rows:

            while coverm_sample.rows[i].contig != eggnog_row.contig:
                i += 1
            
            eggnog_row.abundance = getattr(coverm_sample.rows[i], unit)

# eggnog_file = "0505-0101.emapper.annotations"

# sample_05050101 = Eggnog_sample(eggnog_file)
# sample_05050101.load_sample()

# print(sample_05050101.header)
# for row in sample_05050101.rows:
#      print(row.query)
#      print(row.contig)
#      print(row.max_annot_lvl) # necesito saber con que parte de la notacion taxonomica me quedo
#      print(row) # lo hace ya separado por tab por el __str__ method. 
#      break


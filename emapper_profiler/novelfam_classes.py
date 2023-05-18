#!/usr/bin/python
# -*- coding: utf-8 -*-

from dataclasses import dataclass, field
import numpy as np
import re
import os

@dataclass
class NovelFam_orf(object):

    """
    Dataclass to store novel families' attributes corresponding to one orf
    """

    query: str
    #target: str
    #evalue: str # = field(repr=False) #1.25e-07?? 
    #score: float # = field(repr=False)
    novel_fam: str
    contig: str #= field(init=False) 

    def __str__(self):

        """
        Instance method to define the print format of the instance
        """

        row_list = [self.query, self.novel_fam, self.contig]
        s = '\t'.join([str(item) for item in row_list])
        
        return s



class NovelFam_sample(object):

    """
    Class to store all rows from an eggnog-mapper file. 
    Each row is a eggnog_orf instance
    """

    def __init__(self, filename: str, samplename = None) -> None:
        self.filename = filename
        if samplename == None :
            self.samplename = self.define_samplename() # function to extract samplename from filename
        else: 
            self.samplename = samplename
        self.query_list = []
        self.query_dict = {}

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
                        novel_fam = items[4]
                        query = items[0]
                        contig = re.sub(r'_[0-9]*$', '', query)

                                # line_list = line.split('\t')
                                # args_list = line_list[0:6] + line_list[7:9] + line_list[11:13]
                                # args = ' ,'.join([str(item) for item in args_list])
                                # eggnog_orf = Eggnog_orf(args)

                        novelfam_orf = NovelFam_orf(query, novel_fam, contig)
                        self.rows.append(novelfam_orf)
                        self.query_list.append(query)
                        self.query_dict[novel_fam] = query


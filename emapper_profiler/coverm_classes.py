#!/usr/bin/python
# -*- coding: utf-8 -*-

from dataclasses import dataclass, field
import numpy as np
import re
import os


@dataclass
class Coverm_contig(object):

    """
    Dataclass to store coverm's attributes corresponding to one contig
    """

    contig: str
    contig_length: int
    read_count: int
    reads_per_base: float
    #mean: float
    trimmed_mean: float
    tpm: float
    rpkm: float
    

    def __str__(self):

        """
        Instance method to define the print format of the instance
        """

        row_list = [self.contig, self.contig_length, self.read_count, self.reads_per_base, self.trimmed_mean, self.tpm, self.rpkm]
        s = '\t'.join([str(item) for item in row_list])
        
        return s



class Coverm_sample(object):

    """
    Class to store all rows from a coverm file. 
    Each row is a coverm_contig instance
    """

    def __init__(self, filename: str, samplename = None) -> None:
        self.filename = filename
        self.samplename = samplename # function to extract samplename from filename
        if samplename == None :
            self.samplename = self.define_samplename() # function to extract samplename from filename
        else: 
            self.samplename = samplename
        self.total_rpkm = 0
        self.total_tm = 0

    def define_samplename(self):

        """
        Instance method to extract samplename from the filename
        """
        
        basename = os.path.basename(self.filename)
        self.samplename = re.sub(r'_coverage_values', '', basename) # filename needs to be of format: sample_coverage_values
        
        return self.samplename
    
    def define_clean_header(self, line):
        
        """
        Simplify column names in coverage file
        """
        
        clean_header = []
        header_list = line.split('\t')
        for name in header_list:
            if ' ' in name:
                name = name.split(' ', 1)[1]
                name = re.sub(' ', '_', name)
            
            clean_header.append(name)
        
        self.header = '\t'.join(clean_header)
    
        #return self.header

    def load_sample(self):

        """
        Instance method to load all rows from a file as coverm_contig instances and store them as a 
        coverm_sample instance
        """
        
        self.rows=[]
        with open(self.filename,"r") as file:

            # define coverm file header
            header_line = next(file).strip()
            self.define_clean_header(header_line)

            for line in file:

                items = line.strip().split('\t')

                contig = items[0]
                contig_length = int(items[1])
                read_count = int(items[2])
                reads_per_base = float(items[3])
                #mean = item[4]
                trimmed_mean = float(items[5])
                tpm = float(items[6])
                rpkm = float(items[7])

                        # line_list = line.split('\t')
                        # args_list = line_list[0:6] + line_list[7:9] + line_list[11:13]
                        # args = ' ,'.join([str(item) for item in args_list])
                        # eggnog_orf = Eggnog_orf(args)

                coverm_contig = Coverm_contig(contig, contig_length, read_count, reads_per_base, trimmed_mean, tpm, rpkm)
                self.rows.append(coverm_contig)
                self.total_rpkm += rpkm
                self.total_tm += trimmed_mean




# coverm_file = "0505-0101_coverage_values"

# sample_05050101 = Coverm_sample(coverm_file)
# sample_05050101.load_sample()

# print(sample_05050101.header)
# for row in sample_05050101.rows:
#      print(row.contig)
#      print(row.contig_length)
#      print(row.trimmed_mean) # necesito saber con que parte de la notacion taxonomica me quedo
#      print(row) # lo hace ya separado por tab por el __str__ method. 
#      break




import argparse

def check_arg (args=None) :

    '''
	Description:
        Function collect arguments from command line using argparse
    Input:
        args # command line arguments
    Constant:
        None
    Variables
        parser
    Return
        parser.parse_args() # Parsed arguments
    '''

    parser = argparse.ArgumentParser(prog = 'emapper_profiler', formatter_class=argparse.RawDescriptionHelpFormatter, description= 'Creates 2 tsv files with the ko abundance and kegg_pathway coverage of each sample')

    parser.add_argument('--inputdir','-i', required=True, help='Insert path with _coverage_values and .emapper.annotations files')
    
    parser.add_argument('--outputdir','-o', required=False, default ='results' ,help='The output directory to store the tsv files.')
	
    parser.add_argument('--coverm_suffix', '-c', required=False, default='_coverage_values', help='The suffix of your CoverM files')
    
    parser.add_argument('--version','-v', action='version', version='%(prog)s 0.0.1')

    parser.add_argument('--sample_file', '-s', required=False, help = 'txt file with a list of the samples that need to be processed. If not provided, the sample list will be generated given the input files')

    parser.add_argument('--kegg_dict', '-k', required=False, help = 'Json file of Kegg pathway dictionary') # maybe this needs to be required

    parser.add_argument('--unit', '-u', type=str, choices=['tpm', 'rpkm', 'tm'], default='rpkm', help = 'Specify the output units for the relative abundance results. Default is rpkm')
    
    # on/off flag
    parser.add_argument('--filter_euk', '-e', action='store_true', required=False, help='Option to remove eukaryotes')
    
    parser.add_argument('--novel_fam', '-f', action='store_true', required=False, help='Option to add novel families') 
    
    parser.add_argument('--nf_dir', required=False, help='Input directory that contains the novel families annotations. By default considers a directory called novel_families inside the input direcotyr') 


    # python __main__2.py -i /Users/lucia/Desktop/TFM/scripts/final/data -k /Users/lucia/Desktop/TFM/scripts/parse_KEGGpathway_db/KEGG_pathway_dict.txt -e -u tpm -o results_tpm

    return parser.parse_args()
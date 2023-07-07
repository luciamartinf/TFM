import argparse

def check_arg():

    """
    Function to collect arguments from command line using argparse
    """

    parser = argparse.ArgumentParser(prog = 'emapper_profiler', formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description= 'Creates 2 tsv files with the ko abundance and kegg_pathway coverage of each sample')

    add = parser.add_argument

    add('--inputdir','-i', required=True,
        help='Insert path with _coverage_values and .emapper.annotations files')

    add('--outputdir','-o', default ='results',
        help='The output directory to store the tsv files.')

    add('--coverm_suffix', '-c', default='_coverage_values',
        help='The suffix of your CoverM files')

    add('--version', action='version', version='%(prog)s 0.0.1')

    add('--sample_file', '-s',
        help = 'txt file with a list of the samples that need to be processed. If not provided, the sample list will be generated given the input files')

    add('--kegg_dict', '-k',
        help = 'Directory containing KEGG_pathway and KEGG_kos json dictionaries')

    add('--unit', '-u', type=str, choices=['tpm', 'rpkm', 'tm'], default='rpkm',
        help = 'Specify the output units for the relative abundance results. Default is rpkm')

    # on/off flag
    add('--filter_euk', '-e', action='store_true',
        help='Remove eukaryotes')

    add('--filter_virus', '-v', action='store_true',
        help='Remove viruses')

    add('--novel_fam', '-f', action='store_true',
        help='Add novel families')

    add('--nf_dir',
        help='Input directory that contains the novel families annotations. By default considers a directory called novel_families inside the input directory')


    return parser.parse_args()

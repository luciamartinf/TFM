def check_all_novelfam(eggnog_sample, novelfam_sample, samplename):

    """
    Check how many queries are both in eggnog and novel families outputs
    """

    repeated_queries = []
    #good_reps = []
    good_rep_num = 0

    #print(samplename)

    #print(eggnog_sample.query_list)

    for novelfam_row in novelfam_sample.rows:

        #print(novelfam_row.query)
        query = novelfam_row.query
    
        if query in eggnog_sample.query_list:
            repeated_queries.append(query)
            #print(eggnog_sample.query_dict[query], query)
            if eggnog_sample.query_dict[query] == '-':
                #print('yes')
                #good_reps.append(query)
                good_rep_num += 1
                #print(good_rep_num)
    

    total = len(novelfam_sample.query_list)
    rep = len(repeated_queries)
    percentage = (rep/total)*100
    percent_2 = (good_rep_num/rep)*100


    print(samplename,rep, percentage, percent_2, good_rep_num)

    return repeated_queries

def check_novelfam(eggnog_sample, nf, nf_sample, nf_abundance_dict):

    """
    Check how many queries are both in eggnog and novel families outputs
    """

    query = nf_sample.query_dict[nf]
    
    if query in eggnog_sample.query_list:
        # ko = eggnog_sample.query_dict[query]['ko']
        # cog = eggnog_sample.query_dict[query]['cog']
        # des = eggnog_sample.query_dict[query]['des']
        try:
            nf_abundance_dict[nf]['ko'] += eggnog_sample.all_dict[query]['ko']
            nf_abundance_dict[nf]['cog'] += str(eggnog_sample.all_dict[query]['cog'])
            nf_abundance_dict[nf]['description'] += eggnog_sample.all_dict[query]['description']
        except: 
            nf_abundance_dict[nf]['ko'] = eggnog_sample.all_dict[query]['ko']
            nf_abundance_dict[nf]['cog'] = str(eggnog_sample.all_dict[query]['cog'])
            nf_abundance_dict[nf]['description'] = eggnog_sample.all_dict[query]['description']

    return nf_abundance_dict



def add_nf_abundance(nf_abundance_sample:dict, nf, abundance, total_abun):

    """
    calculate the nf abundance
    """
    try: 
        nf_abundance_sample[nf] += float(abundance)
            
    except:
        nf_abundance_sample[nf] = float(abundance)
        
    total_abun += float(abundance)
    

    return nf_abundance_sample, total_abun

def nf_abundance(novelfam_sample, coverm_sample, samplename, rel_nf_abun:dict, sample_list, eggnog_sample):

    """
    Calculate every nf total abundance in a sample. 
    This function relies on contigs and orf being sorted.
    Takes trimmed_mean value as abundance
    """
    
    i = 0

    total_nf_abun = 0

    nf_abundance_sample = {}
    repeats = {}


    for novelfam_row in novelfam_sample.rows:

        while coverm_sample.rows[i].contig != novelfam_row.contig:
            i += 1
        
        novelfam_row.trimmed_mean = coverm_sample.rows[i].trimmed_mean
        nf_abundance_sample, total_nf_abun = add_nf_abundance(nf_abundance_sample, novelfam_row.novel_fam, novelfam_row.trimmed_mean, total_nf_abun)

        repeats[novelfam_row.novel_fam] = False
        if novelfam_row.query in eggnog_sample.query_list:
            repeats[novelfam_row.novel_fam] = True
            #repeats[novelfam_row.novel_fam].append(samplename)

        

    for nf in nf_abundance_sample.keys():
        try:
            rel_nf_abun[nf][samplename] = nf_abundance_sample[nf]/total_nf_abun
   
        except:
            rel_nf_abun[nf] = {}
            rel_nf_abun[nf]['cog'] = ''
            
            ##
            # rel_nf_abun[nf]['ko'] = '-'
            # rel_nf_abun[nf]['cog'] = '-'
            # rel_nf_abun[nf]['description'] = '-'

            for sample_id in sample_list:
                rel_nf_abun[nf][sample_id] = 0
            
            rel_nf_abun[nf][samplename] = nf_abundance_sample[nf]/total_nf_abun
            
        if repeats[nf]:
            rel_nf_abun[nf]['cog'] += (str(samplename) + ',')

    return rel_nf_abun
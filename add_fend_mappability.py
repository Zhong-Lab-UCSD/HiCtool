"""
Program to add the mappability score to the fragment end (fend) files (one per each chromosome).
add_fend_gc_content.py must be run before.
"""

from time import gmtime, strftime

### Input variables to be updated before running the script

# 1) List of chromosomes of the species you are using
chr_list = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y']

# 2) Path of the fend files generated with add_fend_gc_content.py
restrictionsites_gc_path='insert_path_here'

# 3) Filename used for the fend files generated with add_fend_gc_content.py (same parameter used there)
restrictionsites_gc_filename='restrictionsites_gc'

# 4) Path of the artificial reads files (one per each chromosome)
artificial_reads_path='insert_path_here'

# 5) Path to save the output fend files with mappability score added
output_path='insert_path_here'

# 6) Filename to be used for fend files with mappability score (the chromosome label is automatically added per each file)
restrictionsites_gc_map_filename='restrictionsites_gc_map'

# 7) Number of threads to be used
threads = len(chr_list) # (under the assumption that you have at least 24 threads!)        

print "Adding mappability score in parallel using " + str(threads) + " threads..."
print "Start: " + strftime("%Y-%m-%d %H:%M:%S", gmtime())

def add_mappability_score(chromosome):
    """
        Add for each fragment the mappability score of 500bp upstream and downstream to the restriction site per
        a single chromosome. 
        Parameters:
            chromosome (int): chromosome number (example for chromosome 1: '1').
        
        Output:
        Fend bed file with additional mappability information added in two columns: mappability_upstream and mappability_downstream.
        """
    
    import pybedtools
    import numpy as np
    import os
    import pandas as pd
    
    content = []
    with open(restrictionsites_gc_path + 'chr' + chromosome + '_' + restrictionsites_gc_filename + '.bed') as f:
        for line in f:
            content.append(line.strip().split())
    content = np.array(content)
    
    content[:,4] = np.arange(content.shape[0])
    pos = np.int_(((np.int_(content[:,1]))+3)/5)*5
    
    # Building a bed file from restrictionsites with start coordinate as 500 bp upstream of each position and end coordinate as position.
    content[:,1] = pos-500
    content[:,2] = pos
    content[np.int_(content[:,1]) <= 0,1] = 0
    
    with open(output_path + 'chr' + chromosome + '_upstream.bed', 'w') as f:
        for x in content.tolist():
            f.write("%s\n" % '\t'.join(x))
    
    # Building a bed file from restrictionsites with start coordinate as position and end coordinate as 500 bp downstream of each position.
    content[:,1] = pos
    content[:,2] = pos+500
    
    with open(output_path + 'chr' + chromosome + '_downstream.bed', 'w') as f:
        for x in content.tolist():
            f.write("%s\n" % '\t'.join(x))
    
    a_up = pybedtools.example_bedtool(output_path + 'chr' + chromosome + '_upstream.bed')
    a_down = pybedtools.example_bedtool(output_path + 'chr' + chromosome + '_downstream.bed')
    
    result = pd.read_table(restrictionsites_gc_path + 'chr' + chromosome + '_' + restrictionsites_gc_filename,header=None)
    result.columns = ['chr','start','stop','name','score','strand','gc_up','gc_down']
    result['mappability_upstream'] = '' # adding a field to store the upstrem mappability score
    result['mappability_downstream'] = '' # adding a field to store the downstream mappability score
    
    b = pybedtools.example_bedtool(artificial_reads_path + 'chr' + chromosome + '.txt')

    aup_and_b = a_up.intersect(b,wa=True,wb=True)
    adown_and_b = a_down.intersect(b,wa=True,wb=True)

    aup_and_b.saveas(output_path + 'upstream_intersect_chr' + chromosome + '.bed')
    adown_and_b.saveas(output_path + 'downstream_intersect_chr' + chromosome + '.bed')

    aup_and_b = pd.read_table(output_path + 'upstream_intersect_chr' + chromosome + '.bed', header=None)
    adown_and_b = pd.read_table(output_path + 'downstream_intersect_chr' + chromosome + '.bed', header=None)

    aup_and_b.columns = ['chr','start','stop','name','index','strand','gc_up','gc_down','chr_map','start_map','stop_map','map_up']
    adown_and_b.columns = ['chr','start','stop','name','index','strand','gc_up','gc_down','chr_map','start_map','stop_map','map_down']

    def compute_mappability_score(values):
        return float(sum(values>30))/float(len(values))

    df_up = aup_and_b.groupby('index')['map_up'].apply(compute_mappability_score).to_frame()
    df_down = adown_and_b.groupby('index')['map_down'].apply(compute_mappability_score).to_frame()
    df_up['ind'] = df_up.index
    df_down['ind'] = df_down.index
    df = pd.merge(df_up,df_down,how='outer')

    r1 = 0
    for r in np.array(df.ind):
        if result.ix[r,'strand'] == '+':
            result.ix[r,'mappability_upstream'] = np.array(df.ix[r1,'map_up']).astype(np.str)
            result.ix[r,'mappability_downstream'] = np.array(df.ix[r1,'map_down']).astype(np.str)
        if result.ix[r,'strand'] == '-':
            result.ix[r,'mappability_upstream'] = np.array(df.ix[r1,'map_down']).astype(np.str)
            result.ix[r,'mappability_downstream'] = np.array(df.ix[r1,'map_up']).astype(np.str)
        r1 += 1

    result.score = 1
    result.to_csv(path_or_buf=output_path + 'chr' + chromosome + '_' + restrictionsites_gc_map_filename + '.bed',sep='\t',header=False,index=False)
    
    os.remove(output_path + 'chr' + chromosome + '_upstream.bed')
    os.remove(output_path + 'chr' + chromosome + '_downstream.bed')
    os.remove(output_path + 'upstream_intersect_chr'+ chromosome +'.bed')
    os.remove(output_path + 'downstream_intersect_chr' + chromosome + '.bed')
    print "Mappability score for chromosome " + chromosome + " complete."

# Multiprocessing code execution
from multiprocessing import Pool

if __name__ == '__main__':
    pool = Pool(processes=threads)
    pool.map(add_mappability_score, chr_list)
    
print "End: " + strftime("%Y-%m-%d %H:%M:%S", gmtime())

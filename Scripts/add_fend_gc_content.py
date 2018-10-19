"""
Program to add the gc content to the fragment end (fend) files (one per each chromosome).
"""

from time import gmtime, strftime

### Input variables to be updated before running the script

# 1) List of chromosomes of the species you are using
chr_list = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y']

# 2) Path of the fend files, one per chromosome
restrictionsites_path='insert_path_here'

# 3) Path of the files with the gc content information (one per each chromosome) and named as chr1.txt, chr2.txt, etc.
gc_files_path='insert_path_here'

# 4) Path to save the output fend files with gc content added
output_path='insert_path_here'

# 5) Filename to be used for fend files with gc content (the chromosome label is automatically added per each file)
restrictionsites_gc_filename='restrictionsites_gc'

# 6) Number of threads to be used
threads = len(chr_list) # (under the assumption that you have at least 24 threads!)        

print "Adding GC content information in parallel using " + str(threads) + " threads..."
print "Start: " + strftime("%Y-%m-%d %H:%M:%S", gmtime())

def add_gc_content(chromosome):
    """
        Add for each fragment the GC content of 200bp upstream and downstream to the restriction site per a single
        chromosome. 
        Parameters:
            chromosome (int): chromosome number (example for chromosome 1: '1').

        Output:
        Fend bed file with additional GC content information added in two columns: gc_upstream and gc_downstream.
    """

    import pybedtools
    import numpy as np
    import os
    import pandas as pd
    
    content = []
    with open(restrictionsites_path + 'chr' + chromosome + '.bed') as f:
        for line in f:
            content.append(line.strip().split())
    content = np.array(content)

    content[:,4] = np.arange(content.shape[0])
    pos = np.int_(((np.int_(content[:,1]))+3)/5)*5
    
    # Building a bed file from restrictionsites with start coordinate as 200 bp upstream of each position and end coordinate as position.
    content[:,1] = pos-200
    content[:,2] = pos
    content[np.int_(content[:,1]) <= 0,1] = 0

    with open(output_path + 'chr' + chromosome + '_upstream.bed', 'w') as f:
        for x in content.tolist():
            f.write("%s\n" % '\t'.join(x))

    # Building a bed file from restrictionsites with start coordinate as position and end coordinate as 200 bp downstream of each position.
    content[:,1] = pos
    content[:,2] = pos+200

    with open(output_path + 'chr' + chromosome + '_downstream.bed', 'w') as f:
        for x in content.tolist():
            f.write("%s\n" % '\t'.join(x))    
    
    a_up = pybedtools.example_bedtool(output_path + 'chr' + chromosome + '_upstream.bed')
    a_down = pybedtools.example_bedtool(output_path + 'chr' + chromosome + '_downstream.bed')

    result = pd.read_table(restrictionsites_path + 'chr' + chromosome + '.bed', header=None)
    result.columns = ['chr', 'start','stop','name','score','strand']
    result['gc_upstream'] = '' # adding a field to store the upstream GC content
    result['gc_downstream'] = '' # adding a field to store the upstream GC content
    
    b = pybedtools.example_bedtool(gc_files_path + 'chr' + chromosome + '.txt')

    aup_and_b = a_up.intersect(b,wa=True,wb=True)
    adown_and_b = a_down.intersect(b,wa=True,wb=True)
    
    aup_and_b.saveas(output_path + 'upstream_intersect_chr' + chromosome + '.bed')
    adown_and_b.saveas(output_path + 'downstream_intersect_chr' + chromosome + '.bed')

    aup_and_b = pd.read_table(output_path + 'upstream_intersect_chr' + chromosome + '.bed', header=None)
    adown_and_b = pd.read_table(output_path + 'downstream_intersect_chr' + chromosome + '.bed', header=None)
    
    aup_and_b.columns = ['chr','start','stop','name','index','strand','chr_gc','start_gc','stop_gc','gc_up']
    adown_and_b.columns = ['chr','start','stop','name','index','strand','chr_gc','start_gc','stop_gc','gc_down']
    
    df_up = aup_and_b.groupby('index')['gc_up'].mean().to_frame()
    df_down = adown_and_b.groupby('index')['gc_down'].mean().to_frame()
    df_up.gc_up = df_up.gc_up/100
    df_down.gc_down = df_down.gc_down/100
    df_up['ind'] = df_up.index
    df_down['ind'] = df_down.index
    df = pd.merge(df_up,df_down,how='outer')
    
    r1 = 0
    for r in np.array(df.ind):
        if result.ix[r,'strand'] == '+':
            result.ix[r,'gc_upstream'] = np.array(df.ix[r1,'gc_up']).astype(np.str)
            result.ix[r,'gc_downstream'] = np.array(df.ix[r1,'gc_down']).astype(np.str)
        if result.ix[r,'strand'] == '-':
            result.ix[r,'gc_upstream'] = np.array(df.ix[r1,'gc_down']).astype(np.str)
            result.ix[r,'gc_downstream'] = np.array(df.ix[r1,'gc_up']).astype(np.str)
        r1 += 1

    result.score = 1
    result.to_csv(path_or_buf=output_path + 'chr' + chromosome + '_' + restrictionsites_gc_filename + '.bed',sep='\t',header=False,index=False)

    os.remove(output_path + 'chr' + chromosome + '_upstream.bed')
    os.remove(output_path + 'chr' + chromosome + '_downstream.bed')
    os.remove(output_path + 'upstream_intersect_chr' + chromosome + '.bed')
    os.remove(output_path + 'downstream_intersect_chr' + chromosome + '.bed')
    print "GC content information for chromosome " + chromosome + " complete."
    
    
# Multiprocessing code execution
from multiprocessing import Pool

if __name__ == '__main__':
    pool = Pool(processes=threads)             
    pool.map(add_gc_content, chr_list)
    
print "End: " + strftime("%Y-%m-%d %H:%M:%S", gmtime())
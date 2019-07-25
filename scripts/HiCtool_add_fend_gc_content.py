# Add the GC content of 200bp upstream and downstream to the restriction site per each single chromosome.

# Usage: python2.7 HiCtool_add_fend_gc_content.py [-h] [options]
# Options:
#  -h, --help                   show this help message and exit
#  -c CHROMSIZES_PATH           Path to the folder chromSizes with trailing slash at the end.
#  -s SPECIES                   Species. It has to be one of those present under the chromSizes path. Example: for human hg38 type here "hg38".
#  -r RESTRICTION_SITES_PATH    Path to the folder with the restriction sites bed files.
#  -g GC_FILES_PATH             Path to the folder with the GC files, one per each chromosome.
#  -p THREADS                   Number of parallel threads to use. It has to be less or equal than the number of chromosomes.

# Output files:
#  Fend bed file with additional GC content information added in two columns: gc_upstream and gc_downstream.

from optparse import OptionParser
from time import gmtime, strftime
import pybedtools
import numpy as np
import os
import pandas as pd
from multiprocessing import Pool

parameters = {'chromSizes_path': None,
              'species': None,
              'restrictionsites_path': None,
              'gc_files_path': None,
              'threads': None}

def add_gc_content(chromosome):
    output_path = parameters['restrictionsites_path']
    gc_files_path = parameters['gc_files_path']
    
    # Building a bed file from restrictionsites with start coordinate as 200 bp upstream of each position and end coordinate as position.
    upstream = pd.read_csv(output_path + 'chr' + chromosome + '.bed', sep = '\t', header = None)
    upstream.iloc[:,4] = range(upstream.shape[0])
    upstream.iloc[:,2] = upstream.iloc[:,1]
    upstream.iloc[:,1] = upstream.iloc[:,1] - 200
    upstream.loc[upstream.loc[:,1] < 0, 1] = 0
    a_up = pybedtools.BedTool.from_dataframe(upstream)
    
    # Building a bed file from restrictionsites with start coordinate as position and end coordinate as 200 bp downstream of each position.
    downstream = pd.read_csv(output_path + 'chr' + chromosome + '.bed', sep = '\t', header = None)
    downstream.iloc[:,4] = range(downstream.shape[0])
    downstream.iloc[:,1] = downstream.iloc[:,2]
    downstream.iloc[:,2] = downstream.iloc[:,2] + 200
    a_down = pybedtools.BedTool.from_dataframe(downstream)

    result = pd.read_csv(output_path + 'chr' + chromosome + '.bed', header=None, sep='\t')
    result.columns = ['chr', 'start','stop','name','score','strand']
    result['gc_upstream'] = '' # adding a field to store the upstream GC content
    result['gc_downstream'] = '' # adding a field to store the upstream GC content
    
    b = pybedtools.example_bedtool(gc_files_path + 'chr' + chromosome + '.txt')

    aup_and_b = a_up.intersect(b,wa=True,wb=True).to_dataframe()
    adown_and_b = a_down.intersect(b,wa=True,wb=True).to_dataframe()
    
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
        if result.loc[r,'strand'] == '+':
            result.loc[r,'gc_upstream'] = np.array(df.loc[r1,'gc_up']).astype(np.str)
            result.loc[r,'gc_downstream'] = np.array(df.loc[r1,'gc_down']).astype(np.str)
        if result.loc[r,'strand'] == '-':
            result.loc[r,'gc_upstream'] = np.array(df.loc[r1,'gc_down']).astype(np.str)
            result.loc[r,'gc_downstream'] = np.array(df.loc[r1,'gc_up']).astype(np.str)
        r1 += 1

    result.score = 1
    result.to_csv(path_or_buf=output_path + 'chr' + chromosome + '_gc.bed',sep='\t',header=False,index=False)

    print "GC content information for chr" + chromosome + " complete."
    

if __name__ == '__main__':
    
    usage = 'Usage: python2.7 HiCtool_add_fend_gc_content.py [-h] [options]'
    parser = OptionParser(usage = 'python2.7 %prog -c chromSizes_path -s species -r restriction_sites_path -g gc_files_path -p threads')
    parser.add_option('-c', dest='chromSizes_path', type='string', help='Path to the folder chromSizes with trailing slash at the end.')
    parser.add_option('-s', dest='species', type='string', help='Species. It has to be one of those present under the chromSizes path. Example: for human hg38 type here "hg38".')
    parser.add_option('-r', dest='restrictionsites_path', type='string', help='Path to the folder with the restriction sites bed files.')
    parser.add_option('-g', dest='gc_files_path', type='string', help='Path to the folder with the GC files, one per each chromosome.')
    parser.add_option('-p', dest='threads', type='string', help='Number of parallel threads to use. It has to be less or equal than the number of chromosomes.')
    (options, args) = parser.parse_args( )

    if options.chromSizes_path == None:
        parser.error('-h for help or provide the chromSizes path!')
    else:
        pass
    if options.species == None:
        parser.error('-h for help or provide the species!')
    else:
        pass
    if options.restrictionsites_path == None:
        parser.error('-h for help or provide the restrictionsites bed files path!')
    else:
        pass
    if options.gc_files_path == None:
        parser.error('-h for help or provide the path to the files with GC content information!')
    else:
        pass
    if options.threads == None:
        parser.error('-h for help or provide the number of threads!')
    else:
        pass
    
    parameters['chromSizes_path'] = options.chromSizes_path
    parameters['species'] = options.species
    parameters['restrictionsites_path'] = options.restrictionsites_path
    parameters['gc_files_path'] = options.gc_files_path
    parameters['threads'] = options.threads
    
    if parameters['species'] + ".chrom.sizes" not in os.listdir(parameters['chromSizes_path']):
        available_species = ', '.join([x.split('.')[0] for x in  os.listdir(parameters['chromSizes_path'])])
        parser.error('Wrong species inserted! Check the species spelling or insert an available species: ' + available_species + '. If your species is not listed, please contact Riccardo Calandrelli at <rcalandrelli@eng.ucsd.edu>.')
    
    chromosomes = open(parameters['chromSizes_path'] + parameters['species'] + '.chrom.sizes', 'r')
    chromosomes_list = []
    while True:
        try:
            line2list = next(chromosomes).split('\n')[0].split('\t')
            chromosomes_list.append(line2list[0])
        except StopIteration:
            break
    
    threads = int(parameters['threads'])
    
    if threads > len(chromosomes_list):
        parser.error("Input a number of threads less or equal than the number of chromosomes (" + str(len(chromosomes_list)) + ").")
    else:
        pass
    
    print "Adding GC content information in parallel using " + parameters['threads'] + " threads..."
    print "Start: " + strftime("%Y-%m-%d %H:%M:%S", gmtime())
    pool = Pool(processes=threads)             
    pool.map(add_gc_content, chromosomes_list)
    print "End: " + strftime("%Y-%m-%d %H:%M:%S", gmtime())
    
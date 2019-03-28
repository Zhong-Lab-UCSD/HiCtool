"""
Program to normalize fend data using a parallelized approach.
"""

from time import gmtime, strftime

### Input variables to be updated before running the script

# 1) Object containing learned correction parameters in .hdf5 format obtained with HiCtool_hifive.py
input_file = 'HiC_norm_binning.hdf5'

# 2) List of chromosomes that you are going to parallel process
chr_list = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y']

# 3) Bin size in bp of the contact matrix generated
bin_size = 40000

# 4) 'hg38' or 'mm10' or any other species label in string format
species = 'hg38'

# 5) Dictionary with chromosome sizes of your custom species if you did not use 'hg38' or 'mm10'. Keys: string. Values: int.
# Example: {'1':int1, '2':int2, ... , 'X':intX, 'Y':intY}
chr_sizes = {}

# 6) If true, save the observed contact data
save_obs = True

# 7) If true, save the expected contact data (high dimensional contact matrix)
save_expect = False

# 8) Number of threads to be used
threads = len(chr_list) # (under the assumption that you have at least 24 threads available in this case!)         

print "Normalizing fend contact data in parallel using " + str(threads) + " threads..."
print "Start: " + strftime("%Y-%m-%d %H:%M:%S", gmtime())

def save_matrix(a_matrix, output_file):
    """
    Save an intra-chromosomal contact matrix in the HiCtool compressed format to txt file. 
    1) The upper-triangular part of the matrix is selected (including the
    diagonal).
    2) Data are reshaped to form a vector.
    3) All the consecutive zeros are replaced with a "0" followed by the
    number of times zeros are repeated consecutively.
    4) Data are saved to a txt file.
    Arguments:
        a_matrix (numpy matrix): input contact matrix to be saved
        output_file (str): output file name in txt format
    Output:
        txt file containing the formatted data
    """
    import numpy as np
    n = len(a_matrix)
    iu = np.triu_indices(n)
    vect = a_matrix[iu].tolist()
    with open (output_file,'w') as fout:
        k = len(vect)
        i = 0
        count = 0
        flag = False # flag to set if the end of the vector has been reached
        while i < k and flag == False:
            if vect[i] == 0:
                count+=1
                if (i+count == k):
                    w_out = str(0) + str(count)
                    fout.write('%s\n' %w_out)
                    flag = True
                    break
                while vect[i+count] == 0 and flag == False:
                    count+=1
                    if (i+count == k):
                        w_out = str(0) + str(count)
                        fout.write('%s\n' %w_out)
                        flag = True
                        break
                if flag == False:
                    w_out = str(0) + str(count)
                    fout.write('%s\n' %w_out)
                    i+=count
                    count = 0
            else:
                fout.write('%s\n' %vect[i])
                i+=1


def normalize_chromosome_fend_data_parallel(a_chr):
    """
    --- UPDATED FUNCTION FOR MULTICORE USAGE---
    Normalize the contact data by calculating the corrected reads count for each 
    bin. Observed data, expected fend data and normalized fend data are stored 
    into txt file.
    Parameters:
        a_chr (str): chromosome number (example for chromosome 1: '1').
    """    
    
    import hifive
    import numpy as np
    
    chromosomes = {'hg38':{'1':248956422,
                   '2':242193529,
                   '3':198295559,
                   '4':190214555,
                   '5':181538259,
                   '6':170805979,
                   '7':159345973,
                   '8':145138636,
                   '9':138394717,
                   '10':133797422,
                   '11':135086622,
                   '12':133275309,
                   '13':114364328,
                   '14':107043718,
                   '15':101991189,
                   '16':90338345,
                   '17':83257441,
                   '18':80373285,
                   '19':58617616,
                   '20':64444167,
                   '21':46709983,
                   '22':50818468,
                   'X':156040895,
                   'Y':57227415},
                   'mm10':{'1':195471971,
                   '2':182113224,
                   '3':160039680,
                   '4':156508116,
                   '5':151834684,
                   '6':149736546,
                   '7':145441459,
                   '8':129401213,
                   '9':124595110,
                   '10':130694993,
                   '11':122082543,
                   '12':120129022,
                   '13':120421639,
                   '14':124902244,
                   '15':104043685,
                   '16':98207768,
                   '17':94987271,
                   '18':90702639,
                   '19':61431566,
                   'X':171031299,
                   'Y':91744698}}
    
    chromosome = 'chr' + a_chr
    
    if bin_size >= 1000000:
        bin_size_str = str(bin_size/1000000)
        output_filename = 'HiCtool_' + chromosome + '_' + bin_size_str + 'mb_'
    elif bin_size < 1000000:
        bin_size_str = str(bin_size/1000)
        output_filename = 'HiCtool_' + chromosome + '_' + bin_size_str + 'kb_'
    
    start_pos = 0
    if species == 'hg38' or species == 'mm10':
        end_pos = (chromosomes[species][a_chr]/bin_size)*bin_size
    else:
        end_pos = (chr_sizes[a_chr]/bin_size)*bin_size
        
    # Expected raw (number of possible fend interactions). 
    # These are needed to scale the fend expected data by the mean fend pairs 
    #in each bin.
    
    hic = hifive.HiC(input_file)
    heatmap_raw = hic.cis_heatmap(chrom=chromosome,
                                  start=start_pos,
                                  stop=end_pos,
                                  binsize=bin_size,
                                  arraytype='full',
                                  datatype='raw')
    expected_raw = heatmap_raw[:,:,1]
    n = len(expected_raw)
    scaling_factor = (np.sum(expected_raw)/2)/(n*(n-1)/2) # mean fend pairs in each bin
    
    # Fend data
    hic = hifive.HiC(input_file)
    heatmap_fend = hic.cis_heatmap(chrom=chromosome,
                                   start=start_pos,
                                   stop=end_pos,
                                   binsize=bin_size,
                                   arraytype='full',
                                   datatype='fend')
    
    observed = heatmap_fend[:,:,0] # observed contact data
    if save_obs == True:
        save_matrix(observed, output_filename + 'observed.txt')
    
    # Expected fend (fend corrections)
    expected_fend = heatmap_fend[:,:,1]/scaling_factor # fend correction values
    if save_expect == True:
        save_matrix(expected_fend, output_filename + 'expected_fend.txt')
            
    # In the above calls, all valid possible interactions are queried from 
    # chromosome 'chrom' between 'start' and 'stop' parameters. The 'arraytype' 
    # parameter determines what shape of array data are returned in: 'full' 
    # returns a square, symmetric array of size NxNx2. The 'datatype' parameter
    # specifies which kind of data to extract. The **observed counts** are in 
    # the first index of the last dimension of the returned array (the same 
    # for every 'datatype'), while the **expected counts** are in the second 
    # index of the last dimension.    
    
    # Normalized fend contact matrix
    n = len(expected_fend)
    normalized_fend = np.zeros((n,n))
    for i in xrange(n):
        for j in xrange(n):
            if expected_fend[i][j] == 0:
                normalized_fend[i][j] = 0
            else:
                normalized_fend[i][j] = float(observed[i][j])/float(expected_fend[i][j])
    
    save_matrix(normalized_fend, output_filename + 'normalized_fend.txt')


# Multiprocessing code execution
from multiprocessing import Pool

if __name__ == '__main__':
    pool = Pool(processes=threads)             
    pool.map(normalize_chromosome_fend_data_parallel, chr_list)
    
print "End: " + strftime("%Y-%m-%d %H:%M:%S", gmtime())

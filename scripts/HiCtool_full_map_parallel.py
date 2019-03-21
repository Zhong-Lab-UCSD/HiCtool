"""
Program to generate the global observed contact matrix containing all the matrices (intra and inter) for all the chromosomes using parallel jobs to be normalized using run_ic_mes.sh.

To use this code, an HiC_project_object.hdf5 must be provided (see HiCtool_hifive.py)
"""

### Input variables to be updated before running the script

# 1) Project object file in .hdf5 format obtained with HiCtool_hifive.py
input_file='HiC_project_object.hdf5'

# 2) Bin size in bp of the contact matrix generated
bin_size = 1000000

# 3) 'hg38' or 'mm10' or any other species label in string format
species = 'hg38'

# 4) Dictionary with chromosome sizes of your custom species if you did not use 'hg38' or 'mm10'. Keys: string. Values: int.
# Example: {'1':int1, '2':int2, ... , 'X':intX, 'Y':intY}
custom_species_sizes = {}

# 5) List with the sexual chromosomes of your custom species (if present) as strings.
# Example ['X', 'Y']
sexual_chromosomes = []

# 6) If true, save each single contact matrix
save_each_matrix = False

# 7) Number of threads to be used
threads = 24 # (under the assumption that you have at least 24 threads available in this case!)    

from time import gmtime, strftime
print "Generating the global observed contact matrix in parallel using " + str(threads) + " threads..."
print "Start: " + strftime("%Y-%m-%d %H:%M:%S", gmtime())

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
                   

if species in chromosomes.keys():
        chromosomes_list = [str(i) for i in range(len(chromosomes[species]) - 1)[1:]] + ['X', 'Y']
        chr_dim = []
        for i in chromosomes_list:
            chr_dim.append(chromosomes[species][i]/bin_size) 
        d_chr_dim = {}
        for i in chromosomes_list:
            d_chr_dim[i] = chromosomes[species][i]/bin_size
else:
    if len(sexual_chromosomes) > 0:
        chromosomes_list = [str(i) for i in range(len(custom_species_sizes) - len(sexual_chromosomes) + 1)[1:]]
        chromosomes_list += sexual_chromosomes
    else:
        chromosomes_list = [str(i) for i in range(len(custom_species_sizes) + 1)[1:]]
    chr_dim = []
    for i in chromosomes_list:
        chr_dim.append(custom_species_sizes[i]/bin_size)                
    d_chr_dim = {}
    for i in chromosomes_list:
        d_chr_dim[i] = custom_species_sizes[i]/bin_size


def generate_intrachromosomal_observed_data(a_chr,
                                            bin_size,
                                            input_file='HiC_project_object.hdf5',
                                            species='hg38',
                                            custom_species_sizes={},
                                            save_file=False):
    """
    Generate an observed intrachromosomal contact matrix from HiC_project_object.hdf5
    Parameters:
        a_chr (str): chromosome number (example for chromosome 1: '1').
        bin_size (int): bin size in bp of the contact matrix.
        input_file (str): object containing learned correction parameters in .hdf5 format obtained with
        HiCtool_hifive.py (default: 'HiC_project_object.hdf5').
        species (str): 'hg38' or 'mm10' or any other species label in string format.
        custom_species_sizes (dict): dictionary containing the sizes of the chromosomes
        of your custom species. The keys of the dictionary are chromosomes in string
        format (example for chromosome 1: '1'), the values are chromosome lengths as int.
        save_file (bool): if true, save the observed contact data.
    """
    import hifive
    
    chromosome = 'chr' + a_chr
    
    if bin_size >= 1000000:
        bin_size_str = str(bin_size/1000000)
        output_filename = 'HiCtool_' + chromosome + '_' + bin_size_str + 'mb_' + 'observed_fend'
    elif bin_size < 1000000:
        bin_size_str = str(bin_size/1000)
        output_filename = 'HiCtool_' + chromosome + '_' + bin_size_str + 'kb_' + 'observed_fend'    
    
    if species in chromosomes.keys():
        end_pos = (chromosomes[species][a_chr]/bin_size)*bin_size
    else:
        end_pos = (custom_species_sizes[a_chr]/bin_size)*bin_size
            
    hic = hifive.HiC(input_file)
    heatmap_raw = hic.cis_heatmap(chrom=chromosome,
                                  start=0,
                                  stop=end_pos,
                                  binsize=bin_size,
                                  arraytype='full',
                                  datatype='raw')
    
    observed = heatmap_raw[:,:,0]
    
    if save_file == True:
        save_matrix(observed, output_filename + '.txt')  
    return observed


def generate_interchromosomal_observed_data(chr_row,
                                            chr_col,
                                            bin_size,
                                            input_file='HiC_project_object.hdf5',
                                            species='hg38',
                                            custom_species_sizes={},
                                            save_file=False):
    """
    Generate an observed interchromosomal contact matrix from HiC_project_object.hdf5
    Parameters:
        chr_row (str): chromosome number for the rows (example for chromosome 1: '1').
        chr_col (str): chromosome number for the columns (example for chromosome 1: '1').
        bin_size (int): bin size in bp of the contact matrix.
        input_file (str): object containing learned correction parameters in .hdf5 format obtained with
        HiCtool_hifive.py (default: 'HiC_project_object.hdf5').
        species (str): 'hg38' or 'mm10' or any other species label in string format.
        custom_species_sizes (dict): dictionary containing the sizes of the chromosomes
        of your custom species. The keys of the dictionary are chromosomes in string
        format (example for chromosome 1: '1'), the values are chromosome lengths as int.
        save_file (bool): if true, save the observed contact data.
    """
    import hifive
    
    chromosome_row = 'chr' + chr_row
    chromosome_col = 'chr' + chr_col
    
    if bin_size >= 1000000:
        bin_size_str = str(bin_size/1000000)
        output_filename = 'HiCtool_' + chromosome_row + '_' + chromosome_col + '_' + bin_size_str + 'mb_'
    elif bin_size < 1000000:
        bin_size_str = str(bin_size/1000)
        output_filename = 'HiCtool_' + chromosome_row + '_' + chromosome_col + '_' + bin_size_str + 'kb_'    
    
    if species in chromosomes.keys():
        end_pos_row = (chromosomes[species][chr_row]/bin_size)*bin_size
        end_pos_col = (chromosomes[species][chr_col]/bin_size)*bin_size
    else:
        end_pos_row = (custom_species_sizes[chr_row]/bin_size)*bin_size
        end_pos_col = (custom_species_sizes[chr_col]/bin_size)*bin_size
            
    hic = hifive.HiC(input_file)
    heatmap_raw = hic.trans_heatmap(chromosome_row, chromosome_col, 
                                    start1=0, stop1=end_pos_row, 
                                    start2=0, stop2=end_pos_col,
                                    binsize=bin_size, 
                                    datatype='raw')
    
    observed = heatmap_raw[:,:,0]
    row = observed.shape[0]
    col = observed.shape[1]
    
    if save_file == True:
        row_str = str(row)
        col_str = str(col)
        output_filename = output_filename + row_str + 'x' + col_str + '_'
        save_matrix_rectangular(observed, output_filename + 'observed.txt')
    return observed


def save_matrix_rectangular(a_matrix, output_file):
    """
    Format and save an inter-chromosomal contact matrix in a txt file. 
    1) Data are reshaped to form a vector.
    2) All the consecutive zeros are replaced with a "0" followed by the
    number of times zeros are repeated consecutively.
    3) Data are saved to a txt file.
    Parameters:
    a_matrix (numpy matrix): input contact matrix to be saved
    output_file: output file name in txt format
    Output:
    txt file containing the formatted data
    """
    import numpy as np
    n_row = np.shape(a_matrix)[0]
    n_col = np.shape(a_matrix)[1]
    vect = np.reshape(a_matrix,[1,n_row*n_col]).tolist()[0]
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


def load_matrix_rectangular(input_file,n_row,n_col):
    """
    Load a formatted contact matrix from a txt file and parse it.
    Parameters:
        input_file: input file name in txt format (generated by the function 
        "save_matrix_rectangular")
        n_row (int): number of rows of the matrix
        n_col (int): number of columns of the matrix
    Returns: 
        output_matrix: array containing the parsed values stored 
        in the input txt file to build a contact matrix        
    """
    import numpy as np    
    
    print "Loading " + input_file + "..."
    with open (input_file,'r') as infile:
        matrix_vect = []        
        for i in infile:
            if i[0] == "0" and i[1] != ".":
                for k in xrange(int(i[1:-1])):
                    matrix_vect.append(0)
            else:
                j = i[:-1]            
                matrix_vect.append(float(j))
  
    output_matrix = np.reshape(np.array(matrix_vect),[n_row,n_col])
    return output_matrix
    print "Done!"
     

def save_matrix(a_matrix, output_file):
    """
    Format and save an intra-chromosomal contact matrix in a txt file. 
    1) The upper-triangular part of the matrix is selected (including the
    diagonal).
    2) Data are reshaped to form a vector.
    3) All the consecutive zeros are replaced with a "0" followed by the
    number of times zeros are repeated consecutively.
    4) Data are saved to a txt file.
    Parameters:
    a_matrix (numpy matrix): input contact matrix to be saved
    output_file: output file name in txt format
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


def save_matrix_tab(input_matrix, output_filename):
    """
    Save a contact matrix in a txt file in a tab separated format. Columns are
    separated by tabs, rows are in different lines.
    Parameters:
    input_matrix (numpy matrix): input contact matrix to be saved
    output_filename: output file name in txt format
    Output:
    txt file containing the tab separated data
    """
    with open (output_filename, 'w') as f:
            for i in xrange(len(input_matrix)):
                row = [str(j) for j in input_matrix[i]]
                if i != len(input_matrix) - 1:
                    f.write('\t'.join(row) + '\n')
                else:
                    f.write('\t'.join(row))


def compute_matrix_data_full_observed_parallel(chr_row):
    """
    Function to compute a row of the global matrix, meaning all the contact matrices
    associated to a single chromosome.
    Parameters:
        chr_row (str): chromosome of the row of matrices.
    Returns:
        Txt file of the row of matrices formatted using "save_matrix_rectangular".
    """
    import numpy as np

    chromosome_row = 'chr' + chr_row
    if chromosome_row == 'chr1':
        intra = generate_intrachromosomal_observed_data(chr_row,bin_size,input_file,species,custom_species_sizes,save_each_matrix)
        matrix_full_line = intra
    
    for chr_col in chromosomes_list:
        chromosome_col = 'chr' + chr_col
        
        if chromosome_row == chromosome_col and chromosome_row == 'chr1':
            continue
    
        elif chromosome_row == chromosome_col and chromosome_row != 'chr1':
            intra = generate_intrachromosomal_observed_data(chr_row,bin_size,input_file,species,custom_species_sizes,save_each_matrix)
            n_row = np.shape(intra)[0]
            matrix_full_line = np.concatenate((matrix_full_line,intra),axis=1)
                
        else:
            row = (chromosomes[species][chr_row]/bin_size)*bin_size/bin_size
            col = (chromosomes[species][chr_col]/bin_size)*bin_size/bin_size
            row_str = str(row)
            col_str = str(col)

            matrix_data_full = generate_interchromosomal_observed_data(chr_row,chr_col,bin_size,input_file,species,custom_species_sizes,save_each_matrix)
            n_row = np.shape(matrix_data_full)[0]
            
            if 'matrix_full_line' in locals():
                matrix_full_line = np.concatenate((matrix_full_line,matrix_data_full),axis=1)
            else:
                matrix_full_line = matrix_data_full
    
    save_matrix_rectangular(matrix_full_line, 'matrix_full_line_' + chr_row + '.txt')
    
    
# Multiprocessing code execution
from multiprocessing import Pool
import numpy as np
import os

if __name__ == '__main__':
    pool = Pool(processes=threads)             
    pool.map(compute_matrix_data_full_observed_parallel, chromosomes_list)

matrix_global_observed = load_matrix_rectangular('matrix_full_line_1.txt',d_chr_dim['1'],np.sum(chr_dim))
os.remove('matrix_full_line_1.txt')
for i in chromosomes_list[1:]:
    temp = load_matrix_rectangular('matrix_full_line_' + i + '.txt', d_chr_dim[i], np.sum(chr_dim))
    matrix_global_observed = np.concatenate((matrix_global_observed,temp))
    os.remove('matrix_full_line_' + i + '.txt')

if bin_size >= 1000000:
    bin_size_str = str(bin_size/1000000)
    my_filename = 'HiCtool_' + bin_size_str + 'mb_'
elif bin_size < 1000000:
    bin_size_str = str(bin_size/1000)
    my_filename = 'HiCtool_' + bin_size_str + 'kb_'

save_matrix(matrix_global_observed, my_filename + 'matrix_global_observed.txt')
save_matrix_tab(matrix_global_observed, my_filename + 'matrix_global_observed_tab.txt')
            
with open ('info.txt', 'w') as f:
    f.write('Rows: ' + str(len(matrix_global_observed)) + '\n')
    f.write('Rowsum (average matrix * rows): ' + str(int(np.mean(matrix_global_observed) * len(matrix_global_observed))))
    
print "End: " + strftime("%Y-%m-%d %H:%M:%S", gmtime())

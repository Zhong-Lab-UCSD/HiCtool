# Generate the global observed contact matrix containing all the matrices (intra and inter) for all the chromosomes using parallel jobs to be normalized using run_ic_mes.sh.

# Usage: python2.7 HiCtool_global_map_observed.py [-h] [options]
# Options:
#  -h, --help                           show this help message and exit
#  -i INPUT_FILE                        Project object file in .hdf5 format obtained with HiCtool_hifive.py
#  -o OUTPUT_PATH                       Output path to save the observed contact matrix with trailing slash at the end
#  -b BIN_SIZE                          The bin size (resolution) for the analysis
#  -s SPECIES                           Species. It has to be one of those present under the chromSizes path. Example: for human hg38 type here "hg38"
#  -c CHROMSIZES_PATH                   Path to the folder chromSizes with trailing slash at the end
#  --save_each SAVE_SINGLE_MATRIX       Set to 1 to save each single contact matrix, 0 otherwise
#  -p THREADS                           Number of parallel threads to use. It has to be less or equal than the number of chromosomes

# Output files:
#  Global all-by-all chromosomes observed contact matrix in compressed format and tab separated format.
#  Single contact matrices if SAVE_SINGLE_MATRIX is 1.

from optparse import OptionParser
from time import gmtime, strftime
from multiprocessing import Pool
import numpy as np
import os


parameters = {'input_file': None,
              'output_path': None,
              'bin_size': None,
              'species': None,
              'chromSizes_path': None,
              'save_single_matrix': None,
              'threads': None}  

                
def generate_intrachromosomal_observed_data(a_chr,
                                            bin_size,
                                            input_file,
                                            species='hg38',
                                            save_file=False):
    """
    Generate an observed intrachromosomal contact matrix from HiC_project_object.hdf5.
    Arguments:
        a_chr (str): chromosome number (example for chromosome 1: '1').
        bin_size (int): bin size in bp of the contact matrix.
        input_file (str): object containing learned correction parameters in .hdf5 format obtained with
        HiCtool_hifive.py (default: 'HiC_project_object.hdf5').
        species (str): 'hg38' or 'mm10' or any other species label in string format.
        save_file (bool): if true, save the observed contact data.
    Return: 
        observed intrachromosomal contact matrix in numpy array format.
    Output: 
        observed intrachromosomal contact matrix in HiCtool compressed format if "save_file=True".
    """
    import hifive
    
    chromosome = 'chr' + a_chr
    
    if bin_size >= 1000000:
        bin_size_str = str(bin_size/1000000)
        output_filename = 'HiCtool_' + chromosome + '_' + bin_size_str + 'mb_'
    elif bin_size < 1000000:
        bin_size_str = str(bin_size/1000)
        output_filename = 'HiCtool_' + chromosome + '_' + bin_size_str + 'kb_'
    
    chromosomes = open(parameters['chromSizes_path'] + species + '.chrom.sizes', 'r')
    d_chr_dim = {}
    while True:
        try:
            line2list = next(chromosomes).split('\n')[0].split('\t')
            d_chr_dim[line2list[0]] = int(line2list[1])/bin_size
        except StopIteration:
            break
    
    end_pos = d_chr_dim[a_chr]*bin_size
     
    hic = hifive.HiC(input_file)
    heatmap_raw = hic.cis_heatmap(chrom=chromosome,
                                  start=0,
                                  stop=end_pos,
                                  binsize=bin_size,
                                  arraytype='full',
                                  datatype='raw')
    
    observed = heatmap_raw[:,:,0]
    
    if save_file == True:
        save_matrix(observed, output_filename + 'observed.txt')  
    return observed


def generate_interchromosomal_observed_data(chr_row,
                                            chr_col,
                                            bin_size,
                                            input_file,
                                            species='hg38',
                                            save_file=False):
    """
    Generate an observed interchromosomal contact matrix from HiC_project_object.hdf5
    Arguments:
        chr_row (str): chromosome number for the rows (example for chromosome 1: '1').
        chr_col (str): chromosome number for the columns (example for chromosome 1: '1').
        bin_size (int): bin size in bp of the contact matrix.
        input_file (str): object containing learned correction parameters in hdf5 format obtained with
        HiCtool_hifive.py (default: 'HiC_project_object.hdf5').
        species (str): 'hg38' or 'mm10' or any other species label in string format.
        save_file (bool): if True, save the observed contact data.
    Return: 
        observed interchromosomal contact matrix in numpy array format.
    Output: 
        observed interchromosomal contact matrix in HiCtool compressed format if "save_file=True".
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
    
    chromosomes = open(parameters['chromSizes_path'] + species + '.chrom.sizes', 'r')
    d_chr_dim = {}
    while True:
        try:
            line2list = next(chromosomes).split('\n')[0].split('\t')
            d_chr_dim[line2list[0]] = int(line2list[1])/bin_size
        except StopIteration:
            break
    
    end_pos_row = d_chr_dim[chr_row]*bin_size
    end_pos_col = d_chr_dim[chr_col]*bin_size
    
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
    Save an inter-chromosomal contact matrix in the HiCtool compressed format to txt file.
    1) Data are reshaped to form a vector.
    2) All the consecutive zeros are replaced with a "0" followed by the
    number of times zeros are repeated consecutively.
    3) Data are saved to a txt file.
    Arguments:
        a_matrix (numpy matrix): input contact matrix to be saved
        output_file (str): output file name in txt format
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
    Load an HiCtool compressed rectangular contact matrix from a txt file and parse it.
    Arguments:
        input_file: input file name in txt format (generated by the function 
        "save_matrix_rectangular")
        n_row (int): number of rows of the matrix.
        n_col (int): number of columns of the matrix.
    Return: output_matrix: numpy array the parsed values stored in the input txt file to build a contact matrix.        
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


def save_matrix_tab(input_matrix, output_filename):
    """
    Save a contact matrix in a txt file in a tab separated format. Columns are
    separated by tabs, rows are in different lines.
    Arguments:
        input_matrix (numpy matrix): input contact matrix to be saved
        output_filename (str): output file name in txt format
    Output:
        txt file containing the tab separated data
    """
    with open (output_filename, 'w') as f:
            for i in xrange(len(input_matrix)):
                row = [str(j) for j in input_matrix[i]]
                f.write('\t'.join(row) + '\n')

def compute_matrix_data_full_observed_parallel(chr_row):
    """
    Function to compute an entire row of the global matrix, meaning all the contact matrices
    associated to a single chromosome.
    Arguments:
        chr_row (str): chromosome of the row of matrices.
    Returns: None.
    Output: Txt file of the row of matrices in the compressed HiCtool format.
    """
    import numpy as np
    
    bin_size = int(parameters['bin_size'])
    input_file = parameters['input_file']
    species = parameters['species']
    save_single_matrix = bool(parameters['save_single_matrix'])
    
    chromosomes = open(parameters['chromSizes_path'] + parameters['species'] + '.chrom.sizes', 'r')
    chromosomes_list = []
    while True:
        try:
            line2list = next(chromosomes).split('\n')[0].split('\t')
            chromosomes_list.append(line2list[0])
        except StopIteration:
            break
    
    chromosome_row = 'chr' + chr_row
    if chromosome_row == 'chr' + chromosomes_list[0]:
        intra = generate_intrachromosomal_observed_data(chr_row,bin_size,input_file,species,save_single_matrix)
        matrix_full_line = intra
    
    for chr_col in chromosomes_list:
        chromosome_col = 'chr' + chr_col
        
        if chromosome_row == chromosome_col and chromosome_row == 'chr' + chromosomes_list[0]:
            continue
    
        elif chromosome_row == chromosome_col and chromosome_row != 'chr' + chromosomes_list[0]:
            intra = generate_intrachromosomal_observed_data(chr_row,bin_size,input_file,species,save_single_matrix)
            #n_row = np.shape(intra)[0]
            matrix_full_line = np.concatenate((matrix_full_line,intra),axis=1)
                
        else:
            matrix_data_full = generate_interchromosomal_observed_data(chr_row,chr_col,bin_size,input_file,species,save_single_matrix)
            #n_row = np.shape(matrix_data_full)[0]
            
            if 'matrix_full_line' in locals():
                matrix_full_line = np.concatenate((matrix_full_line,matrix_data_full),axis=1)
            else:
                matrix_full_line = matrix_data_full
    
    save_matrix_rectangular(matrix_full_line, parameters['output_path'] + 'matrix_full_line_' + chr_row + '.txt')
  

if __name__ == '__main__':
    
    usage = 'Usage: python2.7 HiCtool_global_map_observed.py [-h] [options]'
    parser = OptionParser(usage = 'python2.7 %prog -i input_file -o output_path -b bin_size -s species -c chromSizes_path --save_each save_single_matrix -p threads')
    parser.add_option('-i', dest='input_file', type='string', help='Project object file in .hdf5 format obtained with HiCtool_hifive.py.')
    parser.add_option('-o', dest='output_path', type='string', help='Path to save the output files with the trailing slash in the end.')
    parser.add_option('-b', dest='bin_size', type='string', help='Bin size (resolution) of the contact matrix.')
    parser.add_option('-s', dest='species', type='string', help='Species. It has to be one of those present under the chromSizes path. Example: for human hg38 type here "hg38".')  
    parser.add_option('-c', dest='chromSizes_path', type='string', help='Path to the folder chromSizes with trailing slash at the end.')
    parser.add_option('--save_each', dest='save_single_matrix', type='int', default=0, help='Insert 1 to save each single contact matrix, 0 otherwise (default: 0).')  
    parser.add_option('-p', dest='threads', type='int', help='Number of threads to use. It has to be less or equal than the number of chromosomes.')
    (options, args) = parser.parse_args( )
    
    if options.input_file == None:
        parser.error('-h for help or provide the input project object file!')
    else:
        pass
    if options.output_path == None:
        parser.error('-h for help or provide the output path!')
    else:
        pass
    if options.bin_size == None:
        parser.error('-h for help or provide the bin size of the contact matrix!')
    else:
        pass
    if options.chromSizes_path == None:
        parser.error('-h for help or provide the chromSizes path!')
    else:
        pass
    if options.species == None:
        parser.error('-h for help or provide the species!')
    else:
        pass
    if options.threads == None:
        parser.error('-h for help or provide the number of threads!')
    else:
        pass
    
    parameters['input_file'] = options.input_file
    parameters['output_path'] = options.output_path
    parameters['bin_size'] = options.bin_size
    parameters['chromSizes_path'] = options.chromSizes_path
    parameters['species'] = options.species
    parameters['save_single_matrix'] = options.save_single_matrix
    parameters['threads'] = options.threads
    
    if parameters['species'] + ".chrom.sizes" not in os.listdir(parameters['chromSizes_path']):
        available_species = ', '.join([x.split('.')[0] for x in  os.listdir(parameters['chromSizes_path'])])
        parser.error('Wrong species inserted! Check the species spelling or insert an available species: ' + available_species + '. If your species is not listed, please contact Riccardo Calandrelli at <rcalandrelli@eng.ucsd.edu>.')
    
    bin_size = int(parameters['bin_size'])
    chromosomes = open(parameters['chromSizes_path'] + parameters['species'] + '.chrom.sizes', 'r')
    chromosomes_list = []
    chr_dim = []
    d_chr_dim = {}
    while True:
        try:
            line2list = next(chromosomes).split('\n')[0].split('\t')
            chromosomes_list.append(line2list[0])
            chr_dim.append(int(line2list[1])/bin_size)
            d_chr_dim[line2list[0]] = int(line2list[1])/bin_size
        except StopIteration:
            break
    
    threads = parameters['threads']
    
    if threads > len(chromosomes_list):
        parser.error("Input a number of threads less or equal than the number of chromosomes (" + str(len(chromosomes_list)) + ").")
    else:
        pass
    
    output_path = parameters['output_path']
    
    if threads > 1:
        print "Generating the global observed contact matrix in parallel using " + str(threads) + " threads..."
        print "Start: " + strftime("%Y-%m-%d %H:%M:%S", gmtime())
        pool = Pool(processes=threads)             
        pool.map(compute_matrix_data_full_observed_parallel, chromosomes_list)
    else:
        print "Generating the global observed contact matrix with a single thread..."
        print "Start: " + strftime("%Y-%m-%d %H:%M:%S", gmtime())
        for c in chromosomes_list:
            compute_matrix_data_full_observed_parallel(c)
    
    matrix_global_observed = load_matrix_rectangular(output_path + 'matrix_full_line_' + chromosomes_list[0] + '.txt', d_chr_dim[chromosomes_list[0]], np.sum(chr_dim))
    os.remove(output_path + 'matrix_full_line_' + chromosomes_list[0] + '.txt')
    for i in chromosomes_list[1:]:
        temp = load_matrix_rectangular(output_path + 'matrix_full_line_' + i + '.txt', d_chr_dim[i], np.sum(chr_dim))
        matrix_global_observed = np.concatenate((matrix_global_observed,temp))
        os.remove(output_path + 'matrix_full_line_' + i + '.txt')
    
    print "Building global matrix... "
    if bin_size >= 1000000:
        bin_size_str = str(bin_size/1000000) + 'mb'
        my_filename = 'HiCtool_' + bin_size_str + "_"
    elif bin_size < 1000000:
        bin_size_str = str(bin_size/1000) + 'kb'
        my_filename = 'HiCtool_' + bin_size_str + "_"
    
    save_matrix(matrix_global_observed, output_path + my_filename + 'matrix_global_observed.txt')
    save_matrix_tab(matrix_global_observed, output_path + my_filename + 'matrix_global_observed_tab.txt')
    print "Done!"
    with open (output_path + 'info_matrix_global_observed_' + bin_size_str + '.txt', 'w') as f:
        f.write('Rows: ' + str(len(matrix_global_observed)) + '\n')
        f.write('Rowsum (average matrix * rows): ' + str(int(np.mean(matrix_global_observed) * len(matrix_global_observed))))
    print 'Rows: ' + str(len(matrix_global_observed))
    print 'Rowsum (average matrix * rows): ' + str(int(np.mean(matrix_global_observed) * len(matrix_global_observed)))
    print "End: " + strftime("%Y-%m-%d %H:%M:%S", gmtime())

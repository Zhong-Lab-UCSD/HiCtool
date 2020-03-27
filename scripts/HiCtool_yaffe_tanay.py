# Program to:
# - Normalize the contact data with Yaffe-Tanay approach (fend and enrichment ("observed / expected")).
# - Plot normalized fend contact matrix and histogram of contact distribution.
# - Plot enrichement contact matrix and histogram of contact distribution.
# - Plot correlation matrix containing the Pearson correlation coeffiecients calculated from the O/E matrix.
#
# To use this code, fend correction values must be provided (see HiCtool_hifive.py, -m Yaffe-Tanay)

# Usage: python2.7 HiCtool_yaffe_tanay.py [-h] [options]
# Options:
#  -h, --help                show this help message and exit
#  --action                  Action to perform: normalize_fend, plot_map, normalize_enrich, plot_enrich, plot_correlation.
#  -i INPUT_FILE             Input contact matrix for plotting actions, norm binning hdf5 object from HiFive for normalizing actions.
#  -c CHROMSIZES_PATH        Path to the folder chromSizes with trailing slash at the end.
#  -b BIN_SIZE               Bin size (resolution) of the contact matrix.
#  -s SPECIES                Species. It has to be one of those present under the chromSizes path.
#  --processors              Processors to be used to normalize the data in parallel.  
#  --chr                     Single chromosome for normalization or plotting, or list of chromosomes between square brackets for normalization of multiple chromosomes at once.
#  --data_type               Data type to label your data, example: observed, normalized_fend, normalized_enrich.
#  --coord                   List of two integers with start and end coordinates for the chromosome to be plotted.
#  --save_obs                Insert 1 to save the observed data when normalizing, 0 otherwise. Default: 1.
#  --save_expect             Insert 1 to save the expected data when normalizing, 0 otherwise. Default: 0.
#  --my_colormap             Colormap to be used to plot the contact data. You can choose among any colorbar here https://matplotlib.org/examples/color/colormaps_reference.html, or input a list of colors if you want a custom colorbar. Example: [white, red, black]. Colors can be specified also HEX format. Default: [white,red].  
#  --cutoff_type             When plotting the contact data, to select a type of cutoff (percentile or contact_number) or plot the full range of the data (leave it empty). Default: percentile. 
#  --cutoff                  When plotting the contact data, to set a maximum cutoff on the number of contacts for the colorbar based on cutoff_type. Default: 95.
#  --cutoff_max              When plotting the enrichment data, log2 upper cutoff for the colorbar (every enrichment value above cutoff_max is plotted in red).
#  --cutoff_min              When plotting the enrichment data, log2 lower cutoff (negative value) for the colorbar (every enrichment value below cutoff_min is plotted in blue).  
#  --max_color               When plotting the contact data, to set the color of the bins with contact counts over "cutoff". Default: #460000.
#  --my_dpi                  Resolution of the contact map in dpi. Default: 1000.  
#  --plot_histogram          Insert 1 to plot the histogram of the data distribution, 0 otherwise. Default: 0.
#  --topological_domains     Topological domain coordinates file (as generated from HiCtool_TAD_analysis.py) to visualize domains on the heatmap if action is "plot_map".
#  --domain_color            To set the color for topological domains on the heatmap. Default: #0000ff.


from optparse import OptionParser
import numpy as np
import os
import os.path
from os import path
from time import gmtime, strftime
from multiprocessing import Pool

parameters = {'action': None,
              'input_file': None,
              'chromSizes_path': None,
              'chr': None,
              'species': None,
              'processors': None,
              'bin_size': None,
              'save_obs': None,
              'save_expect': None,
              'data_type': None,
              'coord': None,
              'my_colormap': None,
              'cutoff_type': None,
              'cutoff': None,
              'cutoff_max': None,
              'cutoff_min': None,
              'max_color': None,
              'my_dpi': None,
              'plot_histogram': None,
              'topological_domains': None,
              'domain_color': None              
              }


def save_matrix(a_matrix, output_file):
    """
    Save an intra-chromosomal contact matrix in the HiCtool compressed format to txt file.
    1) The upper-triangular part of the matrix is selected (including the diagonal).
    2) Data are reshaped to form a vector.
    3) All the consecutive zeros are replaced with a "0" followed by the
    number of times zeros are repeated consecutively.
    4) Data are saved to a txt file.
    Arguments:
        a_matrix (numpy matrix): input contact matrix to be saved
        output_file (str): output file name in txt format.
    Output: 
        txt file containing the formatted data.
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


def load_matrix(input_file):
    """
    Load an HiCtool compressed square (and symmetric) contact matrix from a txt file and parse it.
    Arguments:
        input_file (str): input file name in txt format (generated by the function 
        "save_matrix").
    Return: 
        numpy array containing the parsed values stored in the input txt file to build a contact matrix.      
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
  
    k = len(matrix_vect)
    matrix_size = int((-1+np.sqrt(1+8*k))/2)
    
    iu = np.triu_indices(matrix_size)
    output_matrix_1 = np.zeros((matrix_size,matrix_size)) # upper triangular plus the diagonal
    output_matrix_1[iu] = matrix_vect
    
    diag_matrix = np.diag(np.diag(output_matrix_1)) # diagonal
    output_matrix_2 = np.transpose(output_matrix_1) # lower triangular plus the diagonal
    output_matrix = output_matrix_1 + output_matrix_2 - diag_matrix
    print "Done!"
    return output_matrix

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

                    
def load_matrix_tab(input_file):
    """
    Load a contact matrix saved in a tab separated format using the function
    "save_matrix_tab".
    Arguments:
        input_file (str): input contact matrix to be loaded.
    Return: 
        numpy array containing the parsed values stored in the input tab separated txt file to build a contact matrix.
    """
    import numpy as np
    
    print "Loading " + input_file + "..."
    with open (input_file, 'r') as infile:
        lines = infile.readlines()
        temp = []
        for line in lines:
            row = [float(i) for i in line.strip().split('\t')]
            temp.append(row)
            
        output_matrix = np.array(temp)
    print "Done!"
    return output_matrix

def load_topological_domains(input_file):
    """
    Function to load the topological domains coordinates from txt file.
    Arguments:
        input_file (str): input file name generated with "calculate_topological_domains" in txt format.
    Return:
        List of lists with topological domain coordinates.
    """
    import csv
    print "Loading topological domain coordinates..."
    with open(input_file, 'r') as f:
        reader = csv.reader(f, dialect='excel', delimiter='\t')
        topological_domains = []
        for row in reader:
            row_int = [int(x) for x in row]
            topological_domains.append(row_int)
        print "Done!"
        return topological_domains
    

def normalize_chromosome_fend_data(a_chr):
    """
    Normalize the contact data by calculating the corrected reads count for each 
    bin. Observed data and expected fend data (correction data) can be saved to txt file.
    Arguments:
        a_chr (str): chromosome number (example for chromosome 1: '1').
    Return:
        Normalized fend contact matrix.
    Outputs:
        Txt file with the normalized enrichment contact matrix saved in the HiCtool compressed format.
        Txt file with the observed contact matrix saved in the HiCtool compressed format if "save_obs=True".
        Txt file with the expected contact matrix saved in the HiCtool compressed format if "save_expect=True".
    """
    import hifive
    import numpy as np
    
    bin_size = parameters['bin_size']
    input_file = parameters['input_file']
    save_obs = bool(parameters['save_obs'])
    save_expect = bool(parameters['save_expect'])
    
    chromosome = 'chr' + a_chr
    print "Normalizing fend data " + chromosome + " ..."
    
    output_path = "yaffe_tanay_" + str(bin_size)
    if not path.exists(output_path):
        os.mkdir(output_path)
    output_filename = chromosome + '_' + str(bin_size) + '_'

    chromosomes = open(parameters['chromSizes_path'] + parameters['species'] + '.chrom.sizes', 'r')
    d_chr_dim = {}
    while True:
        try:
            line2list = next(chromosomes).split('\n')[0].split('\t')
            d_chr_dim[line2list[0]] = int(line2list[1])/bin_size
        except StopIteration:
            break
    
    start_pos = 0
    end_pos = d_chr_dim[a_chr]*bin_size
        
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
    scaling_factor = float(np.sum(expected_raw)/2.0)/float(n*(n-1)/2) # mean fend pairs in each bin
    
    # Fend data
    hic = hifive.HiC(input_file)
    heatmap_fend = hic.cis_heatmap(chrom=chromosome,
                                   start=start_pos,
                                   stop=end_pos,
                                   binsize=bin_size,
                                   arraytype='full',
                                   datatype='fend')
    
    observed = heatmap_fend[:,:,0] # observed contact data extracted from the heatmap object
    if save_obs == True:
        save_matrix(observed, output_path + "/" + output_filename + 'observed.txt') 

    
    # Expected fend (fend corrections)
    expected_fend = heatmap_fend[:,:,1]/scaling_factor # fend correction values
    if save_expect == True:
        save_matrix(expected_fend, output_path + "/" + output_filename + 'expected_fend.txt') 
 
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
    
    save_matrix(normalized_fend, output_path + "/" + output_filename + 'normalized_fend.txt') 
   
    print "Done!"
    return normalized_fend


def plot_chromosome_data(contact_matrix,
                        a_chr,
                        bin_size,
                        coord=None, 
                        species='hg38',
                        data_type='normalized_fend',
                        my_colormap=['white', 'red'],
                        cutoff_type='percentile',
                        cutoff=95,
                        max_color='#460000',
                        my_dpi=1000,
                        plot_histogram=False,
                        topological_domains=None,
                        domain_color='#0000ff'):
    """
    Plot a contact map and histogram of the contact distribution for observed data, normalized fend data, expected fend and enrichment data.
    Arguments:
        contact_matrix (str | obj): txt file of the HiCtool compressed contact matrix or contact matrix object.
        a_chr (str): chromosome number (example for chromosome 1: '1').
        bin_size (int): bin size in bp of the contact matrix.
        coord (list): list with two integers start and end coordinate for the plot in bp.
        species (str): species label in string format.
        data_type (str): type of data to plot either "observed", "normalized_fend", "expected_fend", "expected_enrich".
        my_colormap (str | list): colormap to be used to plot the data. 1) Use a string if you choose among any colorbar here 
        https://matplotlib.org/examples/color/colormaps_reference.html 2) Use a list of strings with colors if you want
        a custom colorbar. Example: ['white', 'red', 'black']. Colors can be specified also in this format: '#000000'.
        cutoff_type (str): to select a type of cutoff ('percentile' or 'contact_number') or plot the full range of the data (set the 
        parameter as None).
        cutoff (int): percentile to set a maximum cutoff on the number of contacts for the colorbar.
        max_color (str): to set the color of the bins with contact counts over "cutoff".
        my_dpi (int): resolution of the contact map in dpi.
        plot_histogram (bool): if True, plot and save to file the histogram of the contact distribution.
        topological_domains (str | obj): topological domain coordinates to visualize domains on the heatmap. 
        They can be passed either as a txt file or object (as generated from HiCtool_TAD_analysis.py) If None, no topological domains.
        domain_color (str): to set the color for topological domains on the heatmap.
    Outputs:
        Contact map saved to pdf format.
        Histogram of the data distribution saved to pdf format if "plot_histogram=True".
    """        
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    import copy

    chromosome = 'chr' + a_chr
    
    output_path = "yaffe_tanay_" + str(bin_size)
    if not path.exists(output_path):
        os.mkdir(output_path)
    output_filename = chromosome + '_' + str(bin_size) + '_' + data_type

    chromosomes = open(parameters['chromSizes_path'] + parameters['species'] + '.chrom.sizes', 'r')
    d_chr_dim = {}
    while True:
        try:
            line2list = next(chromosomes).split('\n')[0].split('\t')
            d_chr_dim[line2list[0]] = int(line2list[1])/bin_size
        except StopIteration:
            break
    
    end_pos = d_chr_dim[a_chr]*bin_size
    
    # Plotting of the data
    if isinstance(contact_matrix, str):
        matrix_data_full = load_matrix(contact_matrix)
        print "Plotting " + contact_matrix + "..."
    else:
        print "Plotting contact matrix " + chromosome + " ..."
        matrix_data_full = copy.deepcopy(contact_matrix)
    
    # Update matrix values to plot topological domains
    if topological_domains != None:
        if bin_size > 40000:
            print "ERROR! To plot topological domains on the heatmap the bin size should be 40000 or lower."
            return
        if isinstance(topological_domains, str):
            domains = load_topological_domains(topological_domains)
        else:
            domains = topological_domains
        output_filename = output_filename + '_domains'
        diag_index = np.diag_indices(len(matrix_data_full))
        for domain in domains:
            temp_start = domain[0]/bin_size
            temp_end = domain[1]/bin_size
            matrix_data_full[temp_start,temp_start:temp_end] = -1
            matrix_data_full[temp_start:temp_end,temp_end-1] = -1
            matrix_data_full[(diag_index[0][temp_start:temp_end],diag_index[1][temp_start:temp_end])] = -1
    
    # Selecting a part
    if coord != None:
        if len(coord) != 2:
            print "ERROR! Both start and end coordinate has to be declared! Otherwise leave both undeclared to plot the entire contact matrix."
            return
        else:
            print "Selecting part [" + str(coord[0]) + "-" + str(coord[1]) + "] ..."
            output_filename += "_" + str(coord[0]) + "_" + str(coord[1]) 
            start_bin = coord[0]/bin_size
            end_bin = coord[1]/bin_size
        
            if coord[0] >= coord[1]:
                print "ERROR! Start coordinate should be lower than end coordinate"
                return
            
            if start_bin >= end_bin:
                print "ERROR! Start coordinate should be much lower than the end coordinate given the bin size"
                return
        
            if coord[1] > end_pos:
                print "ERROR! End coordinate is larger than chromosome size " + str(end_pos) + " bp"
                return
    else: # Full contact matrix
        start_bin = 0
        end_bin = d_chr_dim[a_chr]

    matrix_data_full = matrix_data_full[start_bin:end_bin+1,start_bin:end_bin+1] 
    
    n = len(matrix_data_full)
    output_vect = np.reshape(matrix_data_full,n*n,1)
    non_zero = np.nonzero(output_vect)
    if non_zero[0].size == 0:
        print "ERROR! The portion of chromosome you selected contains no data."
        return
    
    # Heatmap plotting
    def format_e(n):
        a = '%e' % n
        return a.split('e')[0].rstrip('0').rstrip('.') + 'e' + a.split('e')[1]
    
    if bin_size >= 1000000:
        bin_size_str = str(bin_size/1000000)
        my_bin_size = bin_size_str + ' mb'
    elif bin_size < 1000000:
        bin_size_str = str(bin_size/1000)
        my_bin_size = bin_size_str + ' kb'   
        
    if isinstance(my_colormap, list):
        from matplotlib.colors import LinearSegmentedColormap
        my_cmap = LinearSegmentedColormap.from_list('mycmap', my_colormap)
    elif isinstance(my_colormap, str):
        my_cmap = my_colormap
    
    if cutoff_type == 'percentile':
        perc = np.percentile(output_vect[non_zero[0]],cutoff)
    elif cutoff_type == 'contact_number':
        perc = cutoff
        if cutoff >= np.max(matrix_data_full):
            print "Cut-off value greater than the maximum number of contacts! Set a lower one."
            return
    
    plt.close("all")
    #plt.gcf().subplots_adjust(left=0.15)
    #plt.gcf().subplots_adjust(bottom=0.15)
    if cutoff_type == None:
        if topological_domains == None:
            plt.imshow(matrix_data_full, cmap=my_cmap, interpolation='nearest')
            cbar = plt.colorbar()
        else:
            plt.imshow(matrix_data_full, cmap=my_cmap, interpolation='nearest', vmin=0)
            cbar = plt.colorbar()
            cbar.cmap.set_under(domain_color)
        
    elif cutoff_type == 'percentile' or cutoff_type == 'contact_number':
        if topological_domains == '':
            plt.imshow(matrix_data_full, cmap=my_cmap, interpolation='nearest', vmax=perc)
            cbar = plt.colorbar(extend='max')
            cbar.cmap.set_over(max_color)
        else:
            plt.imshow(matrix_data_full, cmap=my_cmap, interpolation='nearest', vmax=perc, vmin=0)
            cbar = plt.colorbar(extend='max')
            cbar.cmap.set_over(max_color)
            cbar.cmap.set_under(domain_color)
    
    plt.title(data_type + ' contact map (' + my_bin_size + ')', fontsize=12)
    plt.xlabel(chromosome + ' coordinate (bp)', fontsize=10)
    plt.ylabel(chromosome + ' coordinate (bp)', fontsize=10)
    cbar.ax.set_ylabel(data_type + ' contact counts', rotation=270, labelpad=20)
    if coord != None:
        ticks = (np.arange(0, n, n/4) * bin_size) + coord[0]
    else:
        ticks = (np.arange(0, n, n/4) * bin_size)
    format_ticks = [format_e(i) for i in ticks.tolist()]
    plt.xticks(np.arange(0, n, n/4), format_ticks)
    plt.yticks(np.arange(0, n, n/4), format_ticks)
    plt.xticks(fontsize=8)
    plt.yticks(fontsize=8)
    plt.tick_params(axis='both', which='both', direction='out', top=False, right=False)
    plt.tight_layout()
    plt.savefig(output_path + "/" + output_filename + '.pdf', format = 'pdf', dpi=my_dpi)
        
    # Plot of the histogram
    if plot_histogram:
        histogram = []
        k = 1
        for i in xrange(n):
            row = matrix_data_full[i][k:]
            for j in row:
                histogram.append(j)
            k += 1
        
        plt.close("all")
        histogram_bins = int(pow(len(histogram),0.3))
        plt.hist(histogram, bins=histogram_bins)
        plt.title(data_type + ' contact counts distribution', fontsize=16)
        plt.xlabel(data_type + ' contact counts', fontsize=14)
        plt.ylabel('Number of bins', fontsize=14)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.tight_layout()
        plt.savefig(output_path + "/" + output_filename + '_histogram.pdf', format = 'pdf')
    print "Done!"


def normalize_chromosome_enrich_data(a_chr):
    """
    Calculate the enrichment data as "observed/expected" where the expected reads
    count is for each bin considering the linear distance between read pairs and the learned
    correction parameters. Observed and expected contact data can be saved
    to txt files.
    Arguments:
        a_chr (str): chromosome number (example for chromosome 1: '1').
    Return: 
        Normalized enrichment contact matrix.
    Outputs:
        Txt file with the normalized enrichment contact matrix saved in the HiCtool compressed format.
        Txt file with the Pearson Correlation contact matrix saved in the HiCtool compressed format.
        Txt file with the observed contact matrix saved in the HiCtool compressed format if "save_obs=True".
        Txt file with the expected contact matrix saved in the HiCtool compressed format if "save_expect=True".
        
    """
    import hifive
    import numpy as np
    
    bin_size = parameters['bin_size']
    input_file = parameters['input_file']
    save_obs = bool(parameters['save_obs'])
    save_expect = bool(parameters['save_expect'])
    
    print "Normalizing enrichment data..."
    chromosome = 'chr' + a_chr
    
    output_path = "yaffe_tanay_" + str(bin_size)
    if not path.exists(output_path):
        os.mkdir(output_path)
    output_filename = chromosome + '_' + str(bin_size) + '_'

    start_pos = 0
    chromosomes = open(parameters['chromSizes_path'] + parameters['species'] + '.chrom.sizes', 'r')
    d_chr_dim = {}
    while True:
        try:
            line2list = next(chromosomes).split('\n')[0].split('\t')
            d_chr_dim[line2list[0]] = int(line2list[1])/bin_size
        except StopIteration:
            break
    
    end_pos = d_chr_dim[a_chr]*bin_size

    # Enrichment data
    hic = hifive.HiC(input_file)
    heatmap_enrich = hic.cis_heatmap(chrom=chromosome,
                                     start=start_pos,
                                     stop=end_pos,
                                     binsize=bin_size,
                                     arraytype='full',
                                     datatype='enrichment')
    
    # Observed data
    observed = heatmap_enrich[:,:,0] # observed contact data extracted from the heatmap object
    if save_obs == True:
        save_matrix(observed, output_path + "/" + output_filename + 'observed.txt') 
     
    # Expected enrichment data (fend corrections and distance property)
    expected_enrich = heatmap_enrich[:,:,1] # expected enrichment contact data extracted from the heatmap object
    if save_expect == True:
        save_matrix(expected_enrich, output_path + "/" + output_filename + 'expected_enrich.txt') 

    # Normalized enrichment contact matrix
    n = len(expected_enrich)
    normalized_enrich = np.zeros((n,n))
    for i in xrange(n):
        for j in xrange(n):
            if expected_enrich[i][j] == 0:
                normalized_enrich[i][j] = -1
            else:
                normalized_enrich[i][j] = float(observed[i][j])/float(expected_enrich[i][j])

    save_matrix(normalized_enrich, output_path + "/" + output_filename + 'normalized_enrich.txt') 
    
    # normalized_enrich[normalized_enrich == -1] = 0
    pcc_contact_matrix = np.corrcoef(normalized_enrich)
    pcc_contact_matrix[np.isnan(pcc_contact_matrix)] = 0
    save_matrix(pcc_contact_matrix, output_path + "/" + output_filename + 'correlation_matrix.txt') 
    
    print "Done!"
    return normalized_enrich


def plot_chromosome_enrich_data(contact_matrix,
                                a_chr,
                                bin_size,
                                coord=None, 
                                species='hg38',
                                cutoff_max=None,
                                cutoff_min=None,
                                my_dpi=1000,
                                plot_histogram=False):
    """
    Plot the log2 of the "observed / expected" contact map and histogram of the enrichment values distribution 
    generated with the function "normalize_chromosome_enrich_data".
    Arguments:
        contact_matrix (str | obj): txt file of the "observed / expected" contact matrix generated with the function 
        "normalize_chromosome_enrich_data" or contact matrix returned by the function "normalize_chromosome_enrich_data".
        a_chr (str): chromosome number (example for chromosome 1: '1').
        bin_size (int): bin size in bp of the contact matrix.
        coord (list): list with two integers start and end coordinate for the plot in bp.
        species (str): 'hg38' or 'mm10' or any other species label in string format.
        cutoff_max (float): log2 upper cutoff for the colorbar (every enrichment value above cutoff_max 
        is plotted in red). Set to None to do not have any cutoff.
        cutoff_min (float): log2 lower cutoff (negative value) for the colorbar (every enrichment value 
        below cutoff_min is plotted in blue). Set to None to do not have any cutoff.
        my_dpi (int): resolution of the contact map in dpi.
        plot_histogram (bool): if True, plot the histogram.
    Outputs:
        Enrichment contact map saved to pdf format.
        Histogram of the enrichment data distribution saved to pdf format if "plot_histogram=True".
    """                                       
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    import math
    from numpy import ma
    from matplotlib import cbook
    from matplotlib.colors import Normalize
    import copy

    chromosome = 'chr' + a_chr
    output_path = "yaffe_tanay_" + str(bin_size)
    if not path.exists(output_path):
        os.mkdir(output_path)
    output_filename = chromosome + '_' + str(bin_size) + '_normalized_enrich'
    
    chromosomes = open(parameters['chromSizes_path'] + parameters['species'] + '.chrom.sizes', 'r')
    d_chr_dim = {}
    while True:
        try:
            line2list = next(chromosomes).split('\n')[0].split('\t')
            d_chr_dim[line2list[0]] = int(line2list[1])/bin_size
        except StopIteration:
            break
    
    end_pos = d_chr_dim[a_chr]*bin_size
    
    # Plotting the enrichment contact data
    if isinstance(contact_matrix, str):
        matrix_data_full = load_matrix(contact_matrix)
        print "Plotting " + contact_matrix + "..."
    else:
        print "Plotting contact matrix " + chromosome + " ..."
        matrix_data_full = copy.deepcopy(contact_matrix)

    # Selecting a part
    if coord != None:
        if len(coord) != 2:
            print "ERROR! Both start and end coordinate has to be declared! Otherwise leave both undeclared to plot the entire contact matrix."
            return
        else:
            print "Selecting part [" + str(coord[0]) + "-" + str(coord[1]) + "] ..."
            output_filename += "_" + str(coord[0]) + "_" + str(coord[1]) 
            start_bin = coord[0]/bin_size
            end_bin = coord[1]/bin_size
        
            if coord[0] >= coord[1]:
                print "ERROR! Start coordinate should be lower than end coordinate"
                return
            
            if start_bin >= end_bin:
                print "ERROR! Start coordinate should be much lower than the end coordinate given the bin size"
                return
        
            if coord[1] > end_pos:
                print "ERROR! End coordinate is larger than chromosome size " + str(end_pos) + " bp"
                return
    else: # Full contact matrix
        start_bin = 0
        end_bin = d_chr_dim[a_chr]
    
    matrix_data_full = matrix_data_full[start_bin:end_bin+1,start_bin:end_bin+1]
    n = len(matrix_data_full)
    
    # Heatmap plotting
    # Defining a class to generate a divergent colorbar with a custom midpoint
    class MidPointNorm(Normalize):    
        def __init__(self, midpoint=0, vmin=None, vmax=None, clip=False):
            Normalize.__init__(self,vmin, vmax, clip)
            self.midpoint = midpoint
    
        def __call__(self, value, clip=None):
            if clip is None:
                clip = self.clip
    
            result, is_scalar = self.process_value(value)
    
            self.autoscale_None(result)
            vmin, vmax, midpoint = self.vmin, self.vmax, self.midpoint
    
            if not (vmin < midpoint < vmax):
                raise ValueError("midpoint must be between maxvalue and minvalue.")       
            elif vmin == vmax:
                result.fill(0)
            elif vmin > vmax:
                raise ValueError("maxvalue must be bigger than minvalue")
            else:
                vmin = float(vmin)
                vmax = float(vmax)
                if clip:
                    mask = ma.getmask(result)
                    result = ma.array(np.clip(result.filled(vmax), vmin, vmax),
                                      mask=mask)
    
                resdat = result.data
    
                #First scale to -1 to 1 range, than to from 0 to 1.
                resdat -= midpoint            
                resdat[resdat>0] /= abs(vmax - midpoint)            
                resdat[resdat<0] /= abs(vmin - midpoint)
    
                resdat /= 2.
                resdat += 0.5
                result = ma.array(resdat, mask=result.mask, copy=False)                
    
            if is_scalar:
                result = result[0]            
            return result  
            
        def inverse(self, value):
            if not self.scaled():
                raise ValueError("Not invertible until scaled")
            vmin, vmax, midpoint = self.vmin, self.vmax, self.midpoint
    
            if cbook.iterable(value):
                val = ma.asarray(value)
                val = 2 * (val-0.5)  
                val[val>0]  *= abs(vmax - midpoint)
                val[val<0] *= abs(vmin - midpoint)
                val += midpoint
                return val
            else:
                val = 2 * (val - 0.5)
                if val < 0: 
                    return  val*abs(vmin-midpoint) + midpoint
                else:
                    return  val*abs(vmax-midpoint) + midpoint
    
    def format_e(n):
        a = '%e' % n
        return a.split('e')[0].rstrip('0').rstrip('.') + 'e' + a.split('e')[1]
        
    x_min,y_min = np.where(matrix_data_full == 0)
    x_max,y_max = np.where(matrix_data_full == -1)
    
    for i in xrange(n):
        for j in xrange(n):
            value = matrix_data_full[i][j]
            if value != -1 and value != 0:
                matrix_data_full[i][j] = math.log(value,2)
    
    if cutoff_max != None:
        x_cutoff_max,y_cutoff_max = np.where(matrix_data_full > cutoff_max)
        matrix_data_full[x_cutoff_max,y_cutoff_max] = cutoff_max
    if cutoff_min != None:
        x_cutoff_min,y_cutoff_min = np.where(matrix_data_full < cutoff_min)
        matrix_data_full[x_cutoff_min,y_cutoff_min] = cutoff_min
    
    threshold_max = np.max(matrix_data_full)
    threshold_min = np.min(matrix_data_full)
    matrix_data_full[x_max,y_max] = threshold_max + 1
    matrix_data_full[x_min,y_min] = threshold_min - 1
    
    if bin_size >= 1000000:
        bin_size_str = str(bin_size/1000000)
        my_bin_size = bin_size_str + ' mb'
    elif bin_size < 1000000:
        bin_size_str = str(bin_size/1000)
        my_bin_size = bin_size_str + ' kb'     
    
    plt.close("all")
    #plt.gcf().subplots_adjust(left=0.15)
    #plt.gcf().subplots_adjust(bottom=0.15)
    norm = MidPointNorm(midpoint=0)
    plt.imshow(matrix_data_full, cmap='seismic', interpolation='nearest', vmax=threshold_max, vmin=threshold_min, norm=norm)
    plt.title('O/E contact map (' + my_bin_size + ')', fontsize=12)
    plt.xlabel(chromosome + ' coordinate (bp)', fontsize=10)
    plt.ylabel(chromosome + ' coordinate (bp)', fontsize=10)
    cbar = plt.colorbar()
    cbar.cmap.set_over('black') # loci where expected enrich values are zero (log not existing)
    cbar.cmap.set_under('gray') # loci where observed values are zero (log equal to minus infinity)
    cbar.ax.set_ylabel('log2(O/E) contact counts', rotation=270, labelpad = 20)
    if coord != None:
        ticks = (np.arange(0, n, n/4) * bin_size) + coord[0]
    else:
        ticks = (np.arange(0, n, n/4) * bin_size)
    format_ticks = [format_e(i) for i in ticks.tolist()]
    plt.xticks(np.arange(0, n, n/4), format_ticks)
    plt.yticks(np.arange(0, n, n/4), format_ticks)
    plt.xticks(fontsize=8)
    plt.yticks(fontsize=8)
    plt.tick_params(axis='both', which='both', direction='out', top=False, right=False)
    plt.tight_layout()
    plt.savefig(output_path + "/" + output_filename + '.pdf', format = 'pdf', dpi = my_dpi)
    
    # Plot the histogram
    if plot_histogram == 1:        
        histogram = []
        k = 1
        for i in xrange(n):
            row = matrix_data_full[i][k:]
            for j in row:
                if j <= threshold_max and j >= threshold_min: 
                    histogram.append(j)
            k += 1        
        
        plt.close("all")
        histogram_bins = int(pow(len(histogram),0.3))
        plt.hist(histogram, bins=histogram_bins)
        plt.title('O/E contact counts distribution', fontsize=16)
        plt.xlabel('log2(O/E) contact counts', fontsize=14)
        plt.ylabel('Number of bins', fontsize=14)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.tight_layout()
        plt.savefig(output_path + "/" + output_filename + '_histogram.pdf', format = 'pdf')
    print "Done!"
    
    
def plot_chromosome_correlation_data(correlation_matrix,
                                     a_chr,
                                     bin_size,
                                     coord=None, 
                                     species='hg38',
                                     my_dpi=1000):
    """
    Plot the Pearson Correlation matrix of the "observed / expected" contact map.
    Arguments:
        correlation_matrix (str): txt file of the pearson correlation matrix generated with the function 
        "normalize_chromosome_enrich_data".
        a_chr (str): chromosome number (example for chromosome 1: '1').
        bin_size (int): bin size in bp of the contact matrix.
        coord (list): list with two integers start and end coordinate for the plot in bp.
        species (str): 'hg38' or 'mm10' or any other species label in string format.
        my_dpi (int): resolution of the contact map in dpi.
    Outputs:
        Pearson correlation matrix saved to pdf format.
    """                                       
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    import math
    from numpy import ma
    from matplotlib import cbook
    from matplotlib.colors import Normalize

    chromosome = 'chr' + a_chr
    output_path = "yaffe_tanay_" + str(bin_size)
    if not path.exists(output_path):
        os.mkdir(output_path)
    output_filename = chromosome + '_' + str(bin_size) + '_correlation_matrix'
    
    chromosomes = open(parameters['chromSizes_path'] + parameters['species'] + '.chrom.sizes', 'r')
    d_chr_dim = {}
    while True:
        try:
            line2list = next(chromosomes).split('\n')[0].split('\t')
            d_chr_dim[line2list[0]] = int(line2list[1])/bin_size
        except StopIteration:
            break
    
    end_pos = d_chr_dim[a_chr]*bin_size
    
    # Plotting the correlation data
    matrix_data_full = load_matrix(correlation_matrix)
    print "Plotting " + correlation_matrix + "..."

    # Selecting a part
    if coord != None:
        if len(coord) != 2:
            print "ERROR! Both start and end coordinate has to be declared! Otherwise leave both undeclared to plot the entire contact matrix."
            return
        else:
            print "Selecting part [" + str(coord[0]) + "-" + str(coord[1]) + "] ..."
            output_filename += "_" + str(coord[0]) + "_" + str(coord[1]) 
            start_bin = coord[0]/bin_size
            end_bin = coord[1]/bin_size
        
            if coord[0] >= coord[1]:
                print "ERROR! Start coordinate should be lower than end coordinate"
                return
            
            if start_bin >= end_bin:
                print "ERROR! Start coordinate should be much lower than the end coordinate given the bin size"
                return
        
            if coord[1] > end_pos:
                print "ERROR! End coordinate is larger than chromosome size " + str(end_pos) + " bp"
                return
    else: # Full contact matrix
        start_bin = 0
        end_bin = d_chr_dim[a_chr]
    
    matrix_data_full = matrix_data_full[start_bin:end_bin+1,start_bin:end_bin+1]
    n = len(matrix_data_full)
    
    # Heatmap plotting
    # Defining a class to generate a divergent colorbar with a custom midpoint
    class MidPointNorm(Normalize):    
        def __init__(self, midpoint=0, vmin=None, vmax=None, clip=False):
            Normalize.__init__(self,vmin, vmax, clip)
            self.midpoint = midpoint
    
        def __call__(self, value, clip=None):
            if clip is None:
                clip = self.clip
    
            result, is_scalar = self.process_value(value)
    
            self.autoscale_None(result)
            vmin, vmax, midpoint = self.vmin, self.vmax, self.midpoint
    
            if not (vmin < midpoint < vmax):
                raise ValueError("midpoint must be between maxvalue and minvalue.")       
            elif vmin == vmax:
                result.fill(0)
            elif vmin > vmax:
                raise ValueError("maxvalue must be bigger than minvalue")
            else:
                vmin = float(vmin)
                vmax = float(vmax)
                if clip:
                    mask = ma.getmask(result)
                    result = ma.array(np.clip(result.filled(vmax), vmin, vmax),
                                      mask=mask)
    
                resdat = result.data
    
                #First scale to -1 to 1 range, than to from 0 to 1.
                resdat -= midpoint            
                resdat[resdat>0] /= abs(vmax - midpoint)            
                resdat[resdat<0] /= abs(vmin - midpoint)
    
                resdat /= 2.
                resdat += 0.5
                result = ma.array(resdat, mask=result.mask, copy=False)                
    
            if is_scalar:
                result = result[0]            
            return result  
            
        def inverse(self, value):
            if not self.scaled():
                raise ValueError("Not invertible until scaled")
            vmin, vmax, midpoint = self.vmin, self.vmax, self.midpoint
    
            if cbook.iterable(value):
                val = ma.asarray(value)
                val = 2 * (val-0.5)  
                val[val>0]  *= abs(vmax - midpoint)
                val[val<0] *= abs(vmin - midpoint)
                val += midpoint
                return val
            else:
                val = 2 * (val - 0.5)
                if val < 0: 
                    return  val*abs(vmin-midpoint) + midpoint
                else:
                    return  val*abs(vmax-midpoint) + midpoint
    
    def format_e(n):
        a = '%e' % n
        return a.split('e')[0].rstrip('0').rstrip('.') + 'e' + a.split('e')[1]
    
    if bin_size >= 1000000:
        bin_size_str = str(bin_size/1000000)
        my_bin_size = bin_size_str + ' mb'
    elif bin_size < 1000000:
        bin_size_str = str(bin_size/1000)
        my_bin_size = bin_size_str + ' kb'     
    
    plt.close("all")
    #plt.gcf().subplots_adjust(left=0.15)
    #plt.gcf().subplots_adjust(bottom=0.15)
    norm = MidPointNorm(midpoint=0, vmin = -1, vmax = 1)
    plt.imshow(matrix_data_full, cmap='seismic', interpolation='nearest', norm=norm)
    plt.title('Correlation matrix (' + my_bin_size + ')', fontsize=12)
    plt.xlabel(chromosome + ' coordinate (bp)', fontsize=10)
    plt.ylabel(chromosome + ' coordinate (bp)', fontsize=10)
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('Pearson correlation', rotation=270, labelpad = 20)
    if coord != None:
        ticks = (np.arange(0, n, n/4) * bin_size) + coord[0]
    else:
        ticks = (np.arange(0, n, n/4) * bin_size)
    format_ticks = [format_e(i) for i in ticks.tolist()]
    plt.xticks(np.arange(0, n, n/4), format_ticks)
    plt.yticks(np.arange(0, n, n/4), format_ticks)
    plt.xticks(fontsize=8)
    plt.yticks(fontsize=8)
    plt.tick_params(axis='both', which='both', direction='out', top=False, right=False)
    plt.tight_layout()
    plt.savefig(output_path + "/" + output_filename + '.pdf', format = 'pdf', dpi = my_dpi)
    print "Done!"


if __name__ == '__main__':
    
    usage = 'Usage: python2.7 HiCtool_yaffe_tanay.py [-h] [options]'
    parser = OptionParser(usage = 'python2.7 %prog --action action -i input_file [options]')
    parser.add_option('--action', dest='action', type='string', help='Action to perform: normalize_fend, plot_map, normalize_enrich, plot_enrich, plot_correlation')
    parser.add_option('-i', dest='input_file', type='string', help='Input contact matrix for plotting actions, norm binning hdf5 object from HiFive for normalizing actions.')
    parser.add_option('-c', dest='chromSizes_path', type='string', help='Path to the folder chromSizes with trailing slash at the end.')
    parser.add_option('-b', dest='bin_size', type='int', help='Bin size (resolution) for the analysis.')
    parser.add_option('-s', dest='species', type='string', help='Species. It has to be one of those present under the chromSizes path.')
    parser.add_option('-p', dest='processors', type='int', default=1, help='Processors to be used to normalize the data in parallel.')
    parser.add_option('--chr', dest='chr', type='str', help='Single chromosome for normalization or plotting, or list of chromosomes between square brackets for normalization of multiple chromosomes at once.')  
    parser.add_option('--data_type', dest='data_type', type='str', help='Data type to label your data, example: observed, normalized_fend, normalized_enrich.')  
    parser.add_option('--save_obs', dest='save_obs', type='int', default=1, help='Insert 1 to save the observed data when normalizing, 0 otherwise. Default: 1.')  
    parser.add_option('--save_expect', dest='save_expect', type='int', default=0, help='Insert 1 to save the expected data when normalizing, 0 otherwise. Default: 0.')  
    parser.add_option('--coord', dest='coord', type='str', help='List of two integers with start and end coordinates for the chromosome to be plotted.')  
    parser.add_option('--my_colormap', dest='my_colormap', type='str', default='[white,red]', help='Colormap to be used to plot the contact data. You can choose among any colorbar here https://matplotlib.org/examples/color/colormaps_reference.html, or input a list of colors if you want a custom colorbar. Example: [white, red, black]. Colors can be specified also HEX format. Default: [white,red]')  
    parser.add_option('--cutoff_type', dest='cutoff_type', type='str', default='percentile', help='When plotting the contact data, to select a type of cutoff (percentile or contact_number) or plot the full range of the data (leave it empty). Default: percentile.')  
    parser.add_option('--cutoff', dest='cutoff', type='str', default='95', help='When plotting the contact data, to set a maximum cutoff on the number of contacts for the colorbar based on cutoff_type. Default: 95.')  
    parser.add_option('--cutoff_max', dest='cutoff_max', type='str', help='When plotting the enrichment data, log2 upper cutoff for the colorbar (every enrichment value above cutoff_max is plotted in red).')  
    parser.add_option('--cutoff_min', dest='cutoff_min', type='str', help='When plotting the enrichment data, log2 lower cutoff (negative value) for the colorbar (every enrichment value below cutoff_min is plotted in blue).')  
    parser.add_option('--max_color', dest='max_color', type='str', default='#460000', help='When plotting the contact data, to set the color of the bins with contact counts over "cutoff". Default: #460000.')  
    parser.add_option('--my_dpi', dest='my_dpi', type='int', default=1000, help='Resolution of the contact map in dpi. Default: 1000.')  
    parser.add_option('--plot_histogram', dest='plot_histogram', type='int', default=0, help='Insert 1 to plot the histogram of the data distribution, 0 otherwise. Default: 0.')  
    parser.add_option('--topological_domains', dest='topological_domains', type='str', help='Topological domain coordinates file (as generated from HiCtool_TAD_analysis.py) to visualize domains on the heatmap if action is "plot_map".')  
    parser.add_option('--domain_color', dest='domain_color', type='str', default='#0000ff', help='To set the color for topological domains on the heatmap. Default: #0000ff.')  
    (options, args) = parser.parse_args( )
    
    if options.action == None:
        parser.error('-h for help or provide the action command (normalize_fend, plot_map, normalize_enrich, plot_enrich, plot_correlation)!')
    else:
        pass
    if options.input_file == None:
        parser.error('-h for help or provide the input file!')
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
    if options.chr == None:
        parser.error('-h for help or provide the input chromosome (or eventually chromosomes if you are normalizing)!')
    else:
        pass
    if options.species == None:
        parser.error('-h for help or provide the species!')
    else:
        pass

    
    parameters['action'] = options.action
    parameters['input_file'] = options.input_file
    parameters['chromSizes_path'] = options.chromSizes_path
    parameters['chr'] = options.chr
    parameters['species'] = options.species
    parameters['processors'] = options.processors
    parameters['bin_size'] = options.bin_size
    parameters['data_type'] = options.data_type
    parameters['save_obs'] = options.save_obs
    parameters['save_expect'] = options.save_expect
    parameters['coord'] = options.coord
    parameters['my_colormap'] = options.my_colormap
    parameters['cutoff_type'] = options.cutoff_type
    parameters['cutoff'] = options.cutoff
    parameters['cutoff_max'] = options.cutoff_max   
    parameters['cutoff_min'] = options.cutoff_min
    parameters['max_color'] = options.max_color
    parameters['my_dpi'] = options.my_dpi
    parameters['plot_histogram'] = options.plot_histogram
    parameters['topological_domains'] = options.topological_domains
    parameters['domain_color'] = options.domain_color
    
    if parameters['species'] + ".chrom.sizes" not in os.listdir(parameters['chromSizes_path']):
        available_species = ', '.join([x.split('.')[0] for x in  os.listdir(parameters['chromSizes_path'])])
        parser.error('Wrong species inserted! Check the species spelling or insert an available species: ' + available_species + '. If your species is not listed, please contact Riccardo Calandrelli at <rcalandrelli@eng.ucsd.edu>.')
    
    chromosomes = open(parameters['chromSizes_path'] + parameters['species'] + '.chrom.sizes', 'r')
    chromosomes_list = []
    chr_dim = []
    d_chr_dim = {}
    while True:
        try:
            line2list = next(chromosomes).split('\n')[0].split('\t')
            chromosomes_list.append(line2list[0])
            chr_dim.append(int(line2list[1])/parameters['bin_size'])
            d_chr_dim[line2list[0]] = int(line2list[1])/parameters['bin_size']
        except StopIteration:
            break
    
    if parameters['action'] == 'normalize_fend':
        
        chr_list = map(str, parameters['chr'].strip('[]').split(','))
        
        if len(chr_list) > 1:
            if parameters['processors'] != None and parameters['processors'] > 1:
                print "Normalizing fend data in parallel for chromosomes " + parameters['chr'] + " using " + str(parameters['processors']) + " threads..."
                print "Start: " + strftime("%Y-%m-%d %H:%M:%S", gmtime())
                pool = Pool(processes=parameters['processors'])             
                pool.map(normalize_chromosome_fend_data, chr_list)
                print "End: " + strftime("%Y-%m-%d %H:%M:%S", gmtime())
            else:
                print "Normalizing fend data for chromosomes " + parameters['chr'] + " using a single core..."
                print "Start: " + strftime("%Y-%m-%d %H:%M:%S", gmtime())
                for c in chr_list:
                    normalize_chromosome_fend_data(c)
                print "End: " + strftime("%Y-%m-%d %H:%M:%S", gmtime())
        else:
            normalize_chromosome_fend_data(chr_list[0])
    
    elif parameters['action'] == 'normalize_enrich':
        
        chr_list = map(str, parameters['chr'].strip('[]').split(','))
        
        if len(chr_list) > 1:
            if parameters['processors'] > 1:
                print "Normalizing enrichment data in parallel for chromosomes " + parameters['chr'] + " using " + str(parameters['processors']) + " threads..."
                print "Start: " + strftime("%Y-%m-%d %H:%M:%S", gmtime())
                pool = Pool(processes=parameters['processors'])             
                pool.map(normalize_chromosome_enrich_data, chr_list)
                print "End: " + strftime("%Y-%m-%d %H:%M:%S", gmtime())
            else:
                print "Normalizing enrichment data for chromosomes " + parameters['chr'] + " using a single core..."
                print "Start: " + strftime("%Y-%m-%d %H:%M:%S", gmtime())
                for c in chr_list:
                    normalize_chromosome_enrich_data(c)
                print "End: " + strftime("%Y-%m-%d %H:%M:%S", gmtime())
        else:
            normalize_chromosome_enrich_data(chr_list[0])
    

    elif parameters['action'] == 'plot_map':
        if options.data_type == None:
            parser.error('-h for help or insert a custom label for the data type (observed, normalized_fend, normalized_enrich)!')
        else:
            pass
        
        my_cmap = map(str, parameters['my_colormap'].strip('[]').split(','))
        if len(my_cmap) == 1:
            my_cmap = my_cmap[0]
        
        chr_list = map(str, parameters['chr'].strip('[]').split(','))
        if len(chr_list) > 1:
            parser.error("Insert a single chromosome to be plotted.")
        
        if parameters['coord'] != None:
            coord = map(int, parameters['coord'].strip('[]').split(','))
        else:
            coord = None
        
        plot_chromosome_data(parameters['input_file'],
                             chr_list[0],
                             parameters['bin_size'],
                             coord,
                             parameters['species'],
                             parameters['data_type'],
                             my_cmap,
                             parameters['cutoff_type'],
                             float(parameters['cutoff']),
                             parameters['max_color'],
                             parameters['my_dpi'],
                             bool(parameters['plot_histogram']),
                             parameters['topological_domains'],
                             parameters['domain_color'])
    
    elif parameters['action'] == 'plot_enrich':
        
        chr_list = map(str, parameters['chr'].strip('[]').split(','))
        if len(chr_list) > 1:
            parser.error("Insert a single chromosome to be plotted.")
            
        if parameters['coord'] != None:
            coord = map(int, parameters['coord'].strip('[]').split(','))
        else:
            coord = None
        
        if parameters['cutoff_max'] != None:
            cutoff_max = float(parameters['cutoff_max'])
        else:
            cutoff_max = None
            
        if parameters['cutoff_min'] != None:
            cutoff_min = float(parameters['cutoff_min'])
        else:
            cutoff_min = None
            
        plot_chromosome_enrich_data(parameters['input_file'],
                                    chr_list[0],
                                    parameters['bin_size'],  
                                    coord,
                                    parameters['species'],                              
                                    cutoff_max,
                                    cutoff_min,
                                    parameters['my_dpi'],
                                    parameters['plot_histogram'])
    
    elif parameters['action']  == "plot_correlation":
        
        chr_list = map(str, parameters['chr'].strip('[]').split(','))
        if len(chr_list) > 1:
            parser.error("Insert a single chromosome to be plotted.")
            
        if parameters['coord'] != None:
            coord = map(int, parameters['coord'].strip('[]').split(','))
        else:
            coord = None
        
        plot_chromosome_correlation_data(parameters['input_file'],
                                         chr_list[0],
                                         parameters['bin_size'],
                                         coord, 
                                         parameters['species'],
                                         parameters['my_dpi'])

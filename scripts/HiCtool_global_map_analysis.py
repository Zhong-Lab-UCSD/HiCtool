# Program to analyze a HiCtool global matrix (all-by-all chromosomes) with several features:
# - Extract single contact maps from the global maps
# - Plot the global map
# - Plot a single map
# - Side-by-side plot of a contact map from different conditions for comparison.

# Usage: python2.7 HiCtool_global_map_analysis.py [-h] [options]
# Options:
#  -h, --help                show this help message and exit
#  --action                  Action to perform: extract_single_map, plot_map, plot_side_by_side_map.
#  -i INPUT_FILE             Input contact matrix file.
#  -c CHROMSIZES_PATH        Path to the folder chromSizes with trailing slash at the end.
#  -b BIN_SIZE               Bin size (resolution) of the contact matrix.
#  -s SPECIES                Species. It has to be one of those present under the chromSizes path.
#  --isGlobal                Insert 1 if the input matrix is a global matrix, 0 otherwise.  
#  --tab_sep                 Insert 1 if the input matrix is in a tab separated format, 0 if it is in compressed format.
#  --chr_row                 Chromosome or list of chromosomes between square brackets in the rows to select specific maps for extraction for plotting.
#  --chr_col                 Chromosome or list of chromosomes between square brackets in the columns to select specific maps for extraction for plotting.
#  --data_type               Data type to label your data, example: observed, normalized, etc.
#  --chr_row_coord           List of two integers with start and end coordinates for the chromosome on the rows to be plotted. It can also be a list of lists of two elements if multiple single maps are plotted.
#  --chr_col_coord           List of two integers with start and end coordinates for the chromosome on the columns to be plotted. It can also be a list of lists of two elements if multiple single maps are plotted.
#  --my_colormap             Colormap to be used to plot the data. You can choose among any colorbar here https://matplotlib.org/examples/color/colormaps_reference.html, or input a list of colors if you want a custom colorbar. Example: [white, red, black]. Colors can be specified also HEX format. Default: [white,red].  
#  --cutoff_type             To select a type of cutoff (percentile or contact_number) or plot the full range of the data (leave it empty). Default: percentile.  
#  --cutoff                  Percentile to set a maximum cutoff on the number of contacts for the colorbar. Default: 95.  
#  --max_color               To set the color of the bins with contact counts over "cutoff". Default: #460000. 
#  --my_dpi                  Resolution of the contact map in dpi. Default: 2000.  
#  --plot_histogram          Insert 1 to plot the histogram of the contact distribution of the single contact matrices, 0 otherwise. Default: 0. 
#  --topological_domains     Topological domain coordinates file (as generated from HiCtool_TAD_analysis.py) to visualize domains on the heatmap (only if a single map is selected).
#  --domain_color            To set the color for topological domains on the heatmap. Default: #0000ff.
#  --samples                 If action is "plot_side_by_side_map", insert here the samples labels between square brackets.

from optparse import OptionParser
import numpy as np
import os

parameters = {'action': None,
              'input_file': None,
              'chromSizes_path': None,
              'isGlobal': None,
              'tab_sep': None,
              'chr_row': None,
              'chr_col': None,
              'species': None,
              'bin_size': None,
              'data_type': None,
              'chr_row_coord': None,
              'chr_col_coord': None,
              'my_colormap': None,
              'cutoff_type': None,
              'cutoff': None,
              'max_color': None,
              'my_dpi': None,
              'plot_histogram': None,
              'topological_domains': None,
              'domain_color': None,
              'samples': None
              }


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


def load_matrix_rectangular(input_file, n_row, n_col):
    """
    Load an HiCtool compressed rectangular contact matrix from a txt file and parse it.
    Arguments:
        input_file (str): input file name in txt format (generated by the function 
        "save_matrix_rectangular")
        n_row (int): number of rows of the matrix.
        n_col (int): number of columns of the matrix.
    Return: 
        output_matrix: numpy array the parsed values stored in the input txt file to build a contact matrix.        
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
    print "Done!"
    return output_matrix
   

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


def extract_single_map(input_global_matrix,
                       tab_sep,
                       chr_row,
                       chr_col,
                       species='hg38',
                       bin_size=1000000,
                       data_type='observed',
                       save_output=True,
                       save_tab=False):
    """
    Extract a single contact matrix for a pair of chromosomes from the global matrix (all-by-all chromosomes).
    Arguments:
        input_global_matrix (object | str): global contact matrix. This can be passed either as
        an object of the workspace or a string of the filename saved to file.
        tab_sep (bool): if "input_global_matrix" is passed with a filename, then this boolean 
        tells if the global matrix was saved in tab separated format (True) or not (False).
        chr_row (str): chromosome in the rows of the output contact matrix.
        chr_col (str): chromosome in the columns of the output contact matrix. If chr_col is 
        equal to chr_row then the intra-chromosomal map is extracted.
        species (str): 'hg38' or 'mm10' or any other species label in string format.
        bin_size (int): bin size in bp of the contact matrix.
        data_type (str): which kind of data type you are extracting: "observed" or "normalized".
        save_output (bool): if True, save the contact matrix in HiCtool compressed txt file.
        save_tab (bool): if True, save the contact matrix in tab separated format.
    Return: 
        Contact matrix in numpy array format.
    Outputs:
        Txt file with the contact matrix in HiCtool compressed format if "save_output=True".
        Txt file with the contact matrix in tab separated format if "save_tab=True".
    """            
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
    
    d_chr_dim_inc = {}
    k=1
    for i in chromosomes_list:
        d_chr_dim_inc[i] = sum(chr_dim[:k])
        k+=1
    
    if isinstance(input_global_matrix,str):
        if tab_sep == False:
            full_matrix = load_matrix(input_global_matrix)
        else:
            full_matrix = load_matrix_tab(input_global_matrix)
    else:
        full_matrix = input_global_matrix
    
        d_chr_dim_inc = {}
    k=1
    for i in chromosomes_list:
        d_chr_dim_inc[i] = sum(chr_dim[:k])
        k+=1
    
    if isinstance(input_global_matrix,str):
        if tab_sep == False:
            full_matrix = load_matrix(input_global_matrix)
        else:
            full_matrix = load_matrix_tab(input_global_matrix)
    else:
        full_matrix = input_global_matrix
    
    if chr_row == '1':
        row_start = 0
    else:
        row_start = d_chr_dim_inc[chromosomes_list[chromosomes_list.index(chr_row)-1]]
    row_end = row_start + d_chr_dim[chr_row]
    
    if chr_col == '1':
        col_start = 0
    else:
        col_start = d_chr_dim_inc[chromosomes_list[chromosomes_list.index(chr_col)-1]]
    col_end = col_start + d_chr_dim[chr_col]
    
    output_matrix = full_matrix[row_start:row_end,col_start:col_end]
    
    if chr_row == chr_col:
        if bin_size >= 1000000:
            bin_size_str = str(bin_size/1000000)
            my_filename = 'HiCtool_' 'chr' + chr_row + '_' + bin_size_str + 'mb_' + data_type + '.txt'
        elif bin_size < 1000000:
            bin_size_str = str(bin_size/1000)
            my_filename = 'HiCtool_' 'chr' + chr_row + '_' + bin_size_str + 'kb_' + data_type + '.txt'
        if save_output == True:
            save_matrix(output_matrix, my_filename)
    else:
        if bin_size >= 1000000:
            bin_size_str = str(bin_size/1000000)
            my_filename = 'HiCtool_' 'chr' + chr_row + '_chr' + chr_col + '_' + bin_size_str + 'mb_' + data_type + '.txt'
        elif bin_size < 1000000:
            bin_size_str = str(bin_size/1000)
            my_filename = 'HiCtool_' 'chr' + chr_row + '_chr' + chr_col + '_' + bin_size_str + 'kb_' + data_type + '.txt'
        if save_output == True:
            save_matrix_rectangular(output_matrix, my_filename)
    
    if save_tab == True:
        save_matrix_tab(output_matrix, my_filename.split('.')[0] + '_tab.txt')
    
    return output_matrix
    

def plot_map(input_matrix,
             isGlobal,
             tab_sep=False,
             chr_row=None,
             chr_col=None,
             bin_size=1000000,
             chr_row_coord=None,
             chr_col_coord=None,
             data_type='observed',
             species='hg38',
             my_colormap=['white', 'red'],
             cutoff_type='percentile',
             cutoff=99,
             max_color='#460000',
             my_dpi=2000,
             plot_histogram=False,
             topological_domains=None,
             domain_color='#0000ff'):
    """
    Plot a contact map, either global or single map. To plot the global matrix leave "chr_row" and
    "chr_col" as None.
    Arguments:
        input_matrix (object | str): contact matrix. This can be passed either as
        an object of the workspace or a string of the filename saved to file.
        isGlobal (bool): set to True if you are passing a global matrix (all-by-all chromosomes), False otherwise.
        tab_sep (bool): if "input_matrix" is passed with a filename, then this boolean 
        tells if the matrix was saved in tab separated format (True) or not (False).
        chr_row (str): chromosome in the rows of the output contact matrix.
        chr_col (str): chromosome in the columns of the output contact matrix. If chr_col is 
        equal to chr_row then the intrachromosomal map is extracted.
        bin_size (int): bin size in bp of the contact matrix.
        chr_row_coord (list): list of two integers with start and end coordinates for the chromosome on the rows to be plotted (only if a single map is selected).
        chr_col_coord (list): list of two integers with start and end coordinates for the chromosome on the columns to be plotted (only if a single map is selected).
        data_type (str): which kind of data type you are extracting ("observed" or "normalized").
        species (str): 'hg38' or 'mm10' or any other species label in string format.
        my_colormap (str | list): colormap to be used to plot the data. 1) Use a string if you choose among any colorbar here 
        https://matplotlib.org/examples/color/colormaps_reference.html 2) Use a list of strings with colors if you want
        a custom colorbar. Example: ['white', 'red', 'black']. Colors can be specified also in this format: '#000000'.
        cutoff_type (str): to select a type of cutoff ('percentile' or 'contact_number') or plot the full range of the data (set the 
        parameter as None).
        cutoff (int): percentile to set a maximum cutoff on the number of contacts for the colorbar.
        max_color (str): to set the color of the bins with contact counts over "cutoff".
        my_dpi (int): resolution of the contact map in dpi.
        plot_histogram (bool): if True, plot the contact data distribution (only if a single map is selected).
        topological_domains (str | obj): topological domain coordinates to visualize domains on the heatmap. 
        They can be passed either as a txt file or object (as generated from HiCtool_TAD_analysis.py) If None, no topological domains (only if a single map is selected).
        domain_color (str): to set the color for topological domains on the heatmap. 
    Outputs:
        Heatmap saved in pdf format.
        Histogram saved in pdf format if "plot_histogram=True".
    """         
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
    import numpy as np
    import copy
    import json
    
    if bin_size > 200000:
        grid_width = 2
    elif bin_size <= 200000 and bin_size > 100000:
        grid_width = 4
    elif bin_size <= 100000 and bin_size > 50000:
        grid_width = 8
    elif bin_size <= 50000:
        grid_width = 16
    
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
    
    d_chr_dim_inc = {}
    k=1
    for i in chromosomes_list:
        d_chr_dim_inc[i] = sum(chr_dim[:k])
        k+=1
    
    d_chr_label_pos = {} # label position for chromosomes in the global matrix
    k = 0 # to consider the pixel occupied by the grid added after
    for i in chromosomes_list:
        d_chr_label_pos[i] = d_chr_dim_inc[i] - d_chr_dim[i]/2 + k
        k += grid_width
    
    label_pos = []
    label_name = []
    for i in chromosomes_list:
        label_pos.append(d_chr_label_pos[i])
        label_name.append('chr' + i)
    label_pos = np.array(label_pos)
    label_name = tuple(label_name)
    
    ### Plot global heatmap
    if chr_row == None and chr_col == None and isGlobal == True:
        if bin_size >= 1000000:
            bin_size_str = str(bin_size/1000000) + 'mb'
            my_filename = 'HiCtool_' + bin_size_str + '_' + data_type + "_heatmap"
        elif bin_size < 1000000:
            bin_size_str = str(bin_size/1000) + 'kb'
            my_filename = 'HiCtool_' + bin_size_str + '_' + data_type + "_heatmap"
        
        if isinstance(input_matrix,str):
            if tab_sep == False:
                matrix_data_full = load_matrix(input_matrix)
            else:
                matrix_data_full = load_matrix_tab(input_matrix)
        else:
            matrix_data_full = copy.deepcopy(input_matrix)
        
        print "Plotting the global matrix..."
        # Adding grid to separate chromosomes
        k=0
        for i in chromosomes_list[:-1]:
            for j in range(grid_width):
                matrix_data_full = np.insert(matrix_data_full, d_chr_dim_inc[i]+k, -1, axis=1)
                matrix_data_full = np.insert(matrix_data_full, d_chr_dim_inc[i]+k, -1, axis=0)
            k += grid_width
        
        row = np.shape(matrix_data_full)[0]
        col = np.shape(matrix_data_full)[1]
        
        output_vect = np.reshape(matrix_data_full,row*col,1)
        negative_indexes = np.where(output_vect==-1)
        output_vect[negative_indexes] = 0
        non_zero = np.nonzero(output_vect)
        
        if isinstance(my_colormap, list):
            my_cmap = LinearSegmentedColormap.from_list('mycmap', my_colormap)
        elif isinstance(my_colormap, str):
            my_cmap = my_colormap
        
        plt.close("all")
        
        if cutoff_type == 'percentile':
            perc = np.percentile(output_vect[non_zero[0]],cutoff)
            plt.imshow(matrix_data_full, cmap=my_cmap, interpolation='nearest', vmax=perc , vmin=0)
            cbar = plt.colorbar(extend='max')
            cbar.cmap.set_over(max_color)
        elif cutoff_type == 'contact':
            perc = cutoff 
            plt.imshow(matrix_data_full, cmap=my_cmap, interpolation='nearest', vmax=perc , vmin=0)
            cbar = plt.colorbar(extend='max')
            cbar.cmap.set_over(max_color)
        elif cutoff_type == None:
            plt.imshow(matrix_data_full, cmap=my_cmap, interpolation='nearest', vmin=0)
            cbar = plt.colorbar()
            
        cbar.cmap.set_under('black')   
        cbar.ax.set_ylabel(data_type + ' contact counts', rotation=270, labelpad=20)
        plt.title(data_type + ' contact map (' + bin_size_str + ')', fontsize=12)
        plt.xticks(label_pos, label_name, rotation='vertical', fontsize = 6)
        plt.yticks(label_pos, label_name, fontsize = 6)
        plt.tick_params(axis='both', which='both', length=0)
        plt.savefig(my_filename + '.pdf', format = 'pdf', dpi=my_dpi)
        print "Done!"
    
    ### Plot a single heatmap
    if (chr_row != None and chr_col == None) or (chr_row == None and chr_col != None):
        print "ERROR! Chromosomes on the rows and on the columns have to be declared."
        return
    if chr_row != None and chr_col != None:
        chr_row_list = map(str, chr_row.strip('[]').split(','))
        chr_col_list = map(str, chr_col.strip('[]').split(','))
        if len(chr_row_list) != len(chr_row_list):
            print 'ERROR! chr_row and chr_col should be of the same length!'
            return
        matrix_data_full_list = []
        
        if topological_domains != None:
            topological_domains_list = map(str, topological_domains.strip('[]').split(','))
            if len(topological_domains_list) != len(chr_row_list):
                print "ERROR! Insert in topological domains the same number of elements than chr_row (or chr_col). Leave empty slots where you do not wish to plot topological domains."
                return
            
        if chr_row_coord != None or chr_col_coord != None:
            if chr_row_coord == None:
                print "ERROR! Coordinates on the chromosomes on the rows should be declared."
                return
            if chr_col_coord == None:
                print "ERROR! Coordinates on the chromosomes on the cols should be declared."
                return
        
        # Plotting one of more single heatmaps from the global map
        if isGlobal == True:
            if tab_sep == True:
                matrix_global = load_matrix_tab(input_matrix)
            elif tab_sep == False:
                matrix_global = load_matrix(input_matrix)
            
            for c_row, c_col in zip(chr_row_list, chr_col_list):
                matrix_data_full_list.append(extract_single_map(matrix_global,tab_sep,c_row,c_col,species,bin_size,data_type,True,False))
        # Plotting one single heatmaps from single file
        else:
            if isinstance(input_matrix,str):
                if tab_sep == True:
                    matrix_data_full = load_matrix_tab(input_matrix)
                else:
                    if chr_row_list[0] == chr_col_list[0]:
                        matrix_data_full = load_matrix(input_matrix)
                    else:
                        matrix_data_full = load_matrix_rectangular(input_matrix)
            else:
                matrix_data_full = copy.deepcopy(input_matrix)
            
            matrix_data_full_list.append(matrix_data_full)
        
        m_index = 0 # to select each single matrix
        for c_row, c_col in zip(chr_row_list, chr_col_list):
            print "Plotting chr" + c_row + "xchr" + c_col + " contact map..."
            chromosome_row = 'chr' + c_row
            chromosome_col = 'chr' + c_col        
            matrix_data_full = matrix_data_full_list[m_index]
            
            if bin_size >= 1000000:
                bin_size_str = str(bin_size/1000000) + 'mb'
                my_filename = 'HiCtool_' + chromosome_row + '_' + chromosome_col + '_' + bin_size_str + '_' + data_type + "_heatmap"
            elif bin_size < 1000000:
                bin_size_str = str(bin_size/1000) + 'kb'
                my_filename = 'HiCtool_' + chromosome_row + '_' + chromosome_col + '_' + bin_size_str + '_' + data_type + "_heatmap"
            
            # Update matrix values to plot topological domains
            if topological_domains != None:
                if topological_domains_list[m_index] != '':
                    if c_row != c_col:
                        print "WARNING! To plot topological domains the matrix in position " + str(m_index) + " should be intrachromosomal. Topological domains for this matrix are skipped."
                    else:
                        if isinstance(topological_domains_list[m_index], str):
                            domains = load_topological_domains(topological_domains_list[m_index])
                        else:
                            domains = topological_domains_list[m_index]
                        my_filename = my_filename + '_domains'
                        diag_index = np.diag_indices(len(matrix_data_full))
                        for domain in domains:
                            temp_start = domain[0]/bin_size
                            temp_end = domain[1]/bin_size
                            matrix_data_full[temp_start,temp_start:temp_end] = -1
                            matrix_data_full[temp_start:temp_end,temp_end-1] = -1
                            matrix_data_full[(diag_index[0][temp_start:temp_end],diag_index[1][temp_start:temp_end])] = -1
                
        
            # Selecting a part of a single heatmap
            if chr_row_coord != None and chr_col_coord != None:
                chr_row_coord_list = json.loads(chr_row_coord)
                chr_col_coord_list = json.loads(chr_col_coord)
                
                if len(chr_row_list) > 1 and len(chr_col_list) > 1:
                    chr_row_coord_temp = chr_row_coord_list[m_index]
                    chr_col_coord_temp = chr_col_coord_list[m_index]
                else:
                    chr_row_coord_temp = chr_row_coord_list
                    chr_col_coord_temp = chr_col_coord_list
                
                if len(chr_row_coord_temp) == 1 or len(chr_col_coord_temp) == 1:
                    print "ERROR! Start and end coordinate for each chromosome in position " + str(m_index) + " should be declared."
                    return
                if len(chr_row_coord_temp) > 2 or len(chr_col_coord_temp) > 2:
                    print "ERROR! Only two coordinates (start and end) for each chromosome in position " + str(m_index) + " should be declared."
                    return
                if len(chr_row_coord_temp) == 2 and len(chr_col_coord_temp) == 2:
                    print "Selecting part of heatmap: [" + str(chr_row_coord_temp[0]) + ":" + str(chr_row_coord_temp[1]) + "-" + str(chr_col_coord_temp[0]) + ":" + str(chr_col_coord_temp[1]) + "] ..."
                    chr_row_bin = map(lambda x: x/bin_size, chr_row_coord_temp)
                    chr_col_bin = map(lambda x: x/bin_size, chr_col_coord_temp)
            
                    if chr_row_coord_temp[0] >= chr_row_coord_temp[1] or chr_col_coord_temp[0] >= chr_col_coord_temp[1]:
                        print "ERROR! Start coordinate for chromosomes in position " + str(m_index) + " should be lower than end coordinate"
                        return
                    if chr_row_bin[0] >= chr_row_bin[1] or chr_col_bin[0] >= chr_col_bin[1]:
                        print "ERROR! Start coordinate for chromosome in position " + str(m_index) + " should be much lower than the end coordinate given the bin size"
                        return
                    if chr_row_bin[1] > d_chr_dim[c_row]:
                        print "ERROR! End coordinate of the chromosome on the rows in position " + str(m_index) + " should be lower than the chromosome size"
                        return
                    if chr_col_bin[1] > d_chr_dim[c_col]:
                        print "ERROR! End coordinate of the chromosome on the columns in position " + str(m_index) + " should be lower than the chromosome size"
                        return
                    
                    matrix_data_full = matrix_data_full[chr_row_bin[0]:chr_row_bin[1]+1,chr_col_bin[0]:chr_col_bin[1]+1]
                         
                
            row = np.shape(matrix_data_full)[0]
            col = np.shape(matrix_data_full)[1]
            
            output_vect = np.reshape(matrix_data_full,row*col,1)
            non_zero = np.nonzero(output_vect)
            
            if isinstance(my_colormap, list):
                my_cmap = LinearSegmentedColormap.from_list('mycmap', my_colormap)
            elif isinstance(my_colormap, str):
                my_cmap = my_colormap
            
            def format_e(n):
                a = '%e' % n
                return a.split('e')[0].rstrip('0').rstrip('.') + 'e' + a.split('e')[1]        
            
            plt.close("all")
            plt.gcf().subplots_adjust(left=0.15)
            plt.gcf().subplots_adjust(bottom=0.15)
            
            if cutoff_type == 'percentile':
                perc = np.percentile(output_vect[non_zero[0]],cutoff)
                cutoff_flag = False
            elif cutoff_type == 'contact':
                perc = cutoff
                if perc > np.max(matrix_data_full):
                    print "WARNING! The contact cutoff you have inserted is above the maximum value of the heatmap --> the colorbar will span the values 0-cutoff but the actual contacts on the heatmap will be below the cutoff. Please insert a cutoff below the maximum value " + str(np.max(matrix_data_full)) + " if you wish to put an upper contact cutoff."
                    cutoff_flag = True
                else:
                    cutoff_flag = False
                
            if cutoff_type != None: 
                if topological_domains != None:
                    if topological_domains_list[m_index] == '':
                        plt.imshow(matrix_data_full, cmap=my_cmap, interpolation='nearest', vmax=perc)
                        if cutoff_flag == False:
                            cbar = plt.colorbar(extend='max')
                            cbar.cmap.set_over(max_color)
                        else:
                            cbar = plt.colorbar()
                    else:
                        plt.imshow(matrix_data_full, cmap=my_cmap, interpolation='nearest', vmax=perc, vmin=0)
                        if cutoff_flag == False:
                            cbar = plt.colorbar(extend='max')
                            cbar.cmap.set_over(max_color)
                        else:
                            cbar = plt.colorbar()
                        cbar.cmap.set_under(domain_color)
                else:
                    plt.imshow(matrix_data_full, cmap=my_cmap, interpolation='nearest', vmax=perc)
                    if cutoff_flag == False:
                            cbar = plt.colorbar(extend='max')
                            cbar.cmap.set_over(max_color)
                    else:
                        cbar = plt.colorbar()
            else:
                if topological_domains != None:
                    if topological_domains_list[m_index] == '':
                        plt.imshow(matrix_data_full, cmap=my_cmap, interpolation='nearest')
                        cbar = plt.colorbar()
                    else:
                        plt.imshow(matrix_data_full, cmap=my_cmap, interpolation='nearest', vmin=0)
                        cbar = plt.colorbar()
                        cbar.cmap.set_under(domain_color)
                else:
                    plt.imshow(matrix_data_full, cmap=my_cmap, interpolation='nearest')
                    cbar = plt.colorbar()
            
            plt.title(data_type + ' contact map (' + bin_size_str + ')', fontsize=12)
            cbar.ax.set_ylabel(data_type + ' contact counts', rotation=270, labelpad=20)
            plt.ylabel(chromosome_row + ' coordinate (bp)', fontsize=10)
            plt.xlabel(chromosome_col + ' coordinate (bp)', fontsize=10)
            if chr_row_coord != None and chr_col_coord != None:
                if chr_row_coord_temp != [] and chr_col_coord_temp != []:
                    ticks_row = (np.arange(0, row, row/4) * bin_size) + chr_row_coord_temp[0]
                    ticks_col = (np.arange(0, col, col/4) * bin_size) + chr_col_coord_temp[0]
                    format_ticks_row = [format_e(i) for i in ticks_row.tolist()]
                    format_ticks_col = [format_e(i) for i in ticks_col.tolist()]
                    plt.yticks(np.arange(0, row, row/4), format_ticks_row)
                    plt.xticks(np.arange(0, col, col/4), format_ticks_col)
                    if chromosome_row != chromosome_col:
                        my_filename += '_' + str(chr_row_coord_temp[0]) + "_" + str(chr_row_coord_temp[1]) + "-" + str(chr_col_coord_temp[0]) + "_" + str(chr_col_coord_temp[1])
                    else:
                        my_filename += '_' + str(chr_row_coord_temp[0]) + "_" + str(chr_row_coord_temp[1])
                else:
                    ticks_row = (np.arange(0, row, row/4) * bin_size)
                    ticks_col = (np.arange(0, col, col/4) * bin_size)
                    format_ticks_row = [format_e(i) for i in ticks_row.tolist()]
                    format_ticks_col = [format_e(i) for i in ticks_col.tolist()]
                    plt.yticks(np.arange(0, row, row/4), format_ticks_row)
                    plt.xticks(np.arange(0, col, col/4), format_ticks_col)
            else:
                ticks_row = (np.arange(0, row, row/4) * bin_size)
                ticks_col = (np.arange(0, col, col/4) * bin_size)
                format_ticks_row = [format_e(i) for i in ticks_row.tolist()]
                format_ticks_col = [format_e(i) for i in ticks_col.tolist()]
                plt.yticks(np.arange(0, row, row/4), format_ticks_row)
                plt.xticks(np.arange(0, col, col/4), format_ticks_col)
            plt.tick_params(axis='both', which='both', direction='out', top=False, right=False)
            plt.xticks(fontsize=8)
            plt.yticks(fontsize=8)
            plt.savefig(my_filename + '.pdf', format = 'pdf', dpi=my_dpi)
            
            # Plot of the histogram
            if plot_histogram:
                histogram = []
                if c_row == c_col:
                    n = len(matrix_data_full)
                    k = 1
                    for i in xrange(n):
                        row = matrix_data_full[i][k:]
                        for j in row:
                            histogram.append(j)
                        k += 1
                else:
                    histogram = matrix_data_full.reshape((1,row*col)).tolist()[0]
                    
                plt.close("all")
                histogram_bins = int(pow(len(histogram),0.3))
                plt.hist(histogram, bins=histogram_bins)
                plt.title(data_type + ' contact counts distribution', fontsize=18)
                plt.xlabel(data_type + ' contact counts', fontsize=16)
                plt.ylabel('Number of bins', fontsize=16)
                plt.xticks(fontsize=16)
                plt.yticks(fontsize=16)
                plt.tight_layout()
                plt.savefig(my_filename + '_histogram.pdf', format = 'pdf')
            
            m_index += 1
            print "Done!"
        

def plot_side_by_side_map(inputFiles,
                          tab_sep,
                          chr_row,
                          chr_col,
                          samples,
                          bin_size,
                          data_type,
                          species='hg38',
                          my_colormap=['white','red'],
                          cutoff_type='percentile',
                          cutoff=95,
                          max_color='#460000',
                          my_dpi=2000):
    """
    Plot the same contact maps from different samples on a side-by-side view, stacked on each other.
    Arguments:
        inputFiles (list): List of global contact matrices for each sample.
        tab_sep (bool): True if the input global matrices are in tab separated format, False otherwise.
        chr_row (list): list of chromosomes in the rows of each single contact matrix.
        chr_col (list): list of chromosomes in the columns of each single contact matrix.
        samples (list): list of the samples labels.
        bin_size (int): bin size in bp of the contact matrix.
        data_type (str): which kind of data type you are extracting ("observed" or "normalized").
        species (str): species label in string format.
        my_colormap (str | list): colormap to be used to plot the data. 1) Use a string if you choose among any colorbar here 
        https://matplotlib.org/examples/color/colormaps_reference.html 2) Use a list of strings with colors if you want
        a custom colorbar. Example: ['white', 'red', 'black']. Colors can be specified also in this format: '#000000'.
        cutoff_type (str): to select a type of cutoff ('percentile' or 'contact_number') or plot the full range of the data (set the 
        parameter as None).
        cutoff (int): percentile to set a maximum cutoff on the number of contacts for the colorbar.
        max_color (str): to set the color of the bins with contact counts over "cutoff".
        my_dpi (int): resolution of the contact map in dpi.
    Output:
        Heatmap saved in pdf format.
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
    import numpy as np
    #import json
    
    # Different width for the grid to separate contact maps at different resolutions in order to be well visualized in the plots
    if bin_size > 200000:
        grid_width = 2
    elif bin_size <= 200000 and bin_size > 100000:
        grid_width = 4
    elif bin_size <= 100000 and bin_size > 50000:
        grid_width = 8
    elif bin_size <= 50000:
        grid_width = 16
    
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
            
    d_chr_dim_inc = {}
    k=1
    for i in chromosomes_list:
        d_chr_dim_inc[i] = sum(chr_dim[:k])
        k+=1
    
    label_pos = []
    label_name = []
    last_pos = 0
    n = 0
    for chr_r,chr_c in zip(chr_row,chr_col):
        if n == 0:
            label_pos.append(d_chr_dim[chr_r]/2)
            last_pos += d_chr_dim[chr_r] + grid_width # to consider the grid
            n+=1
        else:
            label_pos.append(last_pos + d_chr_dim[chr_r]/2)
            last_pos += d_chr_dim[chr_r] + grid_width
             
        if chr_r == chr_c:
            label_name.append('chr' + chr_r)
        else:
            label_name.append('chr' + chr_r + '-chr' + chr_c)
    label_pos = np.array(label_pos)
    label_name = tuple(label_name)
    
    if len(my_colormap) > 1:
        my_cmap = LinearSegmentedColormap.from_list('mycmap', my_colormap)
    elif len(my_colormap) == 1:
        my_cmap = my_colormap[0]
    
    if bin_size >= 1000000:
        bin_size_str = str(bin_size/1000000) + 'mb'
        my_filename = 'HiCtool_' + '_'.join(samples) + '_' + bin_size_str + '_' + data_type + "_heatmap"
    elif bin_size < 1000000:
        bin_size_str = str(bin_size/1000) + 'kb'
        my_filename = 'HiCtool_' + '_'.join(samples) + '_' + bin_size_str + '_' + data_type + "_heatmap"
    
#    if topological_domains != None:
#        topological_domains_list = json.loads(topological_domains) # list of lists
#        if len(topological_domains_list) != len(chr_row_list):
#            print "Insert in topological domains the same number of elements than chr_row (or chr_col). Leave empty lists where you do not wish to plot topological domains."
#            return
    
    # From global heatmap
    # Load each global map at every time point into a dictionary
    inputFiles_dict = dict()
    for i in range(len(inputFiles)):
        if tab_sep == False:
            matrix_data_full = load_matrix(inputFiles[i])
        else:
            matrix_data_full = load_matrix_tab(inputFiles[i])
        inputFiles_dict[samples[i]] = matrix_data_full

    time_steps = np.linspace(0,len(chr_row),11).astype(int).tolist() # to print percentage of completion to console


    print "(1/2) Building the matrix..."
    counter = 0 # to print percentage of completion to console
    init_counter = 0 # counter to initialize the output full matrix (if equal 1) or not
    n_col_list = [] # to save the uumber of columns in each row
    n_col_max = max([d_chr_dim[x] for x in chr_col]) * len(samples) + (len(samples)-1)*grid_width # to consider the grid
    for key, value in d_chr_dim.items():
        if value == max([d_chr_dim[x] for x in chr_col]):
            chr_col_max = key # bigger chromosomes in the columns
    output_full_matrix = np.zeros((1,n_col_max)) # initialize matrix so I can concatenate already a line where the chromosome in the columns is not the biggest
    #row_index = 0 # to select the topological domains for a row
    for i,j in zip(chr_row, chr_col):
        init_counter += 1
        line_dict = dict() # to save the contact matrices per each line
        #col_index = 0 # to select each topological domain file of a row per time point
        for k in samples:
            # Extract the single map
            if i == '1':
                row_start = 0
            else:
                row_start = d_chr_dim_inc[chromosomes_list[chromosomes_list.index(i)-1]]
            row_end = row_start + d_chr_dim[i]
            
            if j == '1':
                col_start = 0
            else:
                col_start = d_chr_dim_inc[chromosomes_list[chromosomes_list.index(j)-1]]
            col_end = col_start + d_chr_dim[j]
            
            input_matrix_array = inputFiles_dict[k] # global map at the k time point
            output_matrix = input_matrix_array[row_start:row_end,col_start:col_end]
            
            # Update matrix values to plot topological domains
#            if topological_domains != None:
#                if topological_domains_list[row_index] != []:
#                    if i != j:
#                        print "WARNING! To plot topological domains the matrices in position " + str(row_index) + " should be intrachromosomal."
#                    else:
#                        domains = load_topological_domains(topological_domains_list[row_index][col_index])
#                        #my_filename = my_filename + '_domains'
#                        diag_index = np.diag_indices(len(output_matrix))
#                        for domain in domains:
#                            temp_start = domain[0]/bin_size
#                            temp_end = domain[1]/bin_size
#                            output_matrix[temp_start,temp_start:temp_end] = -1
#                            output_matrix[temp_start:temp_end,temp_end-1] = -1
#                            output_matrix[(diag_index[0][temp_start:temp_end],diag_index[1][temp_start:temp_end])] = -1
#            
            line_dict[k] = output_matrix # fill in the dictionary with each matrix per time point for this line
            #col_index += 1
        
        # Build the line
        matrix_line = line_dict[samples[0]] # initialize the line with the first matrix
        n_row = np.shape(matrix_line)[0]
        sep_col = np.zeros((n_row,grid_width))-1
        if j == chr_col_max: # I do not add the last sep_col since it's going till the full width of the heatmap
            for t in range(len(samples))[1:]:
                matrix_line = np.concatenate((matrix_line,sep_col,line_dict[samples[t]]), axis=1)
        else: # for the other chromosomes that are shorter I add the last sep_col
            for t in range(len(samples))[1:]:
                if t != len(samples)-1:
                    matrix_line = np.concatenate((matrix_line,sep_col,line_dict[samples[t]]), axis=1)
                else: # another sep_col is added at the end
                    matrix_line = np.concatenate((matrix_line,sep_col,line_dict[samples[t]],sep_col), axis=1)
        
        # Attach the row to the output full matrix
        if init_counter == 1: # initialize the output full matrix
            n_col = np.shape(matrix_line)[1]
            diff_col = n_col_max - n_col
            matrix_line_full = np.concatenate((matrix_line,np.zeros((n_row,diff_col))), axis=1)
            output_full_matrix = np.concatenate((output_full_matrix,matrix_line_full), axis=0)
            n_col_list.append(n_col)
        else:
            n_col = np.shape(matrix_line)[1]
            diff_col = n_col_max - n_col
            matrix_line_full = np.concatenate((matrix_line,np.zeros((n_row,diff_col))), axis=1)
            if n_col < n_col_list[-1]: # the chromosome length is smaller than the previous chromosomes
                sep_row = np.concatenate((np.zeros((grid_width,n_col_list[-1]))-1,np.zeros((grid_width,n_col_max-n_col_list[-1]))), axis=1)
            else:
                sep_row = np.concatenate((np.zeros((grid_width,n_col))-1,np.zeros((grid_width,diff_col))), axis=1)
            n_col_list.append(n_col)
            output_full_matrix = np.concatenate((output_full_matrix,sep_row,matrix_line_full), axis=0)
        
        counter += 1
        for i in range(len(time_steps)):
            if counter == time_steps[i]:
                print str(i*10) + '% completed.'
        #row_index += 1

    print "(1/2) Done!"

    print "(2/2) Plotting..."
    
    output_full_matrix = output_full_matrix[1::]
    row = np.shape(output_full_matrix)[0]
    col = np.shape(output_full_matrix)[1]
    
    output_vect = np.reshape(output_full_matrix,row*col,1)
    negative_indexes = np.where(output_vect==-1)
    output_vect[negative_indexes] = 0
    non_zero = np.nonzero(output_vect)
    
    plt.close("all")
    
    if cutoff_type == 'percentile':
        perc = np.percentile(output_vect[non_zero[0]],cutoff)
        plt.imshow(output_full_matrix, cmap=my_cmap, interpolation='nearest', vmax=perc , vmin=0)
        cbar = plt.colorbar(extend='max')
        cbar.cmap.set_over(max_color)
    elif cutoff_type == 'contact':
        perc = cutoff 
        plt.imshow(output_full_matrix, cmap=my_cmap, interpolation='nearest', vmax=perc , vmin=0)
        cbar = plt.colorbar(extend='max')
        cbar.cmap.set_over(max_color)
    elif cutoff_type == None:
        plt.imshow(output_full_matrix, cmap=my_cmap, interpolation='nearest', vmin=0)
        cbar = plt.colorbar()
    
    cbar.cmap.set_under('black')   
    cbar.ax.set_ylabel(data_type + ' contact counts', rotation=270, labelpad=20)
    plt.title(data_type + ' map (' + bin_size_str + ')', fontsize=12)
    sample_pos_k = [0.5,1.5]
    if len(samples) > 2:
        for i in xrange(len(samples[2:])):
            sample_pos_k.append(sample_pos_k[-1]+1)
    sample_pos = np.array([int(x*d_chr_dim[chr_col_max]) for x in sample_pos_k])
    plt.xticks(sample_pos, tuple(samples), fontsize = 6, rotation=45)
    plt.yticks(label_pos, label_name, fontsize = 6)
    plt.tick_params(axis='both', which='both', length=0)
    plt.savefig(my_filename + '.pdf', format = 'pdf', dpi=my_dpi)


    print "(2/2) Done!"
    
    
if __name__ == '__main__':
    
    usage = 'Usage: python2.7 HiCtool_global_map_analysis.py [-h] [options]'
    parser = OptionParser(usage = 'python2.7 %prog --action action -i input_file [options]')
    parser.add_option('--action', dest='action', type='string', help='Action to perform: extract_single_map, plot_map, plot_side_by_side_map')
    parser.add_option('-i', dest='input_file', type='string', help='Input contact matrix file.')
    #parser.add_option('-o', dest='output_path', type='string', help='Path to save the output files with the trailing slash in the end.')
    parser.add_option('-c', dest='chromSizes_path', type='string', help='Path to the folder chromSizes with trailing slash at the end.')
    parser.add_option('-b', dest='bin_size', type='int', help='Bin size (resolution) of the contact matrix.')
    parser.add_option('-s', dest='species', type='string', help='Species. It has to be one of those present under the chromSizes path.')  
    parser.add_option('--isGlobal', dest='isGlobal', type='int', help='Insert 1 if the input matrix is a global matrix, 0 otherwise.')  
    parser.add_option('--tab_sep', dest='tab_sep', type='int', help='Insert 1 if the input matrix is in a tab separated format, 0 if it is in compressed format.')  
    parser.add_option('--chr_row', dest='chr_row', type='str', help='Chromosome or list of chromosomes between square brackets in the rows to select specific maps for extraction for plotting.')  
    parser.add_option('--chr_col', dest='chr_col', type='str', help='Chromosome or list of chromosomes between square brackets in the columns to select specific maps for extraction for plotting.')  
    parser.add_option('--data_type', dest='data_type', type='str', help='Data type to label your data, example: observed, normalized, etc.')  
    parser.add_option('--chr_row_coord', dest='chr_row_coord', type='str', help='List of two integers with start and end coordinates for the chromosome on the rows to be plotted. It can also be a list of lists of two elements if multiple single maps are plotted.')  
    parser.add_option('--chr_col_coord', dest='chr_col_coord', type='str', help='List of two integers with start and end coordinates for the chromosome on the columns to be plotted. It can also be a list of lists of two elements if multiple single maps are plotted.')  
    parser.add_option('--my_colormap', dest='my_colormap', type='str', default='[white,red]', help='Colormap to be used to plot the data. You can choose among any colorbar here https://matplotlib.org/examples/color/colormaps_reference.html, or input a list of colors if you want a custom colorbar. Example: [white, red, black]. Colors can be specified also HEX format. Default: [white,red]')  
    parser.add_option('--cutoff_type', dest='cutoff_type', type='str', default='percentile', help='To select a type of cutoff (percentile or contact_number) or plot the full range of the data (leave it empty). Default: percentile.')  
    parser.add_option('--cutoff', dest='cutoff', type='str', default='95', help='To set a maximum cutoff on the number of contacts for the colorbar based on cutoff_type. Default: 95.')  
    parser.add_option('--max_color', dest='max_color', type='str', default='#460000', help='To set the color of the bins with contact counts over "cutoff". Default: #460000.')  
    parser.add_option('--my_dpi', dest='my_dpi', type='int', default=2000, help='Resolution of the contact map in dpi. Default: 2000.')  
    parser.add_option('--plot_histogram', dest='plot_histogram', type='int', default=0, help='Insert 1 to plot the histogram of the contact distribution of the single contact matrices, 0 otherwise. Default: 0.')  
    parser.add_option('--topological_domains', dest='topological_domains', type='str', help='Topological domain coordinates file (as generated from HiCtool_TAD_analysis.py) to visualize domains on the heatmap (only if a single map is selected).')  
    parser.add_option('--domain_color', dest='domain_color', type='str', default='#0000ff', help='To set the color for topological domains on the heatmap. Default: #0000ff.')  
    parser.add_option('--samples', dest='samples', type='str', help='If action is "plot_side_by_side_map", insert here the samples labels between square brackets.')  
    (options, args) = parser.parse_args( )
    
    if options.action == None:
        parser.error('-h for help or provide the action command (extract_single_map, plot_map, plot_side_by_side_map)!')
    else:
        pass
    if options.input_file == None:
        parser.error('-h for help or provide the input contact matrix!')
    else:
        pass
#    if options.output_path == None:
#        parser.error('-h for help or provide the output path!')
#    else:
#        pass
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
    if options.tab_sep == None:
        parser.error('-h for help or insert 1 if the contact matrix is in tab separated format, 0 otherwise!')
    else:
        pass
    if options.data_type == None:
        parser.error('-h for help or insert a custom label for the data type (observed, normalized, etc.)!')
    else:
        pass

    
    parameters['action'] = options.action
    parameters['input_file'] = options.input_file
    #parameters['output_path'] = options.output_path
    parameters['chromSizes_path'] = options.chromSizes_path
    parameters['isGlobal'] = options.isGlobal
    parameters['tab_sep'] = options.tab_sep
    parameters['chr_row'] = options.chr_row
    parameters['chr_col'] = options.chr_col
    parameters['species'] = options.species
    parameters['bin_size'] = options.bin_size
    parameters['data_type'] = options.data_type
    parameters['chr_row_coord'] = options.chr_row_coord
    parameters['chr_col_coord'] = options.chr_col_coord
    parameters['my_colormap'] = options.my_colormap
    parameters['cutoff_type'] = options.cutoff_type
    parameters['cutoff'] = options.cutoff
    parameters['max_color'] = options.max_color
    parameters['my_dpi'] = options.my_dpi
    parameters['plot_histogram'] = options.plot_histogram
    parameters['topological_domains'] = options.topological_domains
    parameters['domain_color'] = options.domain_color
    parameters['samples'] = options.samples
    
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
    
    if parameters['action'] == 'extract_single_map':
        if options.chr_row == None:
            parser.error('-h for help or insert a chromosome or a list of chromosomes for the rows of the contact matrices!')
        else:
            pass
        if options.chr_col == None:
            parser.error('-h for help or insert a chromosome or a list of chromosomes for the columns of the contact matrices!')
        else:
            pass
        
        chr_row_list = map(str, parameters['chr_row'].strip('[]').split(','))
        chr_col_list = map(str, parameters['chr_col'].strip('[]').split(','))
        
        if len(chr_row_list) != len(chr_row_list):
            parser.error('chr_row and chr_col should be of the same length!')
        
        if bool(parameters['tab_sep']) == True:
            matrix_global = load_matrix_tab(parameters['input_file'])
        elif bool(parameters['tab_sep']) == False:
            matrix_global = load_matrix(parameters['input_file'])
        
        for c_row, c_col in zip(chr_row_list, chr_col_list):
            print "Extracting chr" + c_row + "xchr" + c_col + " contact matrix..."
            extract_single_map(matrix_global,
                               bool(parameters['tab_sep']),
                               c_row, 
                               c_col,
                               parameters['species'],
                               parameters['bin_size'],
                               parameters['data_type'],
                               True, 
                               True)
            print "Done!"
    
    elif parameters['action'] == 'plot_map':
        if options.isGlobal == None:
            parser.error('-h for help or insert 1 if the contact matrix is a global all-by-all chromosomes, 0 if it is a single contact matrix!')
        else:
            pass
        
        my_cmap = map(str, parameters['my_colormap'].strip('[]').split(','))
        if len(my_cmap) == 1:
            my_cmap = my_cmap[0]
        
        plot_map(parameters['input_file'],
                 bool(parameters['isGlobal']),
                 bool(parameters['tab_sep']),
                 parameters['chr_row'],
                 parameters['chr_col'],
                 parameters['bin_size'],
                 parameters['chr_row_coord'],
                 parameters['chr_col_coord'],
                 parameters['data_type'],
                 parameters['species'],
                 my_cmap,
                 parameters['cutoff_type'],
                 float(parameters['cutoff']),
                 parameters['max_color'],
                 parameters['my_dpi'],
                 bool(parameters['plot_histogram']),
                 parameters['topological_domains'],
                 parameters['domain_color'])
    
    elif parameters['action'] == 'plot_side_by_side_map':
        if options.samples == None:
            parser.error('-h for help or insert the label for each data point between square brackets!')
        else:
            pass
        if options.chr_row == None:
            parser.error('-h for help or insert a chromosome or a list of chromosomes for the rows of the contact matrices!')
        else:
            pass
        if options.chr_col == None:
            parser.error('-h for help or insert a chromosome or a list of chromosomes for the columns of the contact matrices!')
        else:
            pass
        
        input_files = map(str, parameters['input_file'].strip('[]').split(','))
        if len(input_files) <= 1:
            parser.error('Insert at least two input files between square brackets!')
        else:
            pass
        
        chr_row_list = map(str, parameters['chr_row'].strip('[]').split(','))
        chr_col_list = map(str, parameters['chr_col'].strip('[]').split(','))
        samples_list = map(str, parameters['samples'].strip('[]').split(','))
        
        my_cmap = map(str, parameters['my_colormap'].strip('[]').split(','))
        if len(my_cmap) == 1:
            my_cmap = my_cmap[0]
        
        plot_side_by_side_map(input_files,
                              bool(parameters['tab_sep']),
                              chr_row_list,
                              chr_col_list,
                              samples_list,
                              parameters['bin_size'],  
                              parameters['data_type'],
                              parameters['species'],
                              my_cmap,
                              parameters['cutoff_type'],
                              float(parameters['cutoff']),
                              parameters['max_color'],
                              parameters['my_dpi'])
        
        
        
        

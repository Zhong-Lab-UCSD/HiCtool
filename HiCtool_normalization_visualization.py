"""
Program to:
1) Normalize the contact data (fend and enrichment).
2) Plot contact matrix, histogram and colorbar.

To use this code, fend correction values must be provided (see HiCtool_hifive.py)
"""

def save_matrix(a_matrix, output_file):
    """
    Save a square matrix in a txt file. The matrix is reshaped by rows and 
    saved in a list.
    Parameters:
        a_matrix (numpy matrix): input matrix to be saved
        output_file: output file name in txt format
    Output: 
        txt file containing the saved list          
    """
    vect = []
    n = len(a_matrix)
    for row in xrange(n):
        for col in xrange(n):
            vect.append(a_matrix[row,col])
    
    with open (output_file,'w') as fout:
        for i in xrange(n**2):
            fout.write('%s\n' %vect[i])


def load_matrix(input_file):
    """
    Load a list from a txt file and reshape it into a square matrix.
    Parameters:
        input_file: input file name in txt format.
    Returns: 
        output_matrix: array containing all the reshaped values stored 
        in the input txt file.         
    """
    import numpy as np    
    
    with open (input_file,'r') as infile:
        lines = infile.readlines()
        matrix_vect = []        
        for i in lines:
            j = i[:-1]            
            matrix_vect.append(float(j))
        
        matrix_vect_size = len(matrix_vect)
        matrix_size = int(np.sqrt(matrix_vect_size))
        output_matrix = np.reshape(matrix_vect,(matrix_size,matrix_size))
        return output_matrix
        
        
def normalize_chromosome_fend_data(a_chr, bin_size, species='hg38'):
    """
    Normalize the contact data by calculating the corrected reads count for each 
    bin. Observed data, expected fend data and normalized fend data are stored 
    into txt file.
    Parameters:
        a_chr (str): chromosome number (example for chromosome 1: '1').
        bin_size (int): bin size in bp of the contact matrix.
        species (str): hg38 or mm10.
    """
    import hifive
    import numpy as np
    
    chromosomes = {'hg38':{'1':249250621,
                   '2':243199373,
                   '3':198022430,
                   '4':191154276,
                   '5':180915260,
                   '6':171115067,
                   '7':159138663,
                   '8':146364022,
                   '9':141213431,
                   '10':135534747,
                   '11':135006516,
                   '12':133851895,
                   '13':115169878,
                   '14':107349540,
                   '15':102531392,
                   '16':89354753,
                   '17':81195210,
                   '18':77077248,
                   '19':59128983,
                   '20':61025520,
                   '21':48129895,
                   '22':51304566,
                   'X':155270560,
                   'Y':57373566},
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
    start_pos = 0
    end_pos = (chromosomes[species][a_chr]/1000000)*1000000
    
    start_part = str(float(start_pos)/float(1000000))
    end_part = str(float(end_pos)/float(1000000))
    binsize_str = str(float(bin_size)/float(1000000))
    
    # Observed data
    hic = hifive.HiC('HiC_norm_binning.hdf5')
    heatmap_enrich = hic.cis_heatmap(chrom=chromosome,
                                     start=start_pos,
                                     stop=end_pos,
                                     binsize=bin_size,
                                     arraytype='full',
                                     datatype='enrichment')
    observed = heatmap_enrich[:,:,0] # observed contact data extracted from the heatmap object
    save_matrix(observed, 'HiCtool_observed_contact_matrix_' + chromosome + '_' + binsize_str + 'mb_' + start_part + 'mb_' + end_part + 'mb.txt')
    
    # Expected raw (number of possible fend interactions). 
    # These are needed to scale the fend expected data by the mean fend pairs 
    #in each bin.
    
    hic = hifive.HiC('HiC_norm_binning.hdf5')
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
    hic = hifive.HiC('HiC_norm_binning.hdf5')
    heatmap_fend = hic.cis_heatmap(chrom=chromosome,
                                   start=start_pos,
                                   stop=end_pos,
                                   binsize=bin_size,
                                   arraytype='full',
                                   datatype='fend')
    
    # Expected fend (fend corrections)
    expected_fend = heatmap_fend[:,:,1]/scaling_factor # fend correction values
    save_matrix(expected_fend, 'HiCtool_expected_fend_contact_matrix_' + chromosome + '_' + binsize_str + 'mb_' + start_part + 'mb_' + end_part + 'mb.txt')
    
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
    
    save_matrix(normalized_fend, 'HiCtool_normalized_fend_contact_matrix_' + chromosome + '_' + binsize_str + 'mb_' + start_part + 'mb_' + end_part + 'mb.txt')


def normalize_chromosome_enrich_data(a_chr, bin_size, species='hg38'):
    """
    Calculate the enrichment data as observed/expected where the expected reads
    count is for each bin considering the distance between fends and the learned
    correction parameters. Observed, expected and enrichment contact data are saved
    to text files.
    Parameters:
        a_chr (str): chromosome number (example for chromosome 1: '1').
        bin_size (int): bin size in bp of the contact matrix.
        species (str): hg38 or mm10.
    """
    import hifive
    import numpy as np
    
    chromosomes = {'hg38':{'1':249250621,
                   '2':243199373,
                   '3':198022430,
                   '4':191154276,
                   '5':180915260,
                   '6':171115067,
                   '7':159138663,
                   '8':146364022,
                   '9':141213431,
                   '10':135534747,
                   '11':135006516,
                   '12':133851895,
                   '13':115169878,
                   '14':107349540,
                   '15':102531392,
                   '16':89354753,
                   '17':81195210,
                   '18':77077248,
                   '19':59128983,
                   '20':61025520,
                   '21':48129895,
                   '22':51304566,
                   'X':155270560,
                   'Y':57373566},
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
    start_pos = 0
    end_pos = (chromosomes[species][a_chr]/1000000)*1000000
    
    start_part = str(float(start_pos)/float(1000000))
    end_part = str(float(end_pos)/float(1000000))
    binsize_str = str(float(bin_size)/float(1000000))

    # Enrichment data
    hic = hifive.HiC('HiC_norm_binning.hdf5')
    heatmap_enrich = hic.cis_heatmap(chrom=chromosome,
                                     start=start_pos,
                                     stop=end_pos,
                                     binsize=bin_size,
                                     arraytype='full',
                                     datatype='enrichment')
    
    # Observed data
    observed = heatmap_enrich[:,:,0] # observed contact data extracted from the heatmap object
    save_matrix(observed, 'HiCtool_observed_contact_matrix_' + chromosome + '_' + binsize_str + 'mb_' + start_part + 'mb_' + end_part + 'mb.txt')
    
    # Expected enrichment data (fend corrections and distance property)
    expected_enrich = heatmap_enrich[:,:,1] # expected enrichment contact data extracted from the heatmap object
    save_matrix(expected_enrich, 'HiCtool_expected_enrich_contact_matrix_' + chromosome + '_' + binsize_str + 'mb_' + start_part + 'mb_' + end_part + 'mb.txt')
       
    # Normalized enrichment contact matrix
    n = len(expected_enrich)
    normalized_enrich = np.zeros((n,n))
    for i in xrange(n):
        for j in xrange(n):
            if expected_enrich[i][j] == 0:
                normalized_enrich[i][j] = 0
            else:
                normalized_enrich[i][j] = float(observed[i][j])/float(expected_enrich[i][j])
    
    save_matrix(n, normalized_enrich, 'HiCtool_normalized_enrich_contact_matrix_' + chromosome + '_' + binsize_str + 'mb_' + start_part + 'mb_' + end_part + 'mb.txt')
    

def make_cmap(colors):
    '''
    Take a list of tuples which contain RGB values and return a
    cmap with equally spaced colors.
    Arrange your tuples so that the first color is the lowest value for the
    colorbar and the last is the highest.
    '''
    import matplotlib as mpl
    import numpy as np

    bit_rgb = np.linspace(0,1,256)
    position = np.linspace(0,1,len(colors))
    for i in range(len(colors)):
        colors[i] = (bit_rgb[colors[i][0]],
                     bit_rgb[colors[i][1]],
                     bit_rgb[colors[i][2]])
    cdict = {'red':[], 'green':[], 'blue':[]}
    for pos, color in zip(position, colors):
        cdict['red'].append((pos, color[0], color[0]))
        cdict['green'].append((pos, color[1], color[1]))
        cdict['blue'].append((pos, color[2], color[2]))
    cmap = mpl.colors.LinearSegmentedColormap('my_colormap',cdict,256)
    return cmap


def plot_normalized_chromosome_fend_data(a_chr, start_coord, end_coord, bin_size,
                                         species='hg38',
                                         plot_histogram=True,
                                         plot_colorbar=True):
    """
    Plot normalized fend contact map, colorbar and histogram.
    Parameters:
        a_chr (str): chromosome number (example for chromosome 1: '1').
        start_coord (int): start coordinate for the plot in bp.
        end_coord (int): end coordinate for the plot in bp.
        bin_size (int): bin size in bp of the contact matrix.
        species (str): hg38 or mm10.
        plot_histogram (bool): if true, plot the histogram.
        plot_colorbar (bool): if true, plot the colorbar.
    """                                       
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    from PIL import Image
    
    chromosomes = {'hg38':{'1':249250621,
                   '2':243199373,
                   '3':198022430,
                   '4':191154276,
                   '5':180915260,
                   '6':171115067,
                   '7':159138663,
                   '8':146364022,
                   '9':141213431,
                   '10':135534747,
                   '11':135006516,
                   '12':133851895,
                   '13':115169878,
                   '14':107349540,
                   '15':102531392,
                   '16':89354753,
                   '17':81195210,
                   '18':77077248,
                   '19':59128983,
                   '20':61025520,
                   '21':48129895,
                   '22':51304566,
                   'X':155270560,
                   'Y':57373566},
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

    start_pos = 0
    end_pos = (chromosomes[species][a_chr]/1000000)*1000000
    
    start_part = str(float(start_pos)/float(1000000))
    end_part = str(float(end_pos)/float(1000000))
    binsize_str = str(float(bin_size)/float(1000000))
    
    # Plotting of the fend normalized data
    matrix_data_full = load_matrix('HiCtool_normalized_fend_contact_matrix_chr' + a_chr + '_' + binsize_str + 'mb_' + start_part + 'mb_' + end_part + 'mb.txt')
    
    # Selecting a part
    start_bin = start_coord/bin_size
    end_bin = end_coord/bin_size
    
    start_part = str(float(start_coord)/float(1000000))
    end_part = str(float(end_coord)/float(1000000))
    matrix_data_full = matrix_data_full[start_bin:end_bin+1,start_bin:end_bin+1]
    ##################
    
    n = len(matrix_data_full)
    output_vect = np.reshape(matrix_data_full,n*n,1)
    non_zero = np.nonzero(output_vect)
    perc = np.percentile(output_vect[non_zero[0]],98)
    for i in xrange(len(matrix_data_full)):
        for j in xrange(len(matrix_data_full)):
            if matrix_data_full[i][j] > perc:
                matrix_data_full[i][j] = perc
    
    # Heatmap
    max_value = np.max(matrix_data_full)
    min_value = np.min(matrix_data_full)
    norm_matrix_data_full = 255-((matrix_data_full-min_value)/(max_value-min_value))*255 # normalization to have the data between white and red
    
    img = Image.new('RGB',(n,n))
    newData = []
    for i in xrange(n):
        for j in xrange(n):
            value = int(norm_matrix_data_full[i][j])
            newData.append((255,value,value))
    
    img.putdata(newData)
    img.save('HiCtool_normalized_fend_contact_matrix_chr' + a_chr + '_' + binsize_str + 'mb_' + start_part + 'mb_' + end_part + 'mb.png','PNG')
    
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
        plt.title('Contact frequency histogram')
        plt.xlabel('Number of contacts')
        plt.ylabel('Frequency')
        plt.savefig('HiCtool_normalized_fend_histogram_chr' + a_chr + '_' + binsize_str + 'mb_' + start_part + 'mb_' + end_part + 'mb.png')

    # Plot of the colorbar
    if plot_colorbar:
        
        bar_min = np.min(matrix_data_full)
        bar_max = np.max(matrix_data_full)
        colors = []
        for i in xrange(256):
            color = (255,255-i,255-i)
            colors.append(color)
        
        plt.close("all")
        fig = plt.figure(figsize=(1.5, 7))
        ax = fig.add_axes([0.3, 0.08, 0.4, 0.9])
        cmap = make_cmap(colors)
        norm = matplotlib.colors.Normalize(vmin=bar_min, vmax=bar_max)
        matplotlib.colorbar.ColorbarBase(ax, 
                                         cmap=cmap,
                                         norm=norm,
                                         orientation='vertical')
        plt.savefig('HiCtool_normalized_fend_colorbar_chr' + a_chr + '_' + binsize_str + 'mb_' + start_part + 'mb_' + end_part + 'mb.png')
            

def plot_normalized_chromosome_enrich_data(a_chr, start_coord, end_coord, bin_size,
                                           species='hg38',
                                           plot_histogram=True,
                                           plot_colorbar=True):
    """
    Plot the log2 of the normalized enrichment contact map, colorbar and histogram.
    Parameters:
        a_chr (str): chromosome number (example for chromosome 1: '1').
        start_coord (int): start coordinate for the plot in bp.
        end_coord (int): end coordinate for the plot in bp.
        bin_size (int): bin size in bp of the contact matrix.
        species (str): hg38 or mm10.
        plot_histogram (bool): if true, plot the histogram.
        plot_colorbar (bool): if true, plot the colorbar.
    """                                       
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    import math
    from PIL import Image
    
    chromosomes = {'hg38':{'1':249250621,
                   '2':243199373,
                   '3':198022430,
                   '4':191154276,
                   '5':180915260,
                   '6':171115067,
                   '7':159138663,
                   '8':146364022,
                   '9':141213431,
                   '10':135534747,
                   '11':135006516,
                   '12':133851895,
                   '13':115169878,
                   '14':107349540,
                   '15':102531392,
                   '16':89354753,
                   '17':81195210,
                   '18':77077248,
                   '19':59128983,
                   '20':61025520,
                   '21':48129895,
                   '22':51304566,
                   'X':155270560,
                   'Y':57373566},
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

    start_pos = 0
    end_pos = (chromosomes[species][a_chr]/1000000)*1000000
    
    start_part = str(float(start_pos)/float(1000000))
    end_part = str(float(end_pos)/float(1000000))
    binsize_str = str(float(bin_size)/float(1000000))
    
    # Plotting the enrichment contact data
    matrix_data_full = load_matrix('HiCtool_normalized_enrich_contact_matrix_chr' + a_chr + '_' + binsize_str + 'mb_' + start_part + 'mb_' + end_part + 'mb.txt')
    
    # Selecting a part
    start_bin = start_coord/bin_size
    end_bin = end_coord/bin_size
    
    start_part = str(float(start_coord)/float(1000000))
    end_part = str(float(end_coord)/float(1000000))
    matrix_data_full = matrix_data_full[start_bin:end_bin+1,start_bin:end_bin+1]
    ##################
    
    n = len(matrix_data_full)
    output_vect = np.reshape(matrix_data_full,n*n,1)
    non_zero = np.nonzero(output_vect)
    non_zero_values = output_vect[non_zero[0]]
    positive_logs = []
    negative_logs = []
    negative_logs_abs = []
    zero_logs = []
    for i in non_zero_values:
        log_value = math.log(i,2)
        if log_value == 0:
            zero_logs.append(log_value)
        if log_value > 0:
            positive_logs.append(log_value)
        if log_value < 0:
            negative_logs.append(log_value)
            negative_logs_abs.append(abs(log_value))
    
    max_value = np.percentile(positive_logs,99)
    min_value = np.percentile(negative_logs_abs,99)
    min_value = -min_value
    
    # Heatmap
    img = Image.new('RGB',(n,n))
    newData = []
    for i in xrange(n):
        for j in xrange(n):
            value = matrix_data_full[i][j]
            if value==0:
                newData.append((100,100,100))
                continue
            log_value = math.log(value,2)
            if log_value < 0:
                if log_value < min_value: log_value = min_value
                color_value = int(((log_value-min_value)/(abs(min_value)))*255)
                newData.append((color_value,color_value,255))
            if log_value >= 0:
                if log_value > max_value: log_value = max_value
                color_value = int((log_value/max_value)*255)
                newData.append((255,255-color_value,255-color_value))
    
    img.putdata(newData)
    img.save('HiCtool_normalized_enrich_contact_matrix_chr' + a_chr + '_' + binsize_str + 'mb_' + start_part + 'mb_' + end_part + 'mb.png','PNG')
    
    # Plot the histogram
    if plot_histogram:

        for i in xrange(len(positive_logs)):
            if positive_logs[i] > max_value:
                positive_logs[i] = max_value
        for i in xrange(len(negative_logs)):
            if negative_logs[i] < min_value:
                negative_logs[i] = min_value
        
        logs_values = positive_logs + negative_logs + zero_logs
        s = set([x for x in logs_values if logs_values.count(x) > 1])
        for i in s:
            c = logs_values.count(i)
            for j in xrange(c/2):
                logs_values.remove(i)
        
        plt.close("all")
        histogram_bins = int(pow(len(logs_values),0.3))
        plt.hist(logs_values, bins=histogram_bins)
        plt.title('Enrichment histogram')
        plt.xlabel('log2(O/E)')
        plt.ylabel('Frequency')
        plt.savefig('HiCtool_normalized_enrich_histogram_chr' + a_chr + '_' + binsize_str + 'mb_' + start_part + 'mb_' + end_part + 'mb.png')

    # Plot the colorbar
    if plot_colorbar:
        # Positive logs
        bar_min = 0
        bar_max = max_value
        colors = []
        for i in xrange(256):
            color = (255,255-i,255-i)
            colors.append(color)
        
        plt.close("all")
        fig = plt.figure(figsize=(1.5, 3.5))
        ax = fig.add_axes([0.3, 0.08, 0.4, 0.8])
        cmap = make_cmap(colors)
        norm = matplotlib.colors.Normalize(vmin=bar_min, vmax=bar_max)
        matplotlib.colorbar.ColorbarBase(ax, 
                                         cmap=cmap,
                                         norm=norm,
                                         orientation='vertical')
        plt.savefig('HiCtool_normalized_enrich_colorbar_red_chr' + a_chr + '_' + binsize_str + 'mb_' + start_part + 'mb_' + end_part + 'mb.png')
        
        # Negative logs
        bar_min = min_value
        bar_max = 0
        colors = []
        for i in xrange(256):
            color = (i,i,255)
            colors.append(color)
        
        plt.close("all")
        fig = plt.figure(figsize=(1.5, 3.5))
        ax = fig.add_axes([0.3, 0.08, 0.4, 0.8])
        cmap = make_cmap(colors)
        norm = matplotlib.colors.Normalize(vmin=bar_min, vmax=bar_max)
        matplotlib.colorbar.ColorbarBase(ax, 
                                         cmap=cmap,
                                         norm=norm,
                                         orientation='vertical')
        plt.savefig('HiCtool_normalized_enrich_colorbar_blue_chr' + a_chr + '_' + binsize_str + 'mb_' + start_part + 'mb_' + end_part + 'mb.png')

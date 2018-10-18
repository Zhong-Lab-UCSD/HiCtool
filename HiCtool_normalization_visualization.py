"""
Program to:
1) Normalize the contact data (fend and enrichment ("observed / expected")).
2) Plot contact matrix and histogram of contact distribution.

To use this code, fend correction values must be provided (see HiCtool_hifive.py)
"""

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


def load_matrix(input_file):
    """
    Load a formatted contact matrix from a txt file and parse it.
    Parameters:
        input_file: input file name in txt format (generated by the function 
        "save_matrix")
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
  
    k = len(matrix_vect)
    matrix_size = int((-1+np.sqrt(1+8*k))/2)
    
    iu = np.triu_indices(matrix_size)
    output_matrix_1 = np.zeros((matrix_size,matrix_size)) # upper triangular plus the diagonal
    output_matrix_1[iu] = matrix_vect
    
    diag_matrix = np.diag(np.diag(output_matrix_1)) # diagonal
    output_matrix_2 = np.transpose(output_matrix_1) # lower triangular plus the diagonal
    output_matrix = output_matrix_1 + output_matrix_2 - diag_matrix
    
    return output_matrix
    print "Done!"

def normalize_chromosome_fend_data(a_chr, 
                                   bin_size, 
                                   input_file='HiC_norm_binning.hdf5',
                                   species='hg38', 
                                   chr_size=0,
                                   save_obs=True, 
                                   save_expect=False):
    """
    Normalize the contact data by calculating the corrected reads count for each 
    bin. Observed data, expected fend data and normalized fend data are saved 
    into txt file.
    Parameters:
        a_chr (str): chromosome number (example for chromosome 1: '1').
        bin_size (int): bin size in bp of the contact matrix.
        input_file (str): object containing learned correction parameters in .hdf5 format obtained with
        HiCtool_hifive.py (default: 'HiC_norm_binning.hdf5')
        species (str): 'hg38' or 'mm10' or any other species label in string format.
        chr_size (int): chromosome size of your custom species if you did not use 'hg38' or 'mm10'.
        save_obs (bool): if true, save the observed contact data.
        save_expect (bool): if true, save the expected contact data.
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
    
    print "Normalizing fend data..."
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
        end_pos = (chr_size/bin_size)*bin_size
        
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
    print "Done!"

def plot_chromosome_fend_data(contact_matrix,
                              a_chr, 
                              start_coord, 
                              end_coord, 
                              bin_size,
                              species='hg38',
                              chr_size=0,
                              my_colormap=['white', 'red'],
                              cutoff=95,
                              max_color='#460000',
                              my_dpi=1000,
                              plot_histogram=False):
    """
    Plot a contact map and histogram of the contact distribution generated with the function "normalize_chromosome_fend_data".
    Parameters:
        contact_matrix (str): contact matrix generated with the function "normalize_chromosome_fend_data".
        a_chr (str): chromosome number (example for chromosome 1: '1').
        start_coord (int): start coordinate for the plot in bp.
        end_coord (int): end coordinate for the plot in bp.
        bin_size (int): bin size in bp of the contact matrix.
        species (str): 'hg38' or 'mm10' or any other species label in string format.
        chr_size (int): chromosome size of your custom species if you did not use 'hg38' or 'mm10'.
        my_colormap (str | list): colormap to be used to plot the data. 1) Use a string if you choose among any colorbar here 
        https://matplotlib.org/examples/color/colormaps_reference.html 2) Use a list of strings with colors if you want
        a custom colorbar. Example: ['white', 'red', 'black']. Colors can be specified also in this format: '#000000'.
        cutoff (int): percentile to set a maximum cutoff on the number of contacts for the colorbar.
        max_color (str): to set the color of the bins with contact counts over "cutoff".
        my_dpi (int): resolution of the contact map in dpi.
        plot_histogram (bool): if true, plot and save to file the histogram of the contact distribution.
    """                                           
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
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
        output_filename = 'HiCtool_' + chromosome + '_' + bin_size_str + 'mb_normalized_fend'
    elif bin_size < 1000000:
        bin_size_str = str(bin_size/1000)
        output_filename = 'HiCtool_' + chromosome + '_' + bin_size_str + 'kb_normalized_fend'
    
    if species == 'hg38' or species == 'mm10':
        end_pos = (chromosomes[species][a_chr]/bin_size)*bin_size
    else:
        end_pos = (chr_size/bin_size)*bin_size
    
    # Plotting of the fend normalized data
    matrix_data_full = load_matrix(contact_matrix)
    
    # Selecting a part
    start_bin = start_coord/bin_size
    end_bin = end_coord/bin_size
    
    if start_coord >= end_coord:
        print "WARNING! Start coordinate should be lower than end coordinate"
        return
        
    if start_bin >= end_bin:
        print "WARNING! Start coordinate should be much lower than the end coordinate given the bin size"
        return
    
    if end_coord > end_pos:
        if species == 'hg38' or species == 'mm10':
            print "WARNING! End coordinate is larger than chromosome size " + str((chromosomes[species][a_chr]/bin_size)*bin_size) + " bp"
            return
        else:
            print "WARNING! End coordinate is larger than chromosome size " + str((chr_size/bin_size)*bin_size) + " bp"
            return
    
    matrix_data_full = matrix_data_full[start_bin:end_bin+1,start_bin:end_bin+1] 
    
    n = len(matrix_data_full)
    output_vect = np.reshape(matrix_data_full,n*n,1)
    non_zero = np.nonzero(output_vect)
    if non_zero[0].size == 0:
        print "WARNING! The portion of chromosome you selected contains no data."
        return
    
    # Heatmap plotting
    print "Plotting " + contact_matrix + "..."
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
    
    if cutoff == 0:
        plt.close("all")
        plt.imshow(matrix_data_full, cmap=my_cmap, interpolation='nearest')
        plt.title('Normalized contact map (' + my_bin_size + ')', fontsize=14)
        plt.xlabel(chromosome + ' coordinate (bp)', fontsize=12)
        plt.ylabel(chromosome + ' coordinate (bp)', fontsize=12)
        cbar = plt.colorbar()
        cbar.ax.set_ylabel('Normalized contact counts', rotation=270, labelpad=20)
        ticks = (np.arange(0, n, n/4) * bin_size) + start_coord
        format_ticks = [format_e(i) for i in ticks.tolist()]
        plt.xticks(np.arange(0, n, n/4), format_ticks)
        plt.yticks(np.arange(0, n, n/4), format_ticks)
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.savefig(output_filename + '.pdf', format = 'pdf', dpi=my_dpi)
        
    else:
        perc = np.percentile(output_vect[non_zero[0]],cutoff)
        
        plt.close("all")
        plt.imshow(matrix_data_full, cmap=my_cmap, interpolation='nearest', vmax=perc)
        plt.title('Normalized contact map (' + my_bin_size + ')', fontsize=14)
        plt.xlabel(chromosome + ' coordinate (bp)', fontsize=12)
        plt.ylabel(chromosome + ' coordinate (bp)', fontsize=12)
        cbar = plt.colorbar(extend='max')
        cbar.cmap.set_over(max_color)
        cbar.ax.set_ylabel('Normalized contact counts', rotation=270, labelpad=20)
        ticks = (np.arange(0, n, n/4) * bin_size) + start_coord
        format_ticks = [format_e(i) for i in ticks.tolist()]
        plt.xticks(np.arange(0, n, n/4), format_ticks)
        plt.yticks(np.arange(0, n, n/4), format_ticks)
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.savefig(output_filename + '.pdf', format = 'pdf', dpi=my_dpi)
    
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
        plt.title('Normalized contact counts distribution', fontsize=18)
        plt.xlabel('Normalized contact counts', fontsize=16)
        plt.ylabel('Number of bins', fontsize=16)
        plt.xticks(fontsize=16)
        plt.yticks(fontsize=16)
        plt.savefig(output_filename + '_histogram.pdf', format = 'pdf')
    print "Done!"

def normalize_chromosome_enrich_data(a_chr, 
                                     bin_size,
                                     input_file='HiC_norm_binning.hdf5',
                                     species='hg38', 
                                     chr_size=0,
                                     save_obs=True, 
                                     save_expect=False):
    """
    Calculate the enrichment data as "observed/expected" where the expected reads
    count is for each bin considering the linear distance between read pairs and the learned
    correction parameters. Observed, expected and enrichment contact data are saved
    to txt files.
    Parameters:
        a_chr (str): chromosome number (example for chromosome 1: '1').
        bin_size (int): bin size in bp of the contact matrix.
        input_file (str): object containing learned correction parameters in .hdf5 format obtained with
        HiCtool_hifive.py (default: 'HiC_norm_binning.hdf5').
        species (str): 'hg38' or 'mm10' or any other species label in string format.
        chr_size (int): chromosome size of your custom species if you did not use 'hg38' or 'mm10'.
        save_obs (bool): if true, save the observed contact data.
        save_expect (bool): if true, save the expected contact data.
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
    
    print "Normalizing enrichment data..."
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
        end_pos = (chr_size/bin_size)*bin_size

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
        save_matrix(observed, output_filename + 'observed.txt')            
            
    # Expected enrichment data (fend corrections and distance property)
    expected_enrich = heatmap_enrich[:,:,1] # expected enrichment contact data extracted from the heatmap object
    if save_expect == True:  
        save_matrix(expected_enrich, output_filename + 'expected_enrich.txt')
    
    # Normalized enrichment contact matrix
    n = len(expected_enrich)
    normalized_enrich = np.zeros((n,n))
    for i in xrange(n):
        for j in xrange(n):
            if expected_enrich[i][j] == 0:
                normalized_enrich[i][j] = -1
            else:
                normalized_enrich[i][j] = float(observed[i][j])/float(expected_enrich[i][j])
    
    save_matrix(normalized_enrich, output_filename + 'normalized_enrich.txt')
    print "Done!"

def plot_chromosome_enrich_data(contact_matrix,
                                a_chr, 
                                start_coord, 
                                end_coord, 
                                bin_size,
                                species='hg38',
                                chr_size=0,
                                my_dpi=1000,
                                plot_histogram=False):
    """
    Plot the log2 of the "observed / expected" contact map and histogram of the enrichment values distribution 
    generated with the function "normalize_chromosome_enrich_data".
    Parameters:
        contact_matrix (str): "observed / expected" contact matrix generated with the function "normalize_chromosome_enrich_data".
        a_chr (str): chromosome number (example for chromosome 1: '1').
        start_coord (int): start coordinate for the plot in bp.
        end_coord (int): end coordinate for the plot in bp.
        bin_size (int): bin size in bp of the contact matrix.
        species (str): 'hg38' or 'mm10' or any other species label in string format.
        chr_size (int): chromosome size of your custom species if you did not use 'hg38' or 'mm10'.
        my_dpi (int): resolution of the contact map in dpi.
        plot_histogram (bool): if true, plot the histogram.
    """                                       
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    import math
    from numpy import ma
    from matplotlib import cbook
    from matplotlib.colors import Normalize
    
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
        output_filename = 'HiCtool_' + chromosome + '_' + bin_size_str + 'mb_normalized_enrich'
    elif bin_size < 1000000:
        bin_size_str = str(bin_size/1000)
        output_filename = 'HiCtool_' + chromosome + '_' + bin_size_str + 'kb_normalized_enrich'
    
    if species == 'hg38' or species == 'mm10':
        end_pos = (chromosomes[species][a_chr]/bin_size)*bin_size
    else:
        end_pos = (chr_size/bin_size)*bin_size
    
    # Plotting the enrichment contact data
    matrix_data_full = load_matrix(contact_matrix)
    
    # Selecting a part
    start_bin = start_coord/bin_size
    end_bin = end_coord/bin_size
    
    if start_coord >= end_coord:
        print "WARNING! Start coordinate should be lower than end coordinate"
        return
        
    if start_bin >= end_bin:
        print "WARNING! Start coordinate should be much lower than the end coordinate given the bin size"
        return
    
    if end_coord > end_pos:
        if species == 'hg38' or species == 'mm10':
            print "WARNING! End coordinate is larger than chromosome size " + str((chromosomes[species][a_chr]/bin_size)*bin_size) + " bp"
            return
        else:
            print "WARNING! End coordinate is larger than chromosome size " + str((chr_size/bin_size)*bin_size) + " bp"
            return
    
    matrix_data_full = matrix_data_full[start_bin:end_bin+1,start_bin:end_bin+1]
    n = len(matrix_data_full)
    
    # Heatmap plotting
    print "Plotting " + contact_matrix + "..."
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
                result.fill(0) # Or should it be all masked? Or 0.5?
            elif vmin > vmax:
                raise ValueError("maxvalue must be bigger than minvalue")
            else:
                vmin = float(vmin)
                vmax = float(vmax)
                if clip:
                    mask = ma.getmask(result)
                    result = ma.array(np.clip(result.filled(vmax), vmin, vmax),
                                      mask=mask)
    
                # ma division is very slow; we can take a shortcut
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
    norm = MidPointNorm(midpoint=0)
    plt.imshow(matrix_data_full, cmap='seismic', interpolation='nearest', vmax=threshold_max, vmin=threshold_min, norm=norm)
    plt.title('"Observed/expected" contact map (' + my_bin_size + ')', fontsize=14)
    plt.xlabel(chromosome + ' coordinate (bp)', fontsize=12)
    plt.ylabel(chromosome + ' coordinate (bp)', fontsize=12)
    cbar = plt.colorbar()
    cbar.cmap.set_over('black') # loci where expected enrich values are zero (log not existing)
    cbar.cmap.set_under('gray') # loci where observed values are zero (log equal to minus infinity)
    cbar.ax.set_ylabel('log2(O/E) contact counts', rotation=270, labelpad = 20)
    ticks = (np.arange(0, n, n/4) * bin_size) + start_coord
    format_ticks = [format_e(i) for i in ticks.tolist()]
    plt.xticks(np.arange(0, n, n/4), format_ticks)
    plt.yticks(np.arange(0, n, n/4), format_ticks)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.savefig(output_filename + '.pdf', format = 'pdf', dpi = my_dpi)
    
    # Plot the histogram
    if plot_histogram:        
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
        plt.title('"Observed/expected" contact counts distribution', fontsize=18)
        plt.xlabel('log2(O/E) contact counts', fontsize=16)
        plt.ylabel('Number of bins', fontsize=16)
        plt.xticks(fontsize=16)
        plt.yticks(fontsize=16)
        plt.savefig(output_filename + '_histogram.pdf', format = 'pdf')
    print "Done!"
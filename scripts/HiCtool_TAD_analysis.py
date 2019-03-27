"""
Program to:
1) Calculate and plot the observed DI and true DI (Hidden Markov Model).
2) Calculate the topological domains coordinates.

To use this code, normalized contact data must be provided.
"""

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

def save_list(a_list, output_file):
    """
    Save a list in a .txt file.
    Parameters:
        a_list: name of the list to save.
        output_file: output file name in txt format.
    Output: 
        .txt file containing the saved list.       
    """
    with open (output_file,'w') as fout:
        n = len(a_list)
        for i in xrange(n):
            fout.write('%s\n' %a_list[i])

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
        "save_matrix").
    Returns: 
        output_matrix: array containing the parsed values stored 
        in the input txt file to build a contact matrix.      
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
                    

def load_matrix_tab(input_filename):
    """
    Load a contact matrix saved in a tab separated format using the function
    "save_matrix_tab".
    Parameters:
    input_filename (str): input contact matrix to be loaded
    Returns:
    output_matrix: array containing the parsed values stored 
    in the input tab separated txt file to build a contact matrix
    """
    import numpy as np
    
    print "Loading " + input_filename + "..."
    with open (input_filename, 'r') as infile:
        lines = infile.readlines()
        temp = []
        for line in lines:
            row = [float(i) for i in line.strip().split('\t')]
            temp.append(row)
            
        output_matrix = np.array(temp)
    print "Done!"
    return output_matrix


    
def load_DI_values(input_file):
    """
    Load a DI txt file generated with "calculate_chromosome_DI".
    Parameters:
        input_file: input file name in txt format.
    Returns: 
        List of the DI values.        
    """
    import numpy as np
    
    fp = open(input_file,'r+')
    lines = fp.read().split('\n')
    lines = lines[:-1]
    di_values = (np.nan_to_num(np.array(map(float, lines)))).tolist()
    return di_values


def extract_single_map(input_global_matrix,
                       tab_sep,
                       chr_row,
                       chr_col,
                       species='hg38',
                       bin_size=1000000,
                       data_type='observed',
                       custom_species_sizes={},
                       sexual_chromosomes=[],
                       save_output=True,
                       save_tab=False):
    """
    Extract a single contact matrix for a pair of chromosomes from the full matrix.
    Parameters:
        input_global_matrix (object | str): full contact matrix. This can be passed either as
        an object of the workspace or a string of the filename saved to file.
        tab_sep (bool): if "input_global_matrix" is passed with a filename, then this boolean 
        tells if the global matrix was saved in tab separated format (True) or not (False).
        chr_row (str): chromosome in the rows of the output contact matrix.
        chr_col (str): chromosome in the columns of the output contact matrix. If chr_col is 
        equal to chr_row then the intrachromosomal map is extracted.
        species (str): 'hg38' or 'mm10' or any other species label in string format.
        bin_size (int): bin size in bp of the contact matrix.
        data_type (str): which kind of data type you are extracting. "observed" or "normalized".
        custom_species_sizes (dict): dictionary containing the sizes of the chromosomes
        of your custom species. The keys of the dictionary are chromosomes in string
        format (example for chromosome 1: '1'), the values are chromosome lengths as int.
        sexual_chromosomes (list): list of the sexual chromosomes (if present) in your
        custom species (example for chromosome X: 'X').
        save_output (bool): if true, save the contact matrix in formatted txt file.
        save_tab (bool): if true, save the contact matrix in tab separated format.
    Return:
        Contact matrix in numpy array format.
    """            
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
        #dim_row = str(d_chr_dim[chr_row])
        #dim_col = str(d_chr_dim[chr_col])
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


def calculate_chromosome_DI(input_contact_matrix, 
                            a_chr,
                            isGlobal,
                            tab_sep=False,
                            species='hg38',
                            custom_species_sizes={},
                            sexual_chromosomes=[],
                            save_file=True):
    """
    Function to calculate the DI values for a chromosome of and save them 
    in a txt file.
    Parameters:
        input_contact_matrix (str | obj): normalized intra-chromosomal contact matrix at a bin size of 40kb passed as a filename (str)
        or an object. Either a single contact matrix or a global contact matrix can be passed (see following parameters).
        a_chr (str): chromosome number (example for chromosome 1: '1').
        isGlobal (bool): set True if your input matrix is a global matrix.
        tab_sep (bool): set True if your input matrix is in a tab separated format. If the matrix is passed as an
        object, this parameter is not taken into consideration.
        species (str): 'hg38' or 'mm10' or any other species label in string format.
        custom_species_sizes (dict): dictionary containing the sizes of the chromosomes
        of your custom species. The keys of the dictionary are chromosomes in string
        format (example for chromosome 1: '1'), the values are chromosome lengths as int.
        sexual_chromosomes (list): list of the sexual chromosomes (if present) in your
        custom species (example for chromosome X: 'X').
        save_file (bool): if True, saves DI values to txt file.
    Returns:
        List with the DI values.
    """
    
    import copy
    
    if isGlobal == False:
        if isinstance(input_contact_matrix, str):
            if tab_sep == False:
                contact_matrix = load_matrix(input_contact_matrix)
            else:
                contact_matrix = load_matrix_tab(input_contact_matrix)
        else:
            contact_matrix = copy.deepcopy(input_contact_matrix)
    else:
        contact_matrix = extract_single_map(input_global_matrix=input_contact_matrix,
                                            tab_sep=tab_sep,
                                            chr_row=a_chr,
                                            chr_col=a_chr,
                                            species=species,
                                            bin_size=40000,
                                            data_type='',
                                            custom_species_sizes=custom_species_sizes,
                                            sexual_chromosomes=sexual_chromosomes,
                                            save_output=False,
                                            save_tab=False)
        
    print "Calculating DI values..."
    n = contact_matrix.shape[0]

    # Calculation of the DI
    DI = [] # list of the DI for each bin
    len_var = 2000000/40000 # range of upstream or downstream bins to calculate DI

    for locus in xrange(n): # 'locus' refers to a bin
        if locus < len_var:
            A = sum(contact_matrix[locus][:locus])
            B = sum(contact_matrix[locus][locus+1:locus+len_var+1])
        elif locus >= n-len_var:
            A = sum(contact_matrix[locus][locus-len_var:locus])
            B = sum(contact_matrix[locus][locus+1:])
        else:
            A = sum(contact_matrix[locus][locus-len_var:locus])
            B = sum(contact_matrix[locus][locus+1:locus+len_var+1])
    
        E = (A+B)/2 # expected number of reads
        
        if A==0 and B==0:
            di = 0
            DI.append(di)
        else:
            try:
                di = ((B-A)/(abs(B-A)))*((((A-E)**2)/E)+(((B-E)**2)/E))
            except ZeroDivisionError:
                di = 0
            DI.append(di)
    
    if save_file == True:
        save_list(DI,'HiCtool_chr' + a_chr + '_DI.txt')
    
    print "Done!"
    return DI
    
    
def calculate_chromosome_true_DI(input_file_DI,
                                 a_chr,
                                 save_file=True):
    """
    Function to calculate the true DI values for a chromosome and save 
    them in a .txt file.
    Input:
        input_file_DI (str | obj): txt file of the DI values generated with the function "calculate_chromosome_DI" or
        object with the DI values returned by "calculate_chromosome_DI".
        a_chr (str): chromosome number (example for chromosome 1: '1').
        save_file (bool): if True, saves DI values to txt file.
    Returns:
        List with the true DI values (HMM states).
    """
    import numpy as np
    import hmmlearn.hmm as hmm
    
    print "Calculating true DI values..."
    if isinstance(input_file_DI,str):
        A = load_DI_values(input_file_DI)
    else:
        A = input_file_DI   
    
    # Guessed Transition Matrix
    TRANS_GUESS = np.array([[0.4, 0.3, 0.3],
                             [0.3, 0.4, 0.3],
                             [0.3, 0.3, 0.4]])
    
    # Guessed Emission Matrix   
    EMISS_GUESS = np.array([[0.4, 0.3, 0.3],
                             [0.3, 0.4, 0.3],
                             [0.3, 0.3, 0.4]])
    
    # Observed emissions                  
    emissions = []
    zero_threshold = 0.4;
    
    for i in range(0,len(A)):
        if A[i] >= zero_threshold:
            emissions.append(1)
        elif A[i] <= -zero_threshold:
            emissions.append(2)
        else:
            emissions.append(0)
    
    # Hidden Markov Model with discrete emissions
    model = hmm.MultinomialHMM(n_components=3, init_params="")
    model.transmat_ = TRANS_GUESS
    model.emissionprob_ = EMISS_GUESS
    
    input_observations = np.array([emissions]).T
    model.fit(input_observations) # estimate model parameters
    
    # Find most likely state sequence corresponding to input_onservations using the Viterbi Algorithm
    logprob, likelystates_array = model.decode(input_observations, algorithm="viterbi")
    
    likelystates = likelystates_array.tolist()
    
    if save_file == True:
        save_list(likelystates, "HiCtool_chr" + a_chr + "_hmm_states.txt")
    
    print "Done!"
    return likelystates


def load_hmm_states(input_file):
    """
    Load a hmm txt file generated with "calculate_chromosome_true_DI".
    Parameters:
        input_file: input file name in txt format.
    Returns: 
        List of the hmm states.    
    """
    fp = open(input_file,'r+')
    lines = fp.read().split('\n')
    lines = lines[:-1]
    likelystates = map(int,lines)
    return likelystates


def plot_chromosome_DI(input_file_DI, 
                       a_chr,
                       full_chromosome,
                       start_pos=0, 
                       end_pos=0,
                       input_file_hmm='',
                       species='hg38',
                       custom_species_sizes={},
                       plot_legend=True,
                       plot_grid=True):
    """
    Function to plot the DI and true DI values for a chromosome.
    Parameters:
        input_file_DI (str | obj): txt file of the DI values generated with the function "calculate_chromosome_DI" or
        object with the DI values returned by "calculate_chromosome_DI".
        a_chr (str): chromosome number (example for chromosome 1: '1').
        full_chromosome (bool): if True, plot the full chromosome "a_chr". In this case "start_pos" and "end_pos" parameters are not considered.
        start_pos (int): start coordinate for the plot in bp.
        end_pos (int): end coordinate for the plot in bp.
        input_file_hmm (str | obj): txt file of the true DI values generated with the function "calculate_chromosome_true_DI" or
        object with the true DI values returned by "calculate_chromosome_true_DI.
        species (str): species name (hg38, mm10, etc.).
        custom_species_sizes (dict): dictionary containing the sizes of the chromosomes
        of your custom species. The keys of the dictionary are chromosomes in string
        format (example for chromosome 1: '1'), the values are chromosome lengths as int.
        plot_legend (bool): if True, plot the legend.
        plot_grid (bool): if True, plot the grid.
    """    
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib
    matplotlib.use('Agg')
    
    bin_size = 40000
    
    if full_chromosome == True:
        if species in chromosomes.keys():
            start_pos = 0
            end_pos = int(round((chromosomes[species][a_chr])/bin_size))*bin_size
            start_index = 0
            end_index = int(round((chromosomes[species][a_chr])/bin_size)) + 1
        else:
            if len(custom_species_sizes) == 0:
                print "ERROR: insert custom_species_sizes parameter"
            start_pos = 0
            end_pos = int(round((custom_species_sizes[a_chr])/bin_size))*bin_size
            start_index = 0
            end_index = int(round((custom_species_sizes[a_chr])/bin_size)) + 1
    else:
        if end_pos == 0:
            print "ERROR: insert start and end coordinates"
            return
        start_index = int(round(start_pos/bin_size))
        end_index = int(round((end_pos)/bin_size))
    
        if species in chromosomes.keys():
            if end_pos > int(round((chromosomes[species][a_chr])/bin_size))*bin_size and end_pos <= chromosomes[species][a_chr]:
                end_pos = int(round((chromosomes[species][a_chr])/bin_size))*bin_size
            elif end_pos > chromosomes[species][a_chr]:
                print("ERROR: end coordinate exceeds chromosome dimension")
                return
        else:
            if end_pos > int(round((custom_species_sizes[a_chr])/bin_size))*bin_size and end_pos <= custom_species_sizes[a_chr]:
                end_pos = int(round((custom_species_sizes[a_chr])/bin_size))*bin_size
            elif end_pos > custom_species_sizes[a_chr]:
                print("ERROR: end coordinate exceeds chromosome dimension")
                return
    
    if isinstance(input_file_DI,str):
        DI = load_DI_values(input_file_DI)
    else:
        DI = input_file_DI   
        
    DI_part = DI[start_index:end_index]

    x = np.arange(start_pos,end_pos,bin_size)
    width = bin_size/1.5
    pos_DI = np.array(DI_part)
    neg_DI = np.array(DI_part)
    pos_DI[pos_DI <= 0] = np.nan
    neg_DI[neg_DI > 0] = np.nan
    
    if input_file_hmm == '':
        print "Plotting DI values..."

        plt.close("all")
        plt.bar(x, pos_DI, width, color="r", label="Positive DI")
        plt.bar(x, neg_DI, width, color="g", label="Negative DI")
        plt.xlim([x[0]-bin_size*8,x[-1]+bin_size*8])
        plt.ylim([min(DI_part)-25,max(DI_part)+25])
        plt.title("Directionality Index " + species + " [Chr " + a_chr +": " + str(start_pos) + "-" + str(end_pos) + "]")
        plt.xlabel("Base coordinates")
        plt.ylabel("Directionality Index (DI) values")
        plt.grid(plot_grid)
        if plot_legend == True:
            plt.legend()
        plt.savefig("HiCtool_chr" + a_chr + "_DI.pdf", format = 'pdf')
        print "Done!"
    
    else:
        print "Plotting DI and true DI values..."     
        
        if isinstance(input_file_hmm,str):
            likelystates = load_hmm_states(input_file_hmm)
        else:
            likelystates = input_file_hmm
    
        DI_true = []
        for i in range(0,len(likelystates)):
            if likelystates[i] == 1:
                DI_true.append(min(DI_part)-12)
            elif likelystates[i] == 2:
                DI_true.append(min(DI_part)-15)
            else:
                DI_true.append(0)
        
        DI_true_part = DI_true[start_index:end_index]
        
        # Plot
        pos_DI_true = np.array(DI_true_part)
        neg_DI_true = np.array(DI_true_part)
        pos_DI_true[pos_DI_true != min(DI_part)-12] = np.nan
        neg_DI_true[neg_DI_true != min(DI_part)-15] = np.nan
        
        plt.close("all")
        plt.bar(x, pos_DI, width, color="r", label="Positive DI")
        plt.bar(x, neg_DI, width, color="g", label="Negative DI")
        plt.plot(x, pos_DI_true, marker=">", color="r", label="Positive true DI")
        plt.plot(x, neg_DI_true, marker="<", color="g", label="Negative true DI")
        plt.xlim([x[0]-bin_size*8,x[-1]+bin_size*8])
        plt.ylim([min(DI_part)-25,max(DI_part)+25])
        plt.title("Directionality Index " + species + " [Chr " + a_chr +": " + str(start_pos) + "-" + str(end_pos) + "]")
        plt.xlabel("Base coordinates")
        plt.ylabel("Directionality Index (DI) values")
        plt.grid(plot_grid)
        if plot_legend == True:
            plt.legend()
        plt.savefig("HiCtool_chr" + a_chr + "_DI_full.pdf", format = 'pdf')
        print "Done!"

def save_topological_domains(a_matrix, output_file):
    """
    Function to save the topological domains coordinates to text file.
    Each topological domain coordinates (start and end) occupy one row and are
    tab separated.
    Parameters:
        a_matrix (numpy matrix): file to be saved with topological domains coordinates.
        output_file: output file name in txt format.
    """
    def compile_row_string(a_row):
        return str(a_row).strip(']').strip('[').lstrip().replace(' ','\t')
    with open(output_file, 'w') as f:
        for row in a_matrix:
            f.write(compile_row_string(row)+'\n')


def load_topological_domains(input_file):
    """
    Function to load the topological domains coordinates from txt file.
    Parameters:
        input_file: input file name generated with "calculate_topological_domains" in txt format.
    Returns:
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
    

def calculate_chromosome_topological_domains(input_file_hmm,
                                             a_chr):
    """
    Function to calculate the topological domains coordinates of a chromosome. 
    Topological domains are stored in each line with tab separated start and end coordinates.
    Parameters:
        input_file_hmm (str | obj): txt file of the true DI values generated with the function "calculate_chromosome_true_DI" or
        object with the true DI values returned by "calculate_chromosome_true_DI.
        a_chr (str): chromosome number (example for chromosome 1: '1').
    Returns:
        List of lists with topological domain coordinates.
    """
    import numpy as np    
    
    print "Calculating topological domain coordinates..."
    if isinstance(input_file_hmm,str):
        likelystates = load_hmm_states(input_file_hmm)
    else:
        likelystates = input_file_hmm
    
    bin_size = 40000    
    
    # Start coordinates of the domains
    p = []
    for i in range(1,len(likelystates)):
        if (likelystates[i] == 1 and likelystates[i-1] == 2) or (likelystates[i] == 1 and likelystates[i-1] == 0):
            p.append(i * bin_size)
    
    # End coordinates of the domains
    n = []
    for i in range(1,len(likelystates)-1):
        if (likelystates[i] == 2 and likelystates[i+1] == 1) or (likelystates[i] == 2 and likelystates[i+1] == 0):
            n.append(i * bin_size)
    
    if len(p) == 0 or len(n) == 0:
        print "WARNING! No topological domains can be detected in chromosome " + a_chr
        return
    
    p1 = 0
    n1 = 0
    p2 = 1
    n2 = 1
    
    # Step 1: checking if the first negative values are greater than the first positive value.
    while n[n1] < p[p1]:
        n1 = n1 + 1
        n2 = n2 + 1
    
    # Now we have removed all the first negative values before the first positive one.
    topological_domains = []
    while p1 < len(p)-1 and n1 < len(n)-1:
        # Step 2: checking if there are two consecutive positive values.
        while n[n1] > p[p2] and p2 < len(p)-1:
            p2 = p2 + 1
        # Now we have removed the possible gaps between consecutive positive states.
    
        # Step 3: checking if there are two consecutive negative values.
        while n[n2] < p[p2] and n2 < len(n)-1:
            n1 = n1 + 1
            n2 = n2 + 1
        # Now we have removed the possible gaps between consecutive negative states.
    
        # Step 4: identification of the Topological Domain.
        topological_domains.append([p[p1],n[n1]])
        p1 = p2
        n1 = n2
        p2 = p1 + 1
        n2 = n1 + 1
    
    save_topological_domains(np.matrix(topological_domains),"HiCtool_chr" + a_chr + "_topological_domains.txt")
    print "Done!"
    return topological_domains


def compute_full_tad_analysis(input_contact_matrix,
                              a_chr,
                              isGlobal,
                              tab_sep=False,
                              species='hg38',
                              custom_species_sizes={},
                              sexual_chromosomes=[],
                              save_di=False,
                              save_hmm=False):
    """
    Compute DI values, HMM states and topological domain coordinates for a chromosome a_chr.
    Parameters:
        input_contact_matrix (str | obj): normalized intra-chromosomal contact matrix at a bin size of 40kb passed as a filename (str)
        or an object. Either a single contact matrix or a global contact matrix can be passed (see following parameters).
        a_chr: chromosome number (example for chromosome 1: '1').
        isGlobal (bool): set True if your input matrix is a global matrix.
        tab_sep (bool): set True if your input matrix is in a tab separated format. If the matrix is passed as an
        object, this parameter is not taken into consideration.
        species (str): 'hg38' or 'mm10' or any other species label in string format.
        custom_species_sizes (dict): dictionary containing the sizes of the chromosomes
        of your custom species. The keys of the dictionary are chromosomes in string
        format (example for chromosome 1: '1'), the values are chromosome lengths as int.
        sexual_chromosomes (list): list of the sexual chromosomes (if present) in your
        custom species (example for chromosome X: 'X').
        save_di (bool): if True, save the DI values to txt file.
        save_hmm (bool): if True, save the HMM states to txt file.
    Output: 
        txt file containing topological domains coordinates.
    """        
    import copy
    
    if isGlobal == False:
        if isinstance(input_contact_matrix, str):
            if tab_sep == False:
                contact_matrix = load_matrix(input_contact_matrix)
            else:
                contact_matrix = load_matrix_tab(input_contact_matrix)
        else:
            contact_matrix = copy.deepcopy(input_contact_matrix)
    else:
        contact_matrix = extract_single_map(input_global_matrix=input_contact_matrix,
                                            tab_sep=tab_sep,
                                            chr_row=a_chr,
                                            chr_col=a_chr,
                                            species=species,
                                            bin_size=40000,
                                            data_type='',
                                            custom_species_sizes=custom_species_sizes,
                                            sexual_chromosomes=sexual_chromosomes,
                                            save_output=False,
                                            save_tab=False)

    # DI VALUES
    DI = calculate_chromosome_DI(input_contact_matrix=contact_matrix,
                                 a_chr=a_chr,
                                 isGlobal=False,
                                 tab_sep=False,
                                 species=species,
                                 custom_species_sizes=custom_species_sizes,
                                 sexual_chromosomes=sexual_chromosomes,
                                 save_file=save_di)
    
    # HMM STATES
    HMM = calculate_chromosome_true_DI(DI,a_chr,save_hmm)
    
    # TOPOLOGICAL DOMAIN COORDINATES
    tad = calculate_chromosome_topological_domains(HMM,a_chr)
    return tad
"""
Program to:
1) Calculate and plot the observed DI and true DI (Hidden Markov Model).
2) Calculate the topological domains coordinates.

To use this code, normalized contact data must be provided.
"""

def save_list(a_list, output_file):
    """
    Save a list in a .txt file.
    Parameters:
        list_: name of the list to save.
        output_file: output file name in txt format.
    Output: 
        .txt file containing the saved list.       
    """
    with open (output_file,'w') as fout:
        n = len(a_list)
        for i in xrange(n):
            fout.write('%s\n' %a_list[i])
            
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


def calculate_chromosome_DI(a_chr, species="hg38"):
    """
    Function to calculate the DI values for a chromosome of a species and save them 
    in a .txt file.
    Parameters:
        a_chr (str): chromosome number (example for chromosome 1: '1').
        species (str): hg38 or mm10
    Returns:
        List with the DI values.
    """
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

    bin_size = 40000
    start_pos = 0
    end_pos = (chromosomes[species][a_chr]/1000000)*1000000
    
    start_part = str(float(start_pos)/float(1000000))
    end_part = str(float(end_pos)/float(1000000))
    binsize_str = str(float(bin_size)/float(1000000))
                          
    contact_matrix = load_matrix('HiCtool_normalized_fend_contact_matrix_chr' + a_chr + '_' + binsize_str + 'mb_' + start_part + 'mb_' + end_part + 'mb.txt')
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
    
        try:
            di = ((B-A)/(abs(B-A)))*((((A-E)**2)/E)+(((B-E)**2)/E))
        except ZeroDivisionError:
            di = 0
    
        DI.append(di)

    save_list(DI,'HiCtool_chr' + a_chr + '_DI.txt')
    
    return DI
    

def plot_chromosome_DI(a_chr, start_pos, end_pos, species='hg38', show=False):
    """
    Function to plot the DI values for a chromosome.
    Parameters:
        a_chr (str): chromosome number (example for chromosome 1: '1').
        start_pos (int): start coordinate for the plot.
        end_pos (int): end coordinate for the plot.
        species (str): hg38 or mm10
        show (bool): if true, the figure is shown. If false (default),
        the figure is saved to file.
    """    
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib
    matplotlib.use('Agg')
    
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
    
    myfile = "HiCtool_chr" + a_chr + "_DI.txt"
    fp = open(myfile,'r+')
    lines = fp.read().split('\n')
    lines = lines[:-1]
    A = (np.nan_to_num(np.array(map(float, lines)))).tolist()
    
    # Selecting the region of interest for observed DI plot
    bin_size = 40000
    if end_pos > (chromosomes[species][a_chr]/1000000)*1000000 and end_pos < chromosomes[species][a_chr]:
        end_pos = (chromosomes[species][a_chr]/1000000)*1000000
    elif end_pos > chromosomes[species][a_chr]:
        print("ERROR: end coordinate exceeds chromosome dimension")
    start_index = int(round(start_pos/bin_size))
    end_index = int(round((end_pos)/bin_size))
    DI_part = A[start_index:end_index]
    
    # Plot
    x = np.arange(start_pos,end_pos,bin_size)
    width = bin_size/1.5
    pos_DI = np.array(DI_part)
    neg_DI = np.array(DI_part)
    pos_DI[pos_DI <= 0] = np.nan
    neg_DI[neg_DI > 0] = np.nan
    
    plt.close("all")
    plt.bar(x, pos_DI, width, color="r", label="Positive observed DI")
    plt.bar(x, neg_DI, width, color="g", label="Negative observed DI")
    plt.xlim([x[0]-bin_size*8,x[-1]+bin_size*8])
    plt.ylim([min(DI_part)-25,max(DI_part)+25])
    plt.title("Directionality Index " + species + " [Chr " + a_chr +": " + str(start_pos) + "-" + str(end_pos) + "]")
    plt.xlabel("Base coordinates")
    plt.grid(True)
    plt.legend()
    if show==True:
        plt.show()
    else:
        plt.savefig("HiCtool_chr" + a_chr + "_DI.png")


def calculate_chromosome_true_DI(a_chr):
    """
    Function to calculate the true DI values for a chromosome and save 
    them in a .txt file.
    Input:
        a_chr (str): chromosome number (example for chromosome 1: '1').
    Returns:
        List with the true DI values (HMM states).
    """
    import numpy as np
    import hmmlearn.hmm as hmm
    
    myfile = "HiCtool_chr" + a_chr + "_DI.txt"
    fp = open(myfile,'r+')
    lines = fp.read().split('\n')
    lines = lines[:-1]
    A = (np.nan_to_num(np.array(map(float, lines)))).tolist()
    
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
    save_list(likelystates, "HiCtool_chr" + a_chr + "_hmm_states.txt")
    
    return(likelystates)


def plot_chromosome_DI_full(a_chr, start_pos, end_pos, species='hg38', show=False):
    """
    Function to plot the DI values (observed and true DI states) for a chromosome.
    Parameters:
        a_chr (str): chromosome number (example for chromosome 1: '1').
        start_pos (int): start coordinate for the plot.
        end_pos (int): end coordinate for the plot.
        species (str): hg38 or mm10
        show (bool): if false (default), the figure is saved to file. If true, 
        the figure is shown.
    """    
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib
    matplotlib.use('Agg')
    
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
    
    myfile1 = "HiCtool_chr" + a_chr + "_DI.txt"
    fp = open(myfile1,'r+')
    lines = fp.read().split('\n')
    lines = lines[:-1]
    DI = (np.nan_to_num(np.array(map(float, lines)))).tolist()
    
    myfile2 = "HiCtool_chr" + a_chr + "_hmm_states.txt"
    fp = open(myfile2,'r+')
    lines = fp.read().split('\n')
    lines = lines[:-1]
    likelystates = map(int,lines)
    
    # Selecting the region of interest for observed DI plot
    bin_size = 40000
    if end_pos > (chromosomes[species][a_chr]/1000000)*1000000 and end_pos < chromosomes[species][a_chr]:
        end_pos = (chromosomes[species][a_chr]/1000000)*1000000
    elif end_pos > chromosomes[species][a_chr]:
        print("ERROR: end coordinate exceeds chromosome dimension")
    start_index = int(round(start_pos/bin_size))
    end_index = int(round((end_pos)/bin_size))
    DI_part = DI[start_index:end_index]
    
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
    x = np.arange(start_pos,end_pos,bin_size)
    width = bin_size/1.5
    pos_DI = np.array(DI_part)
    neg_DI = np.array(DI_part)
    pos_DI[pos_DI <= 0] = np.nan
    neg_DI[neg_DI > 0] = np.nan
    
    pos_DI_true = np.array(DI_true_part)
    neg_DI_true = np.array(DI_true_part)
    pos_DI_true[pos_DI_true != min(DI_part)-12] = np.nan
    neg_DI_true[neg_DI_true != min(DI_part)-15] = np.nan
    
    plt.close("all")
    plt.bar(x, pos_DI, width, color="r", label="Positive observed DI")
    plt.bar(x, neg_DI, width, color="g", label="Negative observed DI")
    plt.plot(x, pos_DI_true, marker=">", color="r", label="Positive true DI")
    plt.plot(x, neg_DI_true, marker="<", color="g", label="Negative true DI")
    plt.xlim([x[0]-bin_size*8,x[-1]+bin_size*8])
    plt.ylim([min(DI_part)-25,max(DI_part)+25])
    plt.title("Directionality Index " + species + " [Chr " + a_chr +": " + str(start_pos) + "-" + str(end_pos) + "]")
    plt.xlabel("Base coordinates")
    plt.grid(True)
    plt.legend()
    if show==True:
        plt.show()
    else:
        plt.savefig("HiCtool_chr" + a_chr + "_DI_full.png")


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
        input_file: input file name in txt format.
    Returns:
        List of topological domain coordinates.
    """
    import csv
    with open(input_file, 'r') as f:
        reader = csv.reader(f, dialect='excel', delimiter='\t')
        topological_domains = []
        for row in reader:
            row_int = [int(x) for x in row]
            topological_domains.append(row_int)
        return topological_domains


def calculate_chromosome_topological_domains(a_chr):
    """
    Function to calculate the topological domains coordinates of a chromosome.
    Parameters:
        a_chr (str): chromosome number (example for chromosome 1: '1').
    Returns:
        List of topological domains.
    """
    import numpy as np    
    
    myfile = "HiCtool_chr" + a_chr + "_hmm_states.txt"
    fp = open(myfile,'r+')
    lines = fp.read().split('\n')
    lines = lines[:-1]
    likelystates = map(int,lines)
    bin_size = 40000    
    
    # Start coordinates of the domains
    p = []
    for i in range(1,len(likelystates)-1):
        if (likelystates[i] == 1 and likelystates[i-1] == 2) or (likelystates[i] == 1 and likelystates[i-1] == 0):
            p.append(i * bin_size)
    
    # End coordinate of the domains
    n = []
    for i in range(1,len(likelystates)-1):
        if (likelystates[i] == 2 and likelystates[i-1] == 1) or (likelystates[i] == 2 and likelystates[i-1] == 0):
            n.append(i * bin_size)
    
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

    return topological_domains

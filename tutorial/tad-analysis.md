# TAD analysis

This pipeline illustrates the procedure to calculate topologically associated domain (TAD) coordinates and visualize TADs on the heatmap. Topological domain coordinates should be calculated on the normalized data (as we do here), but if you wish you can also use this pipeline to calculate TAD coordinates on the observed data.

- Section 1 is a synthetic pipeline to compute TAD coordinates and plot TADs on the heatmap.
- Section 2 gives additional details about the analysis performed in Section 1.

## Table of Contents

1. [Performing TAD analysis](#1-performing-tad-analysis)
   - [1.1. Data normalized with Yaffe-Tanay](#11-plotting-the-observed-and-true-di-values)
   - [1.2. Data normalized with Hi-Corrector](#12-performing-tad-analysis-with-parallel-processing)
   - [1.3. Plotting TADs over the heatmap]()
2. [Supplementary TAD analysis]
   - [2.1. Calculating and plotting the observed DI]
   - [2.2. Calculating and plotting the true DI]
   - [2.3. Calculating TAD coordinates]
   
## 1. Performing TAD analysis

TAD coordinates are calculated using the shifts of the true Directionality Index (DI) as described by [Dixon et al., (2012)](http://www.nature.com/nature/journal/v485/n7398/abs/nature11082.html). DI values and true DI values (HMM biased states) are calculated following the same paper as well. More details about this are reported in section 2.

### 1.1. Data normalized with Yaffe-Tanay

To calculate TAD coordinates for a chromosome (here chr 6) use the function ``compute_full_tad_analysis`` of [HiCtool_TAD.py](/scripts/HiCtool_TAD.py).
```Python
execfile('HiCtool_TAD.py')

tad_coord = compute_full_tad_analysis('HiCtool_chr6_40kb_normalized_fend.txt', a_chr='6',
                                      isGlobal=False, tab_sep=False, species='hg38',
                                      save_di=True, save_hmm=True)
```
``tad_coord`` is a list of topological domains. Each topological domain is a list of two elements that are the start and end coordinate of the domain. ``save_di`` and ``save_hmm`` set to ``True`` allow to save also the DI values and HMM biased states.

To calculate the **topological domain coordinates for all the chromosomes at once** you may use an approach as the following:
```Python
my_chromosomes = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y']
my_tad_coord = {} # dictionary to store the tad coordinates of different chromosomes
for i in my_chromosomes:
    my_tad_coord[i] = compute_full_tad_analysis('HiCtool_chr' + i + '_40kb_normalized_fend.txt', 
                                                a_chr=i, isGlobal=False, tab_sep=False, 
                                                species='hg38', save_di=True, save_hmm=True)
```

### 1.2. Data normalized with Hi-Corrector

To calculate TAD coordinates for a chromosome (here chr 6) use the function ``compute_full_tad_analysis`` of [HiCtool_TAD.py](/scripts/HiCtool_TAD.py).

Given the big dimension of the global normalized contact matrix, it is suggested to load it to your workspace ahead of the analysis. This will avoid several reloadings of the contact matrix and save time.
```Python
execfile('HiCtool_TAD.py')
global_normalized_40kb = load_matrix_tab('output_ic_mes/output_normalized.txt')

tad_coord = compute_full_tad_analysis(global_normalized_40kb, a_chr='Y', isGlobal=True,
                                      species='hg38', save_di=True, save_hmm=True)
```
``tad_coord`` is a list of topological domains. Each topological domain is a list of two elements that are the start and end coordinate of the domain. ``save_di`` and ``save_hmm`` set to ``True`` allow to save also the DI values and HMM biased states.
   
To calculate the **topological domain coordinates for all the chromosomes at once** you may use an approach as the following:
```Python
my_chromosomes = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y']
my_tad_coord = {} # dictionary to store the tad coordinates of different chromosomes
for i in my_chromosomes:
    my_tad_coord[i] = compute_full_tad_analysis(global_normalized_40kb, a_chr=i, isGlobal=True,
                                                species='hg38', save_di=True, save_hmm=True)
```

### 1.3. Plotting TADs over the heatmap



plot_chromosome_data(fend_normalized_chr6, a_chr='6', bin_size=40000, full_matrix=False, start_coord=50000000, end_coord=54000000, species='hg38', data_type="normalized_fend", my_colormap=['white', 'red'], cutoff_type='percentile', cutoff=95, max_color='#460000', my_dpi=1000, plot_domains=True)












### 1.1. Plotting the observed and true DI values

Observed DI values and HMM states can be also calculated and plotted separately.

To **calculate the DI values** use the function ``calculate_chromosome_DI`` as following:
```Python
#
execfile('HiCtool_tad.py')

# Yaffe and Tanay normalization method
DI_chr6 = calculate_chromosome_DI(input_contact_matrix='HiCtool_chr6_40kb_normalized_fend.txt', 
                                  a_chr='6', isGlobal=False, tab_sep=False)

# Hi-Corrector normalization method (global matrix already loaded above)
DI_chr6 = calculate_chromosome_DI(input_contact_matrix=global_normalized_40kb, 
                                  a_chr='6', isGlobal=True)
```
Previously calculated DI values and saved to file can be loaded using the function ``load_DI_values``:
```Python
DI_chr6 = load_DI_values('HiCtool_chr6_DI.txt')
```
**DI values can be plotted** using the function ``plot_chromosome_DI``:
```Python
execfile('HiCtool_tad.py')

plot_chromosome_DI(DI_chr6, a_chr='6', start_pos=50000000, end_pos=54000000)
```

![](/figures/HiCtool_chr6_DI.png)

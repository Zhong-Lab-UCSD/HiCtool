# TAD analysis

This pipeline illustrates the procedure to calculate topologically associated domain (TAD) coordinates and visualize TADs. Topological domain coordinates should be calculated on the normalized data (as we do here) but, if you wish, you can also use this pipeline to calculate TAD coordinates on the observed data. Default species that you can use are human *hg38* and mouse *mm10*. For other species use the specific function arguments (see [API documentation](https://sysbio.ucsd.edu/public/rcalandrelli/HiCtool_API_documentation.pdf)).

For more information about the Python functions used here check the [API documentation](https://sysbio.ucsd.edu/public/rcalandrelli/HiCtool_API_documentation.pdf).

**Note!** Contact data must be at 40 kb resolution to perform TAD analysis!

- [Section 1 (Performing TAD analysis)](#1-performing-tad-analysis) is a pipeline to compute TAD coordinates and plot TADs on the heatmap.
- [Section 2 (Supplementary TAD analysis)](#2-supplementary-tad-analysis) gives additional details about the analysis performed in Section 1.

## Table of Contents

1. [Performing TAD analysis](#1-performing-tad-analysis)
   - [1.1. Data normalized with Yaffe-Tanay method](#11-data-normalized-with-yaffe-tanay-method)
   - [1.2. Data normalized with Hi-Corrector](#12-data-normalized-with-hi-corrector)
2. [Supplementary TAD analysis](#2-supplementary-tad-analysis)
   - [2.1. Calculating and plotting the observed DI](#21-calculating-and-plotting-the-observed-di)
   - [2.2. Calculating and plotting the true DI](#22-calculating-and-plotting-the-true-di)
   - [2.3. Calculating TAD coordinates](#23-calculating-tad-coordinates)
   
## 1. Performing TAD analysis

TAD coordinates, as well as DI values and true DI values (HMM states), are calculated as described by [Dixon et al., (2012)](http://www.nature.com/nature/journal/v485/n7398/abs/nature11082.html). More details about this are reported in [Supplementary TAD analysis](#2-supplementary-tad-analysis).

**Note!** Contact data must be at 40 kb resolution to perform TAD analysis!

### 1.1. Data normalized with Yaffe-Tanay method

Use these instructions if you normalized your data using the method from [Yaffe and Tanay](https://github.com/Zhong-Lab-UCSD/HiCtool/blob/master/tutorial/normalization-yaffe-tanay.md).

To calculate TAD coordinates for a chromosome (here for chr 6) use the function ``compute_full_tad_analysis`` of [HiCtool_TAD_analysis.py](/scripts/HiCtool_TAD_analysis.py) as following:
```Python
execfile('HiCtool_TAD_analysis.py')

tad_coord = compute_full_tad_analysis('HiCtool_chr6_40kb_normalized_fend.txt', a_chr='6',
                                      isGlobal=False, tab_sep=False, species='hg38',
                                      save_di=True, save_hmm=True)
```
``tad_coord`` is a list of topological domains. Each topological domain is a list of two elements that are the start and end coordinate of the domain. ``save_di`` and ``save_hmm`` set to ``True`` allow to save also the DI values and HMM biased states to txt file.

**Note!** The end coordinate of each domain is saved as the start position of the last bin (40 kb) belonging to each domain.

To calculate the **topological domain coordinates for all the chromosomes at once** you may use an approach as the following:
```Python
my_chromosomes = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y']
my_tad_coord = {} # dictionary to store the tad coordinates of different chromosomes
for i in my_chromosomes:
    my_tad_coord[i] = compute_full_tad_analysis('HiCtool_chr' + i + '_40kb_normalized_fend.txt', 
                                                a_chr=i, isGlobal=False, tab_sep=False, 
                                                species='hg38', save_di=True, save_hmm=True)
```
To plot the topological domains on the heatmap, use the function ``plot_chromosome_data`` of [HiCtool_yaffe_tanay.py](/scripts/HiCtool_yaffe_tanay.py) and pass the topological domain coordinates with the argument ``topological_domains``. Here we plot the heatmap for chr6: 80,000,000-120,000,000. The color of the domains can be changed by using the parameter ``domain_color``.
```Python
execfile('HiCtool_yaffe_tanay.py')
plot_chromosome_data('HiCtool_chr6_40kb_normalized_fend.txt', 
                     a_chr='6', bin_size=40000, full_matrix=False, 
                     start_coord=80000000, end_coord=120000000, 
                     species='hg38', 
                     data_type="normalized_fend", 
                     my_colormap=['white', 'red'], 
                     cutoff_type='percentile', cutoff=95, max_color='#460000', 
                     my_dpi=1000, 
                     topological_domains='HiCtool_chr6_topological_domains.txt',
                     domain_color='#0000ff')
```
![](/figures/HiCtool_chr6_40kb_80-120mb_normalized_fend_domains.png)

Zoom in on a smaller region (chr6: 50,000,000-54,000,000):
```Python
plot_chromosome_data('HiCtool_chr6_40kb_normalized_fend.txt', 
                     a_chr='6', bin_size=40000, full_matrix=False, 
                     start_coord=50000000, end_coord=54000000, 
                     species='hg38', 
                     data_type="normalized_fend", 
                     my_colormap=['white', 'red'], 
                     cutoff_type='percentile', cutoff=95, max_color='#460000', 
                     my_dpi=1000, 
                     topological_domains='HiCtool_chr6_topological_domains.txt',
                     domain_color='#0000ff')
```
![](/figures/HiCtool_chr6_40kb_50-54mb_normalized_fend_domains.png)


### 1.2. Data normalized with Hi-Corrector

Use these instructions if you normalized your data using the [Hi-Corrector](https://github.com/Zhong-Lab-UCSD/HiCtool/blob/master/tutorial/normalization-matrix-balancing.md).

To calculate TAD coordinates for a chromosome (here for chr 6) use the function ``compute_full_tad_analysis`` of [HiCtool_TAD_analysis.py](/scripts/HiCtool_TAD_analysis.py).

**Tip!** Given the big dimension of the global normalized contact matrix, it is suggested to load it to your workspace ahead of the analysis and then always point to this object instead of the matrix file. This will avoid several reloadings of the contact matrix and save time.
```Python
execfile('HiCtool_TAD_analysis.py')
global_normalized_40kb = load_matrix_tab('output_ic_mes/output_normalized.txt')

tad_coord = compute_full_tad_analysis(global_normalized_40kb, a_chr='6', isGlobal=True,
                                      species='hg38', save_di=True, save_hmm=True)
```
``tad_coord`` is a list of topological domains. Each topological domain is a list of two elements that are the start and end coordinate of the domain. ``save_di`` and ``save_hmm`` set to ``True`` allow to save also the DI values and HMM biased states to txt file.

**Note!** The end coordinate of each domain is saved as the start position of the last bin (40 kb) belonging to each domain.
   
To calculate the **topological domain coordinates for all the chromosomes at once** you may use an approach as the following:
```Python
my_chromosomes = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y']
my_tad_coord = {} # dictionary to store the tad coordinates of different chromosomes
for i in my_chromosomes:
    my_tad_coord[i] = compute_full_tad_analysis(global_normalized_40kb, a_chr=i, isGlobal=True,
                                                species='hg38', save_di=True, save_hmm=True)
```
To plot the topological domains on the heatmap, use the function ``plot_map`` of [HiCtool_full_map_analysis.py](/scripts/HiCtool_full_map_analysis.py) and pass the topological domain coordinates with the argument ``topological_domains``. Here we plot the heatmap for chr6: 80,000,000-120,000,000. The color of the domains can be changed by using the parameter ``domain_color``.
```Python
execfile('HiCtool_full_map_analysis.py')
plot_map(input_matrix=global_normalized_40kb, isGlobal=True, tab_sep=True,
         chr_row='6', chr_col='6', bin_size=40000, 
         chr_row_coord=[80000000,120000000], chr_col_coord=[80000000,120000000],
         data_type="normalized", species='hg38',
         my_colormap=['white', 'red'],
         cutoff_type='perc', cutoff=95, max_color='#460000',
         topological_domains='HiCtool_chr6_topological_domains.txt',
         domain_color='#0000ff')
```
![](/figures/HiCtool_chr6_40kb_80-120mb_normalized_domains.png)

Zoom in on a smaller region (chr6: 50,000,000-54,000,000):
```Python
plot_map(input_matrix=global_normalized_40kb, isGlobal=True, tab_sep=True,
         chr_row='6', chr_col='6', bin_size=40000, 
         chr_row_coord=[50000000,54000000], chr_col_coord=[50000000,54000000],
         data_type="normalized", species='hg38',
         my_colormap=['white', 'red'],
         cutoff_type='perc', cutoff=95, max_color='#460000',
         topological_domains='HiCtool_chr6_topological_domains.txt',
         domain_color='#0000ff')
```
![](/figures/HiCtool_chr6_40kb_50-54mb_normalized_domains.png)


## 2. Supplementary TAD analysis

This section illustrates details about the calculation of the DI values, true DI values (HMM states) and topological domains coordinates ([Dixon et al., (2012)](http://www.nature.com/nature/journal/v485/n7398/abs/nature11082.html)).


### 2.1. Calculating and plotting the observed DI

Regions at the periphery of the topological domains are highly biased in their interaction frequencies: the most upstream portion of a topological domain is highly biased towards interacting downstream, and the downstream portion of a topological domain is highly biased towards interacting upstream. To determine the directional bias at any given bin in the genome the Directionality Index (DI) is used, quantifying the degree of upstream or downstream bias of a given bin:

This is the formula used to calculate the DI:

![](/figures/DI_formula.png)

where:

- *A* is the the number of reads that map from a given 40 kb bin to the upstream 2 Mb.
- *B* is the the number of reads that map from a given 40 kb bin to the downstream 2 Mb.
- *E* is the expected number of contacts for each bin and it equals (A+B)/2.

To compute DI, we need the fend normalized contact data at a bin size of 40 kb. Since the bin size is 40 kb (here we used 40 kbpb which stands for 40 kb per bin), hence the detection region of upstream or downstream biases 2 Mb is converted to 50 bins (2 Mb / 40 kbpb = 50 bins).

Observed DI values and HMM states can be also calculated and plotted separately.

To **calculate the DI values** use the function ``calculate_chromosome_DI`` of [HiCtool_TAD_analysis.py](/scripts/HiCtool_TAD_analysis.py) as following:
```Python
execfile('HiCtool_TAD_analysis.py')

# Yaffe and Tanay normalization method
DI_chr6 = calculate_chromosome_DI(input_contact_matrix='HiCtool_chr6_40kb_normalized_fend.txt', 
                                  a_chr='6', isGlobal=False, tab_sep=False,
                                  species='hg38', save_file=True)

# Hi-Corrector normalization method
global_normalized_40kb = load_matrix_tab('output_ic_mes/output_normalized.txt')
DI_chr6 = calculate_chromosome_DI(input_contact_matrix=global_normalized_40kb, 
                                  a_chr='6', isGlobal=True, species='hg38', save_file=True)
```
Previously calculated DI values and saved to file can be loaded using the function ``load_DI_values``:
```Python
DI_chr6 = load_DI_values('HiCtool_chr6_DI.txt')
```
**DI values can be plotted** using the function ``plot_chromosome_DI``:
```Python
plot_chromosome_DI(DI_chr6, a_chr='6', full_chromosome=False,
                   start_pos=50000000, end_pos=54000000, species='hg38',
                   plot_legend=True, plot_grid=True)
```
![](/figures/HiCtool_chr6_DI.png)

### 2.2. Calculating and plotting the true DI

We use a Hidden Markov Model (HMM) based on the Directionality Index to identify biased states (true DI).

For true DI calculation, we consider the **emission sequence** as the observed DI values and the Transition Matrix, Emission Matrix and initial State Sequence as unknown. We have **three emissions** 1, 2, 0 corresponding to a positive (1), negative (2) or zero (0) value of the observed DI. In our analysis, we associate to the emission '0' all the absolute DI values under a threshold of 0.4. So, first we estimate the model parameters and then the most probable sequence of states using the Viterbi algorithm. 

All these steps are performed using the function ``calculate_chromosome_hmm_states`` of [HiCtool_TAD_analysis.py](/scripts/HiCtool_TAD_analysis.py), which calculates true DI values and save the output to a txt file:
```Python
execfile('HiCtool_TAD_analysis.py')
hmm_chr6 = calculate_chromosome_hmm_states(input_file_DI='HiCtool_chr6_DI.txt', 
                                           a_chr='6',
                                           save_file=True)
```
Previously calculated true DI values and saved to file can be loaded using the function ``load_hmm_states``:
```Python
hmm_chr6 = load_hmm_states('HiCtool_chr6_hmm_states.txt')
```
**DI and true DI values can be plotted** using the function ``plot_chromosome_DI`` and passing also ``hmm_chr6`` (or ``'HiCtool_chr6_hmm_states.txt'``) as input parameter:
```Python
plot_chromosome_DI(input_file_DI=DI_chr6, a_chr='6', full_chromosome=False,
                   start_pos=50000000, end_pos=54000000, input_file_hmm=hmm_chr6,
                   species='hg38', plot_legend=True, plot_grid=True)
```
![](/figures/HiCtool_chr6_DI_full.png)

### 2.3. Calculating TAD coordinates

The true DI values allow to infer the locations of the topological domains in the genome. A domain is initiated at the beginning of a single downstream biased HMM state (red color in the above figure). The domain is continuous throughout any consecutive downstream biased state. The domain will then end when the last in a series of upstream biased states (green color in the above figure) is reached, with the domain ending at the end of the last HMM upstream biased state.

To calculate the topological domains coordinates, first we extract all the potential start and end coordinates according to the definition, and then we evaluate a list of conditions to take into account the possible presence of gaps between a series of positive or negative states values. The figure below shows a summary of the procedure:

![](/figures/HiCtool_topological_domains_flowchart.png)

Topological domain coordinates are calculated using the function ``calculate_chromosome_topological_domains`` of [HiCtool_TAD_analysis.py](/scripts/HiCtool_TAD_analysis.py) and passing as input the HMM states:
```Python
execfile('HiCtool_TAD_analysis.py')
domains_chr6 = calculate_chromosome_topological_domains(input_file_hmm='HiCtool_chr6_hmm_states.txt', a_chr='6')
```
Topological domain coordinates are saved in a tab separated file. Each line is a topological domain, with two elements that are the start and end coordinate of the domain.
Previously calculated topological domain coordinates and saved to file can be loaded using the function ``load_topological_domains``:
```Python
domains_chr6 = load_topological_domains('HiCtool_chr6_topological_domains.txt')
```

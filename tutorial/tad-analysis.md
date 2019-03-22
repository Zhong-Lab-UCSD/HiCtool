# TAD analysis

This pipeline illustrates the procedure to calculate topologically associated domain (TAD) coordinates and visualize TADs on the heatmap. Topological domain coordinates should be calculated on the normalized data (as we do here), but if you wish you can also use this pipeline to calculate TAD coordinates on the observed data.

## Table of Contents

1. [Performing TAD analysis](#1-performing-tad-analysis)
   - [1.1. Plotting the observed and true DI values](#11-plotting-the-observed-and-true-di-values)
   - [1.2. Performing TAD analysis with parallel processing](#12-performing-tad-analysis-with-parallel-processing)
2. [Visualizing TADs on the heatmap](#2-visualizing-tads-on-the-heatmap)

## 1. Performing TAD analysis

TAD coordinates are calculated using the shifts of the true Directionality Index (DI) as described by [Dixon et al., (2012)](http://www.nature.com/nature/journal/v485/n7398/abs/nature11082.html). DI values and true DI values (HMM biased states) are calculated following the same paper as well.

To calculate TAD coordinates for a chromosome (here chr 6) use the function ``compute_full_tad_analysis`` of [HiCtool_tad.py](/scripts/HiCtool_tad.py).

- **If your data are normalized using the [Hi-Corrector approach](https://github.com/Zhong-Lab-UCSD/HiCtool/blob/master/tutorial/normalization-matrix-balancing.md)**, it is suggested to load to your workspace the global normalized matrix ahead of the analysis. The global matrix at 40 kb resolution is very big, this will avoid several reloadings of the contact matrix and save time.
    ```Python
   execfile('HiCtool_tad.py')
   global_normalized_40kb = load_matrix_tab('output_ic_mes/output_normalized.txt')

   tad_coord = compute_full_tad_analysis(global_normalized_40kb, a_chr='6', isGlobal=True,
                                         species='hg38', save_di=True, save_hmm=True)
   ```
   ``tad_coord`` is a list of topological domains. Each topological domain is a list of two elements that are the start and end coordinate of the domain. ``save_di`` and ``save_hmm`` set to ``True`` allow to save the DI values and HMM biased states also.
   
   To calculate the **topological domain coordinates for several chromosomes** you may use an approach as the following (in the example for chromosomes 1, 2, 6):
   ```Python
   my_chromosomes = ['1','2','6']
   my_tad_coord = {} # dictionary to store the tad coordinates of different chromosomes
   for i in my_chromosomes:
       my_tad_coord[i] = compute_full_tad_analysis(global_normalized_40kb, a_chr=i, isGlobal=True,
                                                   species='hg38', save_di=True, save_hmm=True)
   ```

- **If your data are normalized using the [Yaffe and Tanay approach](https://github.com/Zhong-Lab-UCSD/HiCtool/blob/master/tutorial/normalization-yaffe-tanay.md)**, you have normalized contact matrix per each single chromosome, therefore there is no need of loading them initially.
   ```Python
   execfile('HiCtool_tad.py')
   
   tad_coord = compute_full_tad_analysis('HiCtool_chr6_40kb_normalized_fend.txt', a_chr='6',
                                         isGlobal=False, tab_sep=False, species='hg38',
                                         save_di=True, save_hmm=True)
   ```
   ``tad_coord`` is a list of topological domains. Each topological domain is a list of two elements that are the start and end coordinate of the domain. ``save_di`` and ``save_hmm`` set to ``True`` allow to save the DI values and HMM biased states also.

   To calculate the **topological domain coordinates for several chromosomes** you may use an approach as the following (in the example for chromosomes 1, 2, 6):
   ```Python
   my_chromosomes = ['1','2','6']
   my_tad_coord = {} # dictionary to store the tad coordinates of different chromosomes
   for i in my_chromosomes:
       my_tad_coord[i] = compute_full_tad_analysis('HiCtool_chr' + i + '_40kb_normalized_fend.txt', 
                                                   a_chr=i, isGlobal=False, tab_sep=False, 
                                                   species='hg38', save_di=True, save_hmm=True)
   ```

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

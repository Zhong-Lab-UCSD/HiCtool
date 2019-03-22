# TAD analysis

This pipeline illustrates the procedure to calculate topologically associated domain (TAD) coordinates and visualize TADs on the heatmap. Topological domain coordinates should be calculated on the normalized data (as we do here), but if you wish you can also use this pipeline to calculate TAD coordinates on the observed data.

## Table of Contents

1. [Performing TAD analysis](#1-performing-tad-analysis)
   - [1.1. Calculating Directionality Index (DI)](#11-calculating-directionality-index)
   - [1.2. Calculating the true DI (HMM states)](#12-creating-the-fend-object)
   - [1.3. Plotting the observed and true DI](#13-plotting-the-observed-and-true-di)
   - [1.4. Performing TAD analysis with parallel processing](#14-performing-tad-analysis-with-parallel-processing)
2. [Visualizing TADs on the heatmap](#2-visualizing-tads-on-the-heatmap)

## 1. Performing TAD analysis

TAD coordinates are calculated using the shifts of the true Directionality Index (DI) as described by [Dixon et al., (2012)](http://www.nature.com/nature/journal/v485/n7398/abs/nature11082.html). DI values and true DI values (HMM biased states) are calculated follow the same paper as well.

To calculate TAD coordinates for a chromosome (here chr 6) use the function ``compute_full_tad_analysis`` of [HiCtool_tad.py](/scripts/HiCtool_tad.py).

- **If your data are normalized using the [Hi-Corrector approach](https://github.com/Zhong-Lab-UCSD/HiCtool/blob/master/tutorial/normalization-matrix-balancing.md)**, it is suggested to load to your workspace the global normalized matrix ahead of the analysis. The global matrix at 40 kb resolution is very big, this will avoid several reloadings of the contact matrix and save time.
    ```Python
   execfile('HiCtool_tad.py')
   global_normalized_40kb = load_matrix_tab('output_ic_mes/output_normalized.txt')

   tad_coord = compute_full_tad_analysis(global_normalized_40kb, a_chr='6', isGlobal=True,
                                         species='hg38', save_di=True, save_hmm=True)
   ```
   ``tad_coord`` is a list of topological domains. Each topological domain is a list of two elements that are the start and end coordinate of the domain. ``save_di`` and ``save_hmm`` set to ``True`` allow to save the DI values and HMM biased states also.
   
   To calculate the topological domain coordinates for several chromosomes you may use an approach as the following (in the example for chromosomes 1, 2, 6):
   ```Python
   my_chromosomes = ['1','2','6']
   my_tad_coord = {} # dictionary to save the tad coordinates of different chromosomes
   for i in my_chromosomes:
   my_tad_coord[i] = compute_full_tad_analysis(global_normalized_40kb, a_chr=i, isGlobal=True,
                                                species='hg38', save_di=True, save_hmm=True)
   ```

- **If your data are normalized using the [Yaffe and Tanay approach](https://github.com/Zhong-Lab-UCSD/HiCtool/blob/master/tutorial/normalization-yaffe-tanay.md)**, you have normalized contact matrix per each single chromosome, therefore there is no need of loading them initially.
   ```Python
   execfile('HiCtool_tad.py')
   
   tad_coord = compute_full_tad_analysis('HiCtool_chr6_40kb_normalized_fend.txt', a_chr='6',
                                         isGlobal=False, tab_sep=False, species='hg38', save_di=True, save_hmm=True)
   ```
   ``tad_coord`` is a list of topological domains. Each topological domain is a list of two elements that are the start and end coordinate of the domain. ``save_di`` and ``save_hmm`` set to ``True`` allow to save the DI values and HMM biased states also.

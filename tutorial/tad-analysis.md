# TAD analysis

This pipeline illustrates the procedure to calculate topologically associated domain (TAD) coordinates and visualize them on the heatmap.

## Table of Contents

1. [Performing TAD analysis]()
   - [1.1. Calculating Directionality Index (DI)](#1-running-hifive-functions)
   - [1.2. Calculating the true DI (HMM states)](#11-creating-the-fend-object)
   - [1.3. Plotting the observed and true DI](#21-single-processor-matrix-generation)
   - [1.4. Performing TAD analysis parallel processing]()
2. [Visualizing TAD on the heatmap](#2-generating-the-global-observed-contact-matrix)


```Python
compute_full_tad_analysis('HiCtool_40kb_matrix_global_normalized_tab.txt',
'6',
isGlobal=True,
tab_sep=True)
```
chr1_intra = extract_single_map(input_global_matrix=global_obs, 
tab_sep=False, 
chr_row='1', chr_col='1', 
bin_size=1000000,
data_type='observed',
save_output=True,
save_tab=True)

sum(sum(norm1==temp))

plot_map(input_global_matrix=global_obs,
tab_sep=False,
bin_size=1000000,
data_type='observed',
species='hg38',
my_colormap=['white', 'red'],
cutoff_type='perc',
cutoff=99,
max_color='#460000')


plot_chromosome_enrich_data(contact_matrix=enrich, 
a_chr='6', 
bin_size=1000000, 
full_matrix=True, 
species='hg38',
cutoff_max=4,
cutoff_min=-4,
plot_histogram=True)

plot_chromosome_data(temp, 
a_chr='6', 
bin_size=40000, 
full_matrix=False, 
start_coord=50000000, end_coord=54000000, 
species='hg38', 
data_type="normalized_fend", 
my_colormap=['white', 'red'], 
cutoff_type='percentile', cutoff=95, max_color='#460000', 
my_dpi=1000, 
plot_histogram=True)

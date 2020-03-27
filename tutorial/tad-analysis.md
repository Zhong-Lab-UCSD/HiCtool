# TAD analysis

This pipeline illustrates the procedure to calculate topologically associated domain (TAD) coordinates and visualize TADs and other features. Topological domain coordinates should be calculated on the normalized data (as we do here).

**Note!** Contact data must be at 40 kb bin size to perform TAD analysis!

## Table of Contents

1. [Performing TAD analysis](#1-performing-tad-analysis)
2. [Plotting TADs and more](#2-plotting-tads-and-more)
   - [2.1. Plotting TADs on the heatmap](#21-plotting-tads-on-the-heatmap)
   - [2.2. Plotting DI values and HMM states](#22-plotting-di-values-and-hmm-states)
   
## 1. Performing TAD analysis

TAD coordinates, as well as DI values and true DI values (HMM states), are calculated as described by [Dixon et al., (2012)](http://www.nature.com/nature/journal/v485/n7398/abs/nature11082.html).

**Note!** Contact data must be at 40 kb bin size to perform TAD analysis!

To perform TAD analysis (calculating DI, HMM states and topological domain coordinates) we use the function ``full_tad_analysis`` of [HiCtool_TAD_analysis.py](/scripts/HiCtool_TAD_analysis.py) and each single intra-chromosomal contact matrix as input files.

Here we take chromosome 6 as example. For this case, you may either have normalized the data using the Hi-Corrector approach (``./normalized_40000/chr6_chr6_40000.txt``) or have used the approach from Yaffe and Tanay (``./yaffe_tanay_40000/chr6_40000_normalized_fend.txt``). In this last case, remember to set ``--tab_sep 0`` below.
```unix
python2.7 ./HiCtool-master/scripts/HiCtool_TAD_analysis.py \
--action full_tad_analysis \
-i ./normalized_40000/chr6_chr6_40000.txt \
-c ./HiCtool-master/scripts/chromSizes/ \
-s hg38 \
--isGlobal 0 \
--tab_sep 1 \
--chr 6 \
--data_type normalized
```
where:

- ``--action``: Action to perform (here ``full_tad_analysis``).
- ``-i``: Input contact matrix file.
- ``-c``: Path to the folder ``chromSizes`` with trailing slash at the end ``/``.
- ``-s``: Species name.
- ``--isGlobal``: 1 if the input matrix is a global matrix, 0 otherwise.
- ``--tab_sep``: 1 if the input matrix is in a tab separated format, 0 in compressed format.
- ``--chr``: Chromosome or chromosomes to perform the TAD analysis in a list between square brackets.
- ``--data_type``: Data type to label your data, here ``normalized``.

This script will produce three output files inside the folder ``tad_analysis``:

- ``HiCtool_chr6_DI.txt`` which contains the DI values.
- ``HiCtool_chr6_hmm_states.txt`` which contains the HMM states extracted from the DI values.
- ``HiCtool_chr6_topological_domains.txt`` which contains topological domains coordinates in a tab separated format with two columns. Each row is a topological domains, first column is start coordinate, second column is end coordinate.

**Note!** The end coordinate of each domain is saved as the start position of the last bin (40 kb) belonging to each domain.

To calculate the **topological domain coordinates for multiple chromosomes** you may use an approach as the following:
```unix
chromosomes=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "X" "Y")

for i in "${chromosomes[@]}"; do
	python2.7 ./HiCtool-master/scripts/HiCtool_TAD_analysis.py \
	--action full_tad_analysis \
	-i "./normalized_40000/chr"$i"_chr"$i"_40000.txt" \
	-c ./HiCtool-master/scripts/chromSizes/ \
	-s hg38 \
	--isGlobal 0 \
	--tab_sep 1 \
	--chr $i \
	--data_type normalized
done
```

## 2. Plotting TADs and more

This section allows you to plot TADs over the heatmap and also DI values and HMM states.

### 2.1. Plotting TADs on the heatmap

To plot the topological domains on the heatmap, use the function ``plot_map`` of [HiCtool_global_map_analysis.py](/scripts/HiCtool_global_map_analysis.py) and pass the topological domain coordinates with the argument ``--topological_domains``. Here we plot the heatmap for chr6: 80,000,000-120,000,000. The color of the domains can be changed by using the parameter ``--domain_color``.

**Note!** To plot topological domains on the heatmap this should be with a bin size of 40kb or lower!

```unix
python2.7 ./HiCtool-master/scripts/HiCtool_global_map_analysis.py \
--action plot_map \
-i ./normalized_40000/chr6_chr6_40000.txt \
-c ./HiCtool-master/scripts/chromSizes/ \
-b 40000 \
-s hg38 \
--isGlobal 0 \
--tab_sep 1 \
--data_type normalized \
--chr_row 6 \
--chr_col 6 \
--my_colormap [white,red] \
--cutoff_type percentile \
--cutoff 95 \
--max_color "#460000" \
--chr_row_coord [80000000,120000000] \
--chr_col_coord [80000000,120000000] \
--topological_domains ./tad_analysis/HiCtool_chr6_topological_domains.txt
```
![](/figures/HiCtool_chr6_chr6_40kb_normalized_domains_80000000_120000000.png)

where:
- ``--action``: Action to perform (here ``plot_map``).
- ``-i``: Input contact matrix file.
- ``-c``: Path to the folder ``chromSizes`` with trailing slash at the end ``/``.
- ``-b``: The bin size (resolution) for the analysis.
- ``-s``: Species name.
- ``--isGlobal``: 1 if the input matrix is a global matrix, 0 otherwise.
- ``--tab_sep``: 1 if the input matrix is in a tab separated format, 0 if it is in compressed format.
- ``--data_type``: Data type to label your data, example: observed, normalized, etc.
- ``--chr_row``: Chromosome to plot the DI values of.
- ``--chr_col``: Same as ``--chr_row`` since this is only for intra-chromosomal maps.
- ``--my_colormap``: Colormap to be used to plot the data. You can choose among any colorbar at https://matplotlib.org/examples/color/colormaps_reference.html, or input a list of colors if you want a custom colorbar. Example: ``[white, red, black]``. Colors can be specified also HEX format. Default: ``[white,red]``.
- ``--cutoff_type``: To select a type of cutoff (``percentile`` or ``contact``) or plot the full range of the data (not declared). Default: ``percentile``.
- ``--cutoff``: To set a maximum cutoff on the number of contacts for the colorbar based on ``--cutoff_type``. Default: ``95``.
- ``--max_color``: To set the color of the bins with contact counts over ``--cutoff``. Default: ``#460000``.
- ``--chr_row_coord``:  List of two integers with start and end coordinates for the chromosome on the rows.
- ``--chr_col_coord``: List of two integers with start and end coordinates for the chromosome on the columns.
- ``--topological_domains``: Topological domain coordinates file to visualize domains on the heatmap.

Zoom on a smaller region (chr6: 50,000,000-54,000,000):
```unix
python2.7 ./HiCtool-master/scripts/HiCtool_global_map_analysis.py \
--action plot_map \
-i ./normalized_40000/chr6_chr6_40000.txt \
-c ./HiCtool-master/scripts/chromSizes/ \
-b 40000 \
-s hg38 \
--isGlobal 0 \
--tab_sep 1 \
--data_type normalized \
--chr_row 6 \
--chr_col 6 \
--my_colormap [white,red] \
--cutoff_type percentile \
--cutoff 95 \
--max_color "#460000" \
--chr_row_coord [50000000,54000000] \
--chr_col_coord [50000000,54000000] \
--topological_domains ./tad_analysis/HiCtool_chr6_topological_domains.txt
```
![](/figures/HiCtool_chr6_chr6_40kb_normalized_domains_50000000_54000000.png)


### 2.2. Plotting DI values and HMM states

**Directionality Index (DI)**

Regions at the periphery of the topological domains are highly biased in their interaction frequencies: the most upstream portion of a topological domain is highly biased towards interacting downstream, and the downstream portion of a topological domain is highly biased towards interacting upstream. To determine the directional bias at any given bin in the genome the Directionality Index (DI) is used, quantifying the degree of upstream or downstream bias of a given bin:

This is the formula used to calculate the DI:

![](/figures/DI_formula.png)

where:

- *A* is the the number of reads that map from a given 40 kb bin to the upstream 2 Mb.
- *B* is the the number of reads that map from a given 40 kb bin to the downstream 2 Mb.
- *E* is the expected number of contacts for each bin and it equals (A+B)/2.

To compute DI, we need the fend normalized contact data at a bin size of 40 kb. Since the bin size is 40 kb (here we used 40 kbpb which stands for 40 kb per bin), hence the detection region of upstream or downstream biases 2 Mb is converted to 50 bins (2 Mb / 40 kbpb = 50 bins).

To **plot the DI values** use the function ``plot_chromosome_DI`` of [HiCtool_TAD_analysis.py](/scripts/HiCtool_TAD_analysis.py) as following (in this case we plot DI values for chromosome 6, from 50 to 54 Mb):
```unix
python2.7 ./HiCtool-master/scripts/HiCtool_TAD_analysis.py \
--action plot_chromosome_DI \
-i ./tad_analysis/HiCtool_chr6_DI.txt \
-c ./HiCtool-master/scripts/chromSizes/ \
-s hg38 \
--chr 6 \
--full_chromosome 0 \
--coord [50000000,54000000]
```
where:
- ``--action``: Action to perform (here ``plot_chromosome_DI``).
- ``-i``: Input DI file.
- ``-c``: Path to the folder ``chromSizes`` with trailing slash at the end ``/``.
- ``-b``: The bin size (resolution) for the analysis.
- ``-s``: Species name.
- ``--chr``: Chromosome.
- ``--full_chromosome``: 1 to plot the DI values for the entire chromosome, 0 otherwise (must use ``--coord``).
- ``--coord``: List with start and end coordinates to plot the DI values.

![](/figures/HiCtool_chr6_DI.png)

**HMM states**

A Hidden Markov Model (HMM) based on the Directionality Index is used to identify biased states (true DI).

For true DI calculation, we consider the **emission sequence** as the observed DI values and the Transition Matrix, Emission Matrix and initial State Sequence as unknown. We have **three emissions** 1, 2, 0 corresponding to a positive (1), negative (2) or zero (0) value of the observed DI. In our analysis, we associate to the emission '0' all the absolute DI values under a threshold of 0.4. So, first we estimate the model parameters and then the most probable sequence of states using the Viterbi algorithm. 

To **plot the DI values and HMM states** use the function ``plot_chromosome_DI`` of [HiCtool_TAD_analysis.py](/scripts/HiCtool_TAD_analysis.py) and add the parameter ``--input_file_hmm`` to input the HMM states file as following:
```unix
python2.7 ./HiCtool-master/scripts/HiCtool_TAD_analysis.py \
--action plot_chromosome_DI \
-i ./tad_analysis/HiCtool_chr6_DI.txt \
-c ./HiCtool-master/scripts/chromSizes/ \
-s hg38 \
--chr 6 \
--full_chromosome 0 \
--coord [50000000,54000000] \
--input_file_hmm ./tad_analysis/HiCtool_chr6_hmm_states.txt
```

![](/figures/HiCtool_chr6_DI_HMM.png)

**TAD coordinates**

The true DI values allow to infer the locations of the topological domains in the genome. A domain is initiated at the beginning of a single downstream biased HMM state (red color in the above figure). The domain is continuous throughout any consecutive downstream biased state. The domain will then end when the last in a series of upstream biased states (green color in the above figure) is reached, with the domain ending at the end of the last HMM upstream biased state.

To calculate the topological domains coordinates, first we extract all the potential start and end coordinates according to the definition, and then we evaluate a list of conditions to take into account the possible presence of gaps between a series of positive or negative states values. The figure below shows a summary of the procedure:

![](/figures/HiCtool_topological_domains_flowchart.png)

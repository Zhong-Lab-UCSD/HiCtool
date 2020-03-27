# Data normalization with the matrix balancing approach of Hi-Corrector

This pipeline illustrates the procedure to generate and normalize a **global Hi-C contact map** (intra- and inter-chromosomal interactions) following the matrix balancing approach of [Hi-Corrector](http://www.nature.com/ng/journal/v43/n11/abs/ng.947.html). Several visualization functionalities are also implemented.

## Table of Contents

1. [Running HiFive functions](#1-running-hifive-functions)
2. [Generating the global observed contact matrix](#2-generating-the-global-observed-contact-matrix)
3. [Normalizing the global contact matrix](#3-normalizing-the-global-contact-matrix)
4. [Visualizing the data](#4-visualizing-the-data)
   - [4.1. Visualizing the global contact data](#41-visualizing-the-global-contact-data)
   - [4.2. Visualizing a single heatmap](#42-visualizing-a-single-heatmap)
   - [4.3. Visualizing maps on a side-by-side view](#43-visualizing-maps-on-a-side-by-side-view)  

## 1. Running HiFive functions

We resort to the HiFive package in order to generate the global observed contact matrix (intra- and inter-chromosomal contact maps all together to form a single, global contact matrix). HiFive allows to remove spurious ligation products, as well as PCR duplicates and non-informative reads before generating the contact matrix.

The Python script [HiCtool_hifive.py](/scripts/HiCtool_hifive.py) is used to run all the three steps needed in order to generate the observed contact data: 

- Creating the Fend object.
- Creating the HiCData object from a Fend object and mapped data in bam format. At this step spurious ligation products (paired-reads whose total distance to their respective restriction sites exceeds 500 bp) are removed. In addition, PCR duplicates are removed and reads with ends mapping to the same fragment and reads with ends mapping to adjacent fragments on opposite strands are also excluded, to consider the possibility of incomplete restriction enzyme digestion and fragment circularization.
- Creating the HiC project object, which stores the information about the observed contact data that we will use for the downstream analysis.

For more information about these functions, please see [HiFive’s API documentation](http://bxlab-hifive.readthedocs.org/en/latest/api.html). To run these steps execute the following command on the Unix console (update parameters properly):
```unix
python2.7 ./HiCtool-master/scripts/HiCtool_hifive.py \
-f restrictionsites.bed \
--b1 HiCfile_pair1.bam \
--b2 HiCfile_pair2.bam \
-e MboI \
-m Hi-Corrector
```
where:

- ``-f`` is the FEND file in bed format from preprocessing.
- ``--b1`` is the first bam file from preprocessing.
- ``--b2`` is the second bam file from preprocessing.
- ``-e`` is the restriction enzyme or enzymes names between square brackets (example ``[MboI,Hinfl]``).
- ``-m`` is the normalization model used (``Hi-Corrector`` in this case).

**The following output files are generated:**

- ``fend_object.hdf5``
- ``HiC_data_object.hdf5``
- ``HiC_project_object.hdf5`` to be used in the following section.


## 2. Generating the global observed contact matrix

This section allows to generate a **global square observed contact matrix** (24-by-24 chromosomes for hg38). The total number of bins of this big matrix will depend on the resolution of the data (for hg38 at 1Mb resolution is 3078x3078). The observed data contain the observed read count per each bin.

Especially at higher resolution, the generation of the global observed contact matrix may be computationally expensive and require long time. Therefore, we implemented a code to allow job parallelization (if your machine allows that). Each row of the contact matrix is computed in parallel, meaning all the contact matrices per each chromosome, and finally they are merged together to generate the global matrix. Each row of the matrix is saved in a temporary file, which is automatically deleted after the job is done.

To calculate and save the global observed contact matrix use the script [HiCtool_global_map_observed.sh](/scripts/HiCtool_global_map_observed.sh) and run this command:
```unix
# Make the bash script executable
chmod u+x ./HiCtool-master/scripts/HiCtool_global_map_observed.sh

# Run the script
./HiCtool-master/scripts/HiCtool_global_map_observed.sh \
-i HiC_project_object.hdf5 \
-h ./HiCtool-master/scripts/ \
-b 1000000 \
-s hg38 \
-p 24
```
where:

- ``-i``: project object file in ``.hdf5`` format obtained with ``HiCtool_hifive.py`` above.
- ``-h``: HiCtool scripts path with the trailing slash at the end ``/``.
- ``-b``: the bin size (resolution) for the analysis.
- ``-s``: species name.
- ``-p``: number of parallel threads to use. It has to be less or equal than the number of chromosomes of your species. If you do not have a multi-core machine, insert 1 here (note that it takes longer time at higher resolutions).

**The following output files are generated inside a folder ``observed_1000000`` (the folder name changes based on the resolution of the analysis):**

- ``HiCtool_observed_global_1000000.txt``, the global observed matrix saved in tab separated format. This matrix will be used in the next section to normalize the data, since Hi-Corrector required the input data in a tab separated format.
- ``chri_chrj_1000000.txt``, all the single observed contact matrices in tab separated format.

## 3. Normalizing the global contact matrix

Here we normalize the data using the sequential implementation from Hi-Corrector (ic_mes) which is "memory efficient and can run on any single computer with limited memory, even for Hi-C datasets of large size. It is designed to overcome the memory limit by loading a portion of data into the memory at each time." ([Li, Wenyuan, et al. "Hi-Corrector: a fast, scalable and memory-efficient package for normalizing large-scale Hi-C data." Bioinformatics 31.6 (2014): 960-962](https://academic.oup.com/bioinformatics/article/31/6/960/215261)).

The Hi-Corrector source code ([see here](https://github.com/Zhong-Lab-UCSD/HiCtool#installation)) is already inside ``/HiCtool-master/scripts/``. To normalize the data, run the following command:
```unix
# Make the bash script executable
chmod u+x ./HiCtool-master/scripts/HiCtool_normalize_global_matrix.sh

# Run the script
./HiCtool-master/scripts/HiCtool_normalize_global_matrix.sh \
-h ./HiCtool-master/scripts/ \
-i ./observed_1000000/HiCtool_observed_global_1000000.txt \
-b 1000000 \
-q 100 \
-m 32000 \
-s hg38
```
where:

- ``-h``: HiCtool scripts path with the trailing slash at the end ``/``.
- ``-i``: the observed global contact matrix in tab delimited format.
- ``-b``: the bin size (resolution) for the analysis.
- ``-q``: maximum number of iterations performed in the algorithm.
- ``-m``: the memory size in Megabytes (MB). The bigger memory you allocate for the normalization process, the faster it is. Even 100 Mb is fine for 1 Mb resolution map, suggested at least 16000 Mb (16 GB) for 40 kb resolution.
- ``-s``: species name.
- ``-u`` (optional): the row sum after normalization. The iterative correction algorithm can allow users to specify the row sum after the normalization, because this method is a matrix scaling approach that normalizes the matrix to be a doubly stochastic matrix (rows and columns sums equal to 1). Then we can multiple each element of this normalized matrix by the given value of this parameter, say 10.0 or 100.0 or whatever you choose. In such a way, the row sums of normalized matrix becomes this number (10.0 or 100.0 or whatever you choose). As **default**, we use a row sum value calculated as "the average number of contacts of the observed matrix multiplied by the number of rows" to make the normalized data counts "comparable" with the observed ones. **If you wish to compare several samples, declare the same value for this parameter for each sample.**

**The following output files are generated inside a folder ``normalized_1000000`` (the folder name changes based on the resolution of the analysis):**

- ``HiCtool_normalized_global_1000000.txt``, the normalized global matrix saved in tab separated format.
- ``chri_chrj_1000000.txt``, all the single normalized contact matrices in tab separated format.

Another folder ``output_ic_mes`` with 2 files inside is also generated, deriving from the ic_mes algorithm:

- ``output.log``: a log file.
- ``output.bias``: a bias file used by the software to normalize the data.

## 4. Visualizing the data

To plot the contact maps we use the function ``plot_map`` of [HiCtool_global_map_analysis.py](/scripts/HiCtool_global_map_analysis.py). 

**Note that at higher resolution, loading and plotting a global contact matrix will most likely produce an error given to memory limitation.** Therefore, for high resolution matrices (what high resolution means it depends on your machine specifications), you may only be able to visualize single contact matrices.

### 4.1. Visualizing the global contact data

You can visualize either the observed or the normalized data. Here we plot both the global maps at 1 Mb resolution as calculated above.
```unix
# Observed data
python2.7 ./HiCtool-master/scripts/HiCtool_global_map_analysis.py \
--action plot_map \
-i ./observed_1000000/HiCtool_observed_global_1000000.txt \
-c ./HiCtool-master/scripts/chromSizes/ \
-b 1000000 \
-s hg38 \
--isGlobal 1 \
--tab_sep 1 \
--data_type observed \
--my_colormap [white,red] \
--cutoff_type percentile \
--cutoff 99 \
--max_color "#460000"

```
![](/figures/HiCtool_1mb_observed.png)

```unix
# Normalized data
python2.7 ./HiCtool-master/scripts/HiCtool_global_map_analysis.py \
--action plot_map \
-i ./normalized_1000000/HiCtool_normalized_global_1000000.txt \
-c ./HiCtool-master/scripts/chromSizes/ \
-b 1000000 \
-s hg38 \
--isGlobal 1 \
--tab_sep 1 \
--data_type normalized \
--my_colormap [white,red] \
--cutoff_type percentile \
--cutoff 99 \
--max_color "#460000"
```
![](/figures/HiCtool_1mb_normalized.png)

where:
- ``--action``: action to perform (here ``plot_map``).
- ``-i``: Input global contact matrix file.
- ``-c``: Path to the folder ``chromSizes`` with trailing slash at the end ``/``.
- ``-b``: The bin size (resolution) for the analysis.
- ``-s``: Species name.
- ``--isGlobal``: Set to 1 if the input matrix is a global matrix, 0 otherwise.
- ``--tab_sep``: Set to 1 if the input matrix is in a tab separated format, 0 if it is in compressed format.
- ``--data_type``: Data type to label your data, example: observed, normalized, etc.
- ``--my_colormap``: Colormap to be used to plot the data. You can choose among any colorbar at https://matplotlib.org/examples/color/colormaps_reference.html, or input a list of colors if you want a custom colorbar. Example: ``[white, red, black]``. Colors can be specified also HEX format. Default: ``[white,red]``.
- ``--cutoff_type``: To select a type of cutoff (``percentile`` or ``contact``) or plot the full range of the data (not declared). Default: ``percentile``.
- ``--cutoff``: To set a maximum cutoff on the number of contacts for the colorbar based on ``--cutoff_type``. Default: ``95``.
- ``--max_color``: To set the color of the bins with contact counts over ``--cutoff``. Default: ``#460000``.

The resolution in DPI of the output PDF file can be changed using ``--my_dpi``. Default is 2000. Be aware that very high DPI levels could not be feasible due to memory limitations.

To emphasize the inter-chromosomal contacts in the global matrix you may use a lower cut-off setting ``--cutoff 95``:

![](/figures/HiCtool_1mb_normalized_95th.png)


### 4.2. Visualizing a single heatmap

- **Visualizing a single heatmap by passing as input the global contact matrix**

A single contact matrix can be plotted by passing the global contact matrix and as arguments the chromosome(s) in the rows (``--chr_row``) and in the columns (``--chr_col``) as a list between square brackets.

**Tip!** Loading the global map to extract the single matrices to plot may require time especially at higher resolutions. Therefore, this solution is useful when you want to plot several matrices at once using ``--chr_row`` and ``--chr_col`` as below. Note that at very high resolution loading the global matrix may result in a memory error, thus you should plot each heatmap passing a single contact matrix.

- **Visualizing a single heatmap by passing as input the single contact matrix**

A single heatmap can be plotted also using the single contact matrix as input. **At higher resolution, this is suggested to avoid any memory error.**

***

To plot the **intra-chromosomal heatmap** of chromosome 6 and **inter-chromosomal heatmap** (chr6-chr3) from the global contact matrix, run the following:
```unix
# Observed data
python2.7 ./HiCtool-master/scripts/HiCtool_global_map_analysis.py \
--action plot_map \
-i ./observed_1000000/HiCtool_observed_global_1000000.txt \
-c ./HiCtool-master/scripts/chromSizes/ \
-b 1000000 \
-s hg38 \
--isGlobal 1 \
--tab_sep 1 \
--data_type observed \
--chr_row [6,6] \
--chr_col [6,3] \
--my_colormap [white,red] \
--cutoff_type percentile \
--cutoff 99 \
--max_color "#460000"

# Normalized data
python2.7 ./HiCtool-master/scripts/HiCtool_global_map_analysis.py \
--action plot_map \
-i ./normalized_1000000/HiCtool_normalized_global_1000000.txt \
-c ./HiCtool-master/scripts/chromSizes/ \
-b 1000000 \
-s hg38 \
--isGlobal 1 \
--tab_sep 1 \
--data_type normalized \
--chr_row [6,6] \
--chr_col [6,3] \
--my_colormap [white,red] \
--cutoff_type percentile \
--cutoff 99 \
--max_color "#460000"
```
Observed (chr 6)           |  Normalized (chr 6)
:-------------------------:|:-------------------------:
![](/figures/HiCtool_chr6_chr6_1mb_observed.png)  |  ![](/figures/HiCtool_chr6_chr6_1mb_normalized.png)


Observed (chr6-chr3)            |  Normalized (chr6-chr3)
:-------------------------:|:-------------------------:
![](/figures/HiCtool_chr6_chr3_1mb_observed.png)  |  ![](/figures/HiCtool_chr6_chr3_1mb_normalized.png)

A histogram of the contact data distribution can be plotted by setting ``--plot_histogram 1``. This can be used only when single heatmaps are plotted.

To plot the **observed inter-chromosomal heatmap chr6-chr3 passing the single contact matrix as input** run the following (in this case you can only plot one matrix at the time):
```unix
# Observed data
python2.7 ./HiCtool-master/scripts/HiCtool_global_map_analysis.py \
--action plot_map \
-i ./observed_1000000/chr6_chr3_1000000.txt \
-c ./HiCtool-master/scripts/chromSizes/ \
-b 1000000 \
-s hg38 \
--isGlobal 0 \
--tab_sep 1 \
--data_type observed \
--chr_row 6 \
--chr_col 3 \
--my_colormap [white,red] \
--cutoff_type percentile \
--cutoff 99 \
--max_color "#460000"
```

In addition, only a **region of the heatmap** can be plotted by setting the parameters ``--chr_row_coord`` and ``--chr_col_coord``. These are lists with two integers indicating the start and end coordinate of the chromosome on the rows and on the columns respectively. If several single maps are inputed at once, these parameters can be lists of lists, each with coordinates corresponding to a single heatmap (see below).
```unix
python2.7 ./HiCtool-master/scripts/HiCtool_global_map_analysis.py \
--action plot_map \
-i ./normalized_1000000/HiCtool_normalized_global_1000000.txt \
-c ./HiCtool-master/scripts/chromSizes/ \
-b 1000000 \
-s hg38 \
--isGlobal 1 \
--tab_sep 1 \
--data_type normalized \
--chr_row [6,6] \
--chr_col [6,3] \
--chr_row_coord [[0,80000000],[0,50000000]] \
--chr_col_coord [[0,80000000],[0,80000000]] \
--my_colormap [white,red] \
--cutoff_type percentile \
--cutoff 99 \
--max_color "#460000"
```
Normalized (chr6) 0-80 Mb         |  Normalized (chr6-chr3) 0-50Mb; 0-80Mb
:-------------------------:|:-------------------------:
![](/figures/HiCtool_chr6_chr6_1mb_0-80mb_normalized.png)  |  ![](/figures/HiCtool_chr6_chr3_1mb_0-50_0-80mb_normalized.png)


### 4.3. Visualizing maps on a side-by-side view

The side-by-side maps visualization may be helpful to visually compare different datasets. For this part, we will use the sample used so far (B-lymphoblastoids GM12878: GSM1551550) and the HEK293T sample GSM1081530. The latter dataset has been processed following the same pipeline than the former.

To plot side-by-side contact maps for two or more samples, we use the function ``plot_side_by_side_map`` of [HiCtool_global_map_analysis.py](/scripts/HiCtool_global_map_analysis.py). Here we plot all the intra-chromosomal maps. Only global maps can be passed as input ``-i`` (either observed or normalized) as list between square brackets.
```unix
chromosomes=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y]

python2.7 ./HiCtool-master/scripts/HiCtool_global_map_analysis.py \
--action plot_side_by_side_map \
-i [./observed_1000000/HiCtool_observed_global_1000000.txt,./HEK293T/observed_1000000/HiCtool_observed_global_1000000.txt] \
-c ./HiCtool-master/scripts/chromSizes/ \
-b 1000000 \
-s hg38 \
--tab_sep 1 \
--data_type observed \
--chr_row $chromosomes \
--chr_col $chromosomes \
--samples [GM12878,HEK293T]
```
where:
- ``-i``: List of input global contact matrices between square brackets.
- ``chr_row``: List of chromosomes on the rows of the contact maps to be plotted between square brackets.
- ``chr_col``: List of chromosomes on the columns of the contact maps to be plotted between square brackets.
- ``--samples``: Sample labels between square brackets.

![](/figures/HiCtool_GM12878_HEK293T_1mb_observed.png)

# Data normalization with the matrix balancing approach of Hi-Corrector

This pipeline illustrates the procedure to normalize a **global Hi-C contact map** (intra- and inter-chromosomal interactions) following the matrix balancing approach of [Hi-Corrector](http://www.nature.com/ng/journal/v43/n11/abs/ng.947.html).

## Table of Contents

1. [Running HiFive functions](#1-running-hifive-functions)
   - [1.1. Creating the Fend object](#11-creating-the-fend-object)
   - [1.2. Creating the HiCData object](#12-creating-the-hicdata-object)
   - [1.3. Creating the HiC project object](#13-creating-the-hic-project-object)
2. [Generating the global observed contact matrix](#2-generating-the-global-observed-contact-matrix)
   - [2.1. Single-processor matrix generation](#21-single-processor-matrix-generation)
3. [Normalizing the global contact matrix](#3-normalizing-the-global-contact-matrix)
4. [Visualizing the data](#4-visualizing-the-data)
   - [4.1. Visualizing the global contact data](#41-visualizing-the-global-contact-data)
   - [4.2. Visualizing a single heatmap](#42-visualizing-a-single-heatmap)

## 1. Running HiFive functions

We resort to the HiFive package in order to generate the global observed contact matrix (intra- and inter-chromosomal contact maps all together to form a single, global contact matrix). HiFive allows to remove spurious ligation products, as well as PCR duplicates and non-informative reads (see below).

The function ``run_hifive`` of [HiCtool_hifive.py](/scripts/HiCtool_hifive.py) is used to run all the three steps needed in order to obtain the observed contact data, whose outputs are .hdf5 files. For more information about these functions, please see [HiFive’s API documentation](http://bxlab-hifive.readthedocs.org/en/latest/api.html). To run these steps execute the following code on the Python or iPython console:
```Python
execfile('HiCtool_hifive.py')
run_hifive(fend_file='restrictionsites_gc_map_valid.bed',
           bam_file1='HiCfile_pair1.bam', 
           bam_file2='HiCfile_pair2.bam',
           restriction_enzyme='MboI',
           model='Hi-Corrector')
```
If you computed your custom FEND file, you can use the ``restrictionsites.bed`` file, without need of adding GC content and mappability score.

More details about all the steps performed here are illustrated in the following steps 1.1.-1.3. If not interested, go to [section 2](#2-generating-the-global-observed-contact-matrix).

***

### 1.1. Creating the Fend object

A Fragment-end (Fend) object (``hdf5`` format) contains information about the fragments created by digestion of a genome by a specific restriction enzyme (RE). In our script, this information is supplied in the form of a BED-formatted file (``restrictionsites_gc_map_valid.bed``) containing information about the fragment ends like coordinates, GC content and mappability score (see [preprocessing, step 5](/tutorial/data-preprocessing.md#5-creating-the-fragment-end-fend-bed-file)). In this case, the information of GC content and mappability score is not used. If you computed your custom FEND file, you can use the ``restrictionsites.bed`` file, without need of adding GC content and mappability score.

To create a Fend object use the function ``hifive.fend``:
```Python
import hifive

fend = hifive.Fend('fend_object.hdf5', mode='w')
fend.load_fends('restrictionsites_gc_map_valid.bed', re_name='MboI', format='bed')
fend.save()
```
### 1.2. Creating the HiCData object

**HiC dataset** (``hdf5`` format) created from a Fend object and mapped data in bam format. The two bam files are passed to the function as elements of a list. The ``maxinsert`` parameter (``<int>``) is a cutoff for removing paired-reads whose total distance to their respective restriction sites exceeds this value. According to Yaffe and Tanay, we choose a value of 500 bp to remove spurious ligation products. In addition, when the HiCData object is created, **PCR duplicates** are removed and reads with ends mapping to the same fragment and reads with ends mapping to adjacent fragments on opposite strands are also excluded, to consider the possibility of incomplete restriction enzyme digestion and fragment circularization.

To create a HiCData object use the function ``hifive.HiCData``:
```Python
import hifive

data = hifive.HiCData('HiC_data_object.hdf5', mode='w')
data.load_data_from_bam('fend_object.hdf5',
                        ['HiCfile_pair1.bam', 'HiCfile_pair2.bam'],
                        maxinsert=500,
                        skip_duplicate_filtering=False)
data.save()
```
### 1.3. Creating the HiC project object

The **HiC project object** (``hdf5`` format) is a compressed format object which stores the information about the observed contact data that we will use for the downstream analysis.

To create a HiC project object use the function ``hifive.HiC``:
```Python
import hifive

hic = hifive.HiC('HiC_project_object.hdf5', 'w')
hic.load_data('HiC_data_object.hdf5')
hic.save()
```
where ``HiC_project_object.hdf5`` specifies the location to save the HiC object to and ``HiC_data_object.hdf5`` is the data object.

## 2. Generating the global observed contact matrix

This section will allow you to generate a **global square observed contact matrix** (24-by-24 chromosomes for hg38). The total number of bins of this big matrix will depend on the resolution of the data, and it can be estimated as the entire genome length over the resolution (for hg38 at 1Mb resolution is 3078x3078). The observed data contain the observed read count per each bin.

Especially at higher resolution, the generation of the global observed contact matrix may be computationally expensive and require long time. Therefore, we implemented a code to allow job parallelization. Each row of the contact matrix is computed in parallel, meaning all the contact matrices per each chromosome, and finally they are merged together to generate the global matrix. Each row of the matrix is saved in a temporary file, which is automatically deleted after the job is done.

If you do not have a multi-processor machine available, go to [section 2.1.](#21-single-processor-matrix-generation) to run the code using a single processor.

To calculate and save the global observed contact matrix in parallel use the script [HiCtool_full_map_parallel.py](/scripts/HiCtool_full_map_parallel.py). **Open the script, update the parameters on the top and save.** Then, just execute the script in the Python or iPython console:
```Python
execfile('HiCtool_full_map_parallel.py')
```
The global matrix will be saved as default using a compressed format but **to normalize the data using Hi-Corrector a tab separated format is required**. Therefore, also the global contact matrix in tab separated format will be saved, and this will be the input file to the normalization algorithm of Hi-Corrector. Besides the global contact matrix file, another file named **info.txt** will be saved. This contains information that are required to be inserted as Hi-Corrector input as well.

After having generated the global observed contact matrix, it is possible to extract a single contact matrix (either intra- or inter-chromosomal) using the function ``extract_single_map`` of [HiCtool_full_map.py](/scripts/HiCtool_full_map.py) as following:
```Python
execfile('HiCtool_full_map.py')

chr1_intra = extract_single_map(input_global_matrix=my_global_matrix, 
                                tab_sep=False, 
                                chr_row='1', chr_col='1', 
                                bin_size=1000000,
                                data_type='observed',
                                save_output=True,
                                save_tab=True)

chr1_2_inter = extract_single_map(input_global_matrix=my_global_matrix, 
                                  tab_sep=False, 
                                  chr_row='1', chr_col='2', 
                                  bin_size=1000000,
                                  data_type='observed',
                                  save_output=True,
                                  save_tab=True)
```

### 2.1. Single-processor matrix generation

To generate the global observed contact matrix with a single processor, open the Python or iPython console and use the function ``compute_matrix_data_full_observed`` of [HiCtool_full_map.py](/scripts/HiCtool_full_map.py):
```Python
execfile('HiCtool_full_map.py')
my_global_matrix = compute_matrix_data_full_observed(input_file='HiC_project_object.hdf5',
                                                     bin_size=1000000, 
                                                     species='hg38', 
                                                     save_each_matrix=False, 
                                                     save_tab=True)
```
It is possible also to use a custom species (see API documentation). Note that the global matrix will be saved as default using a compressed format but **to normalize the data using Hi-Corrector a tab separated format is required**. The parameter ``save_tab`` set as ``True`` will save the global contact matrix in tab separated format, and this will be the input file to the normalization algorithm of Hi-Corrector. Besides the global contact matrix file, another file named **info.txt** will be saved. This contains information that are required to be inserted as Hi-Corrector input as well.

## 3. Normalizing the global contact matrix

Here we normalize the data using the sequential implementation from Hi-Corrector (ic_mes) which is "memory efficient and can run on any single computer with limited memory, even for Hi-C datasets of large size. It is designed to overcome the memory limit by loading a portion of data into the memory at each time." ([Li, Wenyuan, et al. "Hi-Corrector: a fast, scalable and memory-efficient package for normalizing large-scale Hi-C data." Bioinformatics 31.6 (2014): 960-962](https://academic.oup.com/bioinformatics/article/31/6/960/215261)).

First, you should have downloaded Hi-Corrector source code ([see here](https://github.com/Zhong-Lab-UCSD/HiCtool#installation)). Then to normalize the data, download the bash script [run_ic_mes.sh](/scripts/run_ic_mes.sh), copy it inside your working directory, go into your working directory and run the following unix commands:
```unix
chmod u+x run_ic_mes.sh
./run_ic_mes.sh -q 100 \
                -m "100" \
                -r 3078 \
                -s 17237 \
                -h "/path_to_Hi-Corrector/Hi-Corrector1.2" \
                -i "HiCtool_1mb_matrix_global_observed_tab.txt"
```
Where the options are:

- ``-q`` int: maximum number of iterations performed in the algorithm.
- ``-m`` str: the memory size. Its unit is Megabytes (MB).
- ``-r`` int: the number of rows or columns of the input chromatin contact frequency matrix to be normalized (provided in the **info.txt** file generated in [section 2](#2-generating-the-global-observed-contact-matrix)).
- ``-s`` int: the row sum after normalization. The iterative correction algorithm can allow users to specify the row sum after the normalization, because this method is a matrix scaling approach that normalizes the matrix to be a doubly stochastic matrix (rows and columns sums equal to 1). Then we can multiple each element of this normalized matrix by the given value of this parameter, say 10.0 or 100.0 or whatever you choose. In such a way, the row sums of normalized matrix becomes this number (10.0 or 100.0 or whatever you choose). In **info.txt** we provide a row sum value that you could use calculated as "the average number of contacts of the observed matrix multiplied by the number of rows" to make the normalized data counts "comparable" with the observed ones. The choice is arbitrary.
- ``-h`` str: the path to the Hi-Corrector source code that you have downloaded.
- ``-i`` str: the observed global contact matrix in tab delimited format.

This command creates a **folder named ``output_ic_mes``** with 3 files inside:

- ``output.log``: a log file
- ``output.bias``: a bias file used by the software to normalize the data
- ``output_normalized.txt``: the **global normalized contact matrix** in tab separated format

After having normalized the data, it is possible to extract a single normalized contact matrix (either intra- or inter-chromosomal) using the function ``extract_single_map`` of [HiCtool_full_map.py](/scripts/HiCtool_full_map.py) as following:
```Python
execfile('HiCtool_full_map.py')

chr1_intra_norm = extract_single_map(input_global_matrix="output_ic_mes/output_normalized.txt", 
                                     tab_sep=True, 
                                     chr_row='1', chr_col='1', 
                                     bin_size=1000000,
                                     data_type='normalized',
                                     save_output=True,
                                     save_tab=True)  

chr1_2_inter_norm = extract_single_map(input_global_matrix="output_ic_mes/output_normalized.txt", 
                                       tab_sep=True, 
                                       chr_row='1', chr_col='2', 
                                       bin_size=1000000,
                                       data_type='normalized',
                                       save_output=True,
                                       save_tab=True)
```

## 4. Visualizing the data

To plot the contact maps use the function ``plot_map`` of [HiCtool_full_map.py](/scripts/HiCtool_full_map.py).
```Python
execfile('HiCtool_full_map.py')
```

### 4.1. Visualizing the global contact data

You can visualize either the observed or the normalized data. Here we plot both the global maps at 1 Mb resolution as calculated above.
```Python
# Observed data
plot_map(input_global_matrix='HiCtool_1mb_matrix_global_observed.txt',
         tab_sep=False,
         bin_size=1000000,
         data_type='observed',
         species='hg38',
         my_colormap=['white', 'red'],
         cutoff_type='perc',
         cutoff=99,
         max_color='#460000')
```
![](/figures/HiCtool_1mb_observed.png)

```Python
# Normalized data
plot_map(input_global_matrix='output_ic_mes/output_normalized.txt',
         tab_sep=True,
         bin_size=1000000,
         data_type='normalized',
         species='hg38',
         my_colormap=['white', 'red'],
         cutoff_type='perc',
         cutoff=99,
         max_color='#460000')
```
![](/figures/HiCtool_1mb_normalized.png)

### 4.2. Visualizing a single heatmap

A single contact matrix can be plotted by passing as argument the chromosome in the rows (``chr_row``) and in the columns (``chr_col``). 

To plot the **intra-chromosomal heatmap** of chromosome 6, run the following:
```Python
# Observed contact heatmap
plot_map(input_global_matrix='HiCtool_1mb_matrix_global_observed.txt',
         tab_sep=False,
         chr_row='6', chr_col='6', 
         bin_size=1000000, 
         data_type="observed",
         species='hg38',
         my_colormap=['white', 'red'],
         cutoff_type='perc',
         cutoff=99,
         max_color='#460000')

# Normalized contact heatmap
plot_map(input_global_matrix='output_ic_mes/output_normalized.txt',
         tab_sep=True,
         chr_row='6', chr_col='6', 
         bin_size=1000000, 
         data_type="normalized",
         species='hg38',
         my_colormap=['white', 'red'],
         cutoff_type='perc',
         cutoff=99,
         max_color='#460000')
```
Observed (chr 6)           |  Normalized (chr 6)
:-------------------------:|:-------------------------:
![](/figures/HiCtool_chr6_chr6_1mb_170x170_observed.png)  |  ![](/figures/HiCtool_chr6_chr6_1mb_170x170_normalized.png)

An **inter-chromosomal heatmap** can be also plotted (chr6-chr3) by setting the parameters ``chr_row`` and ``chr_col`` (we plot also the histogram of the contact distribution):
```Python
# Observed contact heatmap
plot_map(input_global_matrix='HiCtool_1mb_matrix_global_observed.txt',
         tab_sep=False,
         chr_row='6', chr_col='3', 
         bin_size=1000000, 
         data_type="observed",
         species='hg38',
         my_colormap=['white', 'red'],
         cutoff_type='perc',
         cutoff=99,
         max_color='#460000',
         plot_histogram=True)

# Normalized contact heatmap
plot_map(input_global_matrix='output_ic_mes/output_normalized.txt',
         tab_sep=True,
         chr_row='6', chr_col='3', 
         bin_size=1000000, 
         data_type="normalized",
         species='hg38',
         my_colormap=['white', 'red'],
         cutoff_type='perc',
         cutoff=99,
         max_color='#460000',
         plot_histogram=True)
```
Observed (chr6-chr3)            |  Normalized (chr6-chr3)
:-------------------------:|:-------------------------:
![](/figures/HiCtool_chr6_chr3_1mb_170x198_observed.png)  |  ![](/figures/HiCtool_chr6_chr3_1mb_170x198_normalized.png)
![](/figures/HiCtool_chr6_chr3_1mb_170x198_observed_histogram.png)  |  ![](/figures/HiCtool_chr6_chr3_1mb_170x198_normalized_histogram.png)

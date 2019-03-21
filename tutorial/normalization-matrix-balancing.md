# Data normalization with the matrix balancing approach of Hi-Corrector

This pipeline illustrates the procedure to normalize a **global Hi-C contact map** (intra- and inter-chromosomal interactions) following the matrix balancing approach of [Hi-Corrector](http://www.nature.com/ng/journal/v43/n11/abs/ng.947.html).

## Table of Contents

1. [Running HiFive functions](#1-running-hifive-functions)
   - [1.1. Creating the Fend object](#11-creating-the-fend-object)
   - [1.2. Creating the HiCData object](#12-creating-the-hicdata-object)
   - [1.3. Creating the HiC project object](#13-creating-the-hic-project-object)
2. [Generating the global observed contact matrix](#2-generating-the-global-observed-contact-matrix)
3. [Normalizing the global contact matrix](#3-normalizing-the-global-contact-matrix)
   - [3.1. Multi-processing normalization](#31-multi-processing-normalization)
4. [Visualizing the data](#4-visualizing-the-data)
   - [3.1. Visualizing the global contact data](#41-visualizing-the-global-contact-data)
   - [3.2. Visualizing a single heatmap](#42-visualizing-a-single-heatmap)

## 1. Running HiFive functions

We resort to the HiFive package in order to generate the global observed contact matrix (intra- and inter-chromosomal contact maps all together to form a single, global contact matrix). HiFive allows to remove spurious ligation products, as well as PCR duplicates and non-informative reads (see below).

The script [HiCtool_hifive.py](/scripts/HiCtool_hifive.py) is used to run all the three steps needed in order to obtain the observed contact data, whose outputs are .hdf5 files. For more information about these functions, please see [HiFive’s API documentation](http://bxlab-hifive.readthedocs.org/en/latest/api.html). To run these steps execute the following code on the Python or iPython console:
```Python
execfile('HiCtool_hifive.py')
run_hifive('restrictionsites_gc_map_valid.bed',
           'HiCfile_pair1.bam', 
           'HiCfile_pair2.bam',
           'MboI',
           'Hi-Corrector')
```
More details about all the steps performed here are illustrated in the following steps 1.1.-1.3. If not interested, go to [section 2](#2-generating-the-global-observed-contact-matrix).

### 1.1. Creating the Fend object

A Fragment-end (Fend) object (hdf5 format) contains information about the fragments created by digestion of a genome by a specific restriction enzyme (RE). In our script, this information is supplied in the form of a BED-formatted file (```restrictionsites_gc_map_valid.bed```) containing information about the fragment ends like coordinates, GC content and mappability score (see [preprocessing, step 5](/tutorial/data-preprocessing.md#5-creating-the-fragment-end-fend-bed-file)). In this case, the information of GC content and mappability score is not used.

To create a Fend object use the function ```hifive.fend```:
```Python
import hifive

fend = hifive.Fend('fend_object.hdf5', mode='w')
fend.load_fends('restrictionsites_gc_map_valid.bed', re_name='MboI', format='bed')
fend.save()
```
### 1.2. Creating the HiCData object

**HiC dataset** (hdf5 format) created from a Fend object and mapped data in bam format. The two bam files are passed to the function as elements of a list. The ‘maxinsert’ parameter (int.) is a cutoff for removing paired-reads whose total distance to their respective restriction sites exceeds this value. According to Yaffe and Tanay, we choose a value of 500 bp to remove spurious ligation products. In addition, when the HiCData object is created, **PCR duplicates** are removed and reads with ends mapping to the same fragment and reads with ends mapping to adjacent fragments on opposite strands are also excluded, to consider the possibility of incomplete restriction enzyme digestion and fragment circularization.

To create a HiCData object use the function ```hifive.HiCData```:
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

The **HiC project object** (hdf5 format) is a compressed format object which stores the information about the observed contact data that we will use for the downstream analysis.

To create a HiC project object use the function ```hifive.HiC```:
```Python
import hifive

hic = hifive.HiC('HiC_project_object.hdf5', 'w')
hic.load_data('HiC_data_object.hdf5')
hic.save()
```
where ```HiC_project_object.hdf5``` specifies the location to save the HiC object to and ```HiC_data_object.hdf5``` is the data object.

## 2. Generating the global observed contact matrix

This section will allow you to generate a **global square contact matrix** (24-by-24 chromosomes for hg38). The total number of bins of this big matrix will depend on the resolution of the data, and it can be estimated as the entire genome length ratio the resolution (for hg38 at 1Mb resolution is 3078x3078). The observed data contain the observed read count per each bin.

To generate the global observed contact matrix, open the Python or iPython console and use the function ```compute_matrix_data_full_observed``` of [HiCtool_full_map.py](/scripts/HiCtool_full_map.py):
```Python
execfile('HiCtool_full_map.py')
my_global_matrix = compute_matrix_data_full_observed(input_file='HiC_project_object.hdf5',
                                                                bin_size=40000, 
                                                                species='hg38', 
                                                                save_each_matrix=False, 
                                                                save_tab=True)
```
The function contains easy parameters to be set. It is possible also to use a custom species (see API documentation). Note that the global matrix will be saved as default using a compressed format but to normalize the data using Hi-Corrector a tab separated format is required. The parameter ```save_tab```, if set as ```True```, will save the global contact matrix in tab separated format, and this will be the input file to the normalization algorithm of Hi-Corrector. Besides the global contact matrix file, another file named **info.txt** will be saved. This contains information that are required to be inserted as Hi-Corrector input as well.

After having generated the global contact matrix, it is possible to extract a single contact matrix (either intra or interchromosomal) using the function ``extract_single_map`` as following::

>>> chr1_intra = extract_single_map(input_global_matrix=my_global_matrix, tab_sep=False, chr_row='1', chr_col='1', bins_size=40000)

>>> chr1_2_inter = extract_single_map(input_global_matrix=my_global_matrix, tab_sep=False, chr_row='1', chr_col='2', bins_size=40000)

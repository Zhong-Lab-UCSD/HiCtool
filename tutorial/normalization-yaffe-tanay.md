# Data normalization with explicit-factor correction model of Yaffe and Tanay

This pipeline illustrates the procedure to normalize and visualize Hi-C **intra-chromosomal contact data only** following the explicit-factor model of [Yaffe and Tanay](http://www.nature.com/ng/journal/v43/n11/abs/ng.947.html). Default species that you can use are human *hg38* and mouse *mm10*. For other species use the specific function arguments (see [API documentation](https://sysbio.ucsd.edu/public/rcalandrelli/HiCtool_API_documentation.pdf)).

For more information about the Python functions used here check the [API documentation](https://sysbio.ucsd.edu/public/rcalandrelli/HiCtool_API_documentation.pdf).

## Table of contents

1. [Running HiFive functions](#1-running-hifive-functions)
   - [1.1. Creating the Fend object](#11-creating-the-fend-object)
   - [1.2. Creating the HiCData object](#12-creating-the-hicdata-object)
   - [1.3. Creating the HiC project object](#13-creating-the-hic-project-object)
   - [1.4. Filtering HiC fends](#14-filtering-hic-fends)
   - [1.5. Estimating the HiC distance function](#15-estimating-the-hic-distance-function)
   - [1.6. Learning the correction model](#16-learning-the-correction-model)
2. [Normalizing the data](#2-normalizing-the-data)
   - [2.1. Normalized fend data](#21-normalized-fend-data)
   - [2.2. Normalized enrichment data](#22-normalized-enrichment-data)
   - [2.3. Multi-processing normalization](#23-multi-processing-normalization)
3. [Visualizing the data](#3-visualizing-the-data)
   - [3.1. Visualizing the contact data](#31-visualizing-the-contact-data)
   - [3.2. Visualizing the enrichment data](#32-visualizing-the-enrichment-data)

## 1. Running HiFive functions

The function ``run_hifive`` of [HiCtool_hifive.py](/scripts/HiCtool_hifive.py) is used to run all the six steps needed in order to normalize the data, whose outputs are ``hdf5`` files. For more information about these functions, please see [HiFive’s API documentation](http://bxlab-hifive.readthedocs.org/en/latest/api.html). To run these steps execute the following code on the Python or iPython console:
```Python
execfile('HiCtool_hifive.py')
run_hifive(fend_file='restrictionsites_gc_map_valid.bed',
           bam_file_1='HiCfile_pair1.bam', bam_file_2='HiCfile_pair2.bam',
           restriction_enzyme='MboI',
           model='Yaffe-Tanay')
```
More details about all the steps performed here are illustrated in the following steps 1.1.-1.6. If not interested, jump to [section 2](#2-normalizing-the-data).

***

### 1.1. Creating the Fend object

A Fragment-end (Fend) object in HiFive (``hdf5`` format) contains information about the fragments created by digestion of a genome by a specific restriction enzyme (RE). In our script, this information is supplied in the form of a BED-formatted file (``restrictionsites_gc_map_valid.bed``) containing information about the fragment ends like coordinates, GC content and mappability score (see [preprocessing, step 5](/tutorial/data-preprocessing.md#5-creating-the-fragment-end-fend-bed-file)).

To create a Fend object use the function ``hifive.fend``:
```Python
import hifive

fend = hifive.Fend('fend_object.hdf5', mode='w')
fend.load_fends('restrictionsites_gc_map_valid.bed', re_name='MboI', format='bed')
fend.save()
```
### 1.2. Creating the HiCData object

**HiC dataset** (``hdf5`` format) created from a Fend object and mapped data in bam format. The two BAM files are passed to the function as elements of a list. The ``maxinsert`` parameter (``<int>``) is a cutoff for removing paired-reads whose total distance to their respective restriction sites exceeds this value. According to Yaffe and Tanay, we choose a value of 500 bp to remove spurious ligation products. In addition, when the HiCData object is created, **PCR duplicates** are removed and reads with ends mapping to the same fragment and reads with ends mapping to adjacent fragments on opposite strands are also excluded, to consider the possibility of incomplete restriction enzyme digestion and fragment circularization.

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

The **HiC project object** (``hdf5`` format) links the HiCData object with information about which fends to include in the analysis, model parameters and learned model values. This is the standard way of working with Hi-C data in HiFive and this object will be used for learning the correction model and downstream analysis.

To create a HiC project object use the function ``hifive.HiC``:
```Python
import hifive

hic = hifive.HiC('HiC_project_object.hdf5', 'w')
hic.load_data('HiC_data_object.hdf5')
hic.save()
```
where ``HiC_project_object.hdf5`` specifies the location to save the HiC object to and ``HiC_data_object.hdf5`` is the data object.

### 1.4. Filtering HiC fends

At this step we filter out fragments that do not have at least one interaction before learning correction parameters. Fragment ends within a distance of 500 kb are filtered out in step 1.6. This will allow to normalize data without confounding technical biases with features associated to biological-relevant structures (Yaffe and Tanay).

To filter out fends use the function ``hifive.filter_fends``:
```Python
import hifive

hic = hifive.HiC('HiC_project_object.hdf5')
hic.filter_fends(mininteractions=1, mindistance=0, maxdistance=0)
hic.save()
```

### 1.5. Estimating the HiC distance function

Estimation of the **distance-dependence relationship** from the data prior to normalization, in order to avoid biases that may result due to restriction site distribution characteristics or the influence of distance/signal relationship.

Restriction sites over the genome are unevenly distributed and this results in a large set of distances between fragments and their neighbors. Since the interaction frequency is strongly inversely-related to inter-fragment distance, this means that fragments surrounded by shorter ones will show higher nearby interactions than those with longer adjacent fragments, due to the uneven distribution of the restriction sites position.

To estimate the HiC distance function use ``hifive.find_distance_parameters``:
```Python
import hifive

hic = hifive.HiC('HiC_project_object.hdf5')
hic.find_distance_parameters(numbins=90, minsize=200, maxsize=0)
hic.save()
```
### 1.6. Learning the correction model

Algorithm to learn the correction model for Hi-C data. For the normalization, we take into account of fragments length, inter-fragment distance, GC content and mappability score biases, according to the information included in the Fend object. We also consider a minimum distance of 500 kb between fragments to take into account of the effect of biological biases (TSSs and CTCF bound sites) while learning the correction parameters.

To normalize the data using the binning algorithm ([Yaffe E. and Tanay A., 2011](http://www.ncbi.nlm.nih.gov/pubmed/22001755)) use ``hic.find_binning_fend_corrections``:
```Python
import hifive

hic = hifive.HiC('HiC_project_object.hdf5')
hic.find_binning_fend_corrections(max_iterations=1000,
                                  mindistance=500000,
                                  maxdistance=0,
                                  num_bins=[20,20,20,20],
                                  model=['len','gc','mappability','distance'],
                                  parameters=['even','even','even','even'],
                                  usereads='cis',
                                  learning_threshold=1.0)
hic.save('HiC_norm_binning.hdf5')
```

## 2. Normalizing the data

For the normalization, observed data and correction parameters to remove biases to obtain the corrected read counts are required. Therefore, the observed contact matrix and the fend expected contact matrix are calculated. In addition, the enrichment expected contact matrix is calculated to compute the observed over expected enrichment values, considering also the distance between fends.

For each chromosome, the following five matrices are computed at a bin size of 40 kb (the bin size can be changed with a function parameter). Every contact matrix is AUTOMATICALLY saved in txt format using the function ``save_matrix``.

Data are compressed in a format to reduce storage occupation and improving saving and loading time ([see here for more details](/tutorial/HiCtool_compressed_format.md)).

- The **observed data** contain the observed reads count for each bin.
- The **fend expected data** contain the learned correction value to remove biases related to fends for each bin.
- The **enrichment expected data** contain the expected reads count for each bin, considering the linear distance between read pairs and the learned correction parameters.
- The **normalized fend data** contain the corrected read count for each bin.
- The **normalized enrichment data** ("observed over expected" matrix) contain the enrichment value (O/E) for each bin.

First execute the script [HiCtool_yaffe_tanay.py](/scripts/HiCtool_yaffe_tanay.py) in the Python or iPython console:
```Python
execfile('HiCtool_yaffe_tanay.py')
```
### 2.1. Normalized fend data

To calculate and save the **normalized intra-chromosomal contact matrix** for a chromosome ``a_chr``, use the function ``normalize_chromosome_fend_data``:
```Python
fend_normalized_chr6 = normalize_chromosome_fend_data(a_chr='6', bin_size=40000, 
                                                      input_file='HiC_norm_binning.hdf5', 
                                                      species='hg38',
                                                      save_obs=True, save_expect=False)
```
Data are compressed in a format to reduce storage occupation and improving saving and loading time ([see here for more details](/tutorial/HiCtool_compressed_format.md)). To load a previously generated contact matrix use the function ```load_matrix```:
```Python
my_contact_matrix = load_matrix('my_contact_matrix.txt')
```
where ``'my_contact_matrix.txt'`` is a contact matrix file saved using ``normalize_chromosome_fend_data`` .

### 2.2. Normalized enrichment data

To calculate and save the **"observed/expected" intra-chromosomal contact matrix** for a chromosome ``a_chr`` use the function ``normalize_chromosome_enrich_data`` (see API Documentation):
```Python
enrich_normalized_chr6 = normalize_chromosome_enrich_data(a_chr='6', bin_size=40000, 
                                                          input_file='HiC_norm_binning.hdf5', 
                                                          species='hg38',
                                                          save_obs=True, save_expect=False)
```
**Note!**
If you need only the normalized contact matrices, there is no need to calculate also the enrichment data. If you do not need the expected data, do not save it since they are the biggest files and the process may take time.

Data are compressed in a format to reduce storage occupation and improving saving and loading time ([see here for more details](/tutorial/HiCtool_compressed_format.md)). To load a previously generated contact matrix use the function ```load_matrix```:
```Python
my_contact_matrix = load_matrix('my_contact_matrix.txt')
```
where ``'my_contact_matrix.txt'`` is a contact matrix file saved using ``normalize_chromosome_enrich_data``.

### 2.3. Multi-processing normalization

To calculate and save the normalized contact matrices in parallel for multiple chromosomes, use the script [HiCtool_normalize_fend_parallel.py](/scripts/HiCtool_normalize_fend_parallel.py). **Open the script, update the parameters on the top and save.** Then, just execute the script:
```Python
execfile('HiCtool_normalize_fend_parallel.py')
```
To calculate and save the "observed/expected" contact matrices in parallel use the script [HiCtool_normalize_enrich_parallel.py](/scripts/HiCtool_normalize_enrich_parallel.py). **Open the script, update the parameters on the top and save.** Then, just execute the script:
```Python
execfile('HiCtool_normalize_enrich_parallel.py')
```
## 3. Visualizing the data

To plot the contact maps, first execute the script [HiCtool_yaffe_tanay.py](/tutorial/HiCtool_yaffe_tanay.py) in the Python or iPython console:
```Python
execfile('HiCtool_yaffe_tanay.py')
```
### 3.1. Visualizing the contact data

This part is to plot heatmaps and histograms of the contact data. 

To plot and save the heatmap and histogram use the function ```plot_chromosome_data```:
```Python
plot_chromosome_data('HiCtool_chr6_40kb_normalized_fend.txt', 
                     a_chr='6', bin_size=40000, full_matrix=False, 
                     start_coord=50000000, end_coord=54000000, 
                     species='hg38', 
                     data_type="normalized_fend", 
                     my_colormap=['white', 'red'], 
                     cutoff_type='percentile', cutoff=95, max_color='#460000', 
                     my_dpi=1000, 
                     plot_histogram=True)
```
Instead of ``'HiCtool_chr6_40kb_normalized_fend.txt'``, the object containing the contact matrix calculated above ``fend_normalized_chr6`` can be passed as well.

This function can be used also to plot observed data, expected fend and enrichment data by simply passing a different input matrix as first parameter.

Heatmap             |  Histogram
:-------------------------:|:-------------------------:
![](/figures/HiCtool_chr6_40kb_normalized_fend.png)  |  ![](/figures/HiCtool_chr6_40kb_normalized_fend_histogram.png)


**Additional example of the contact matrix for chromosome 6 at 1 Mb resolution**

In order to change the heatmap resolution, first data have to be normalized at the desired resolution set with the parameter ``bin_size`` of ``normalize_chromosome_fend_data`` ([see section 2.1.](#21-normalized-fend-data)):
```Python
fend_normalized_chr6 = normalize_chromosome_fend_data(a_chr='6', 
                                                      bin_size=1000000, 
                                                      input_file='HiC_norm_binning.hdf5', 
                                                      species='hg38',
                                                      save_obs=True, 
                                                      save_expect=False)
```
Then, we plot the entire heatmap (we also change here the color map to white and blue):
```Python
plot_chromosome_data(fend_normalized_chr6, 
                     a_chr='6', 
                     bin_size=1000000, 
                     full_matrix=True, 
                     species='hg38', 
                     data_type="normalized_fend", 
                     my_colormap=['white', 'blue'], 
                     cutoff_type='percentile', cutoff=95, max_color='#460000', 
                     my_dpi=1000, 
                     plot_histogram=False)
```
![](/figures/HiCtool_chr6_1mb_normalized_fend.png)

### 3.2. Visualizing the enrichment data

This part is to plot the heatmap and histogram for the enrichment normalized data ("observed over expected"). The **log2 of the data** is plotted to quantify the positive enrichment (red) and the negative enrichment (blue). Loci (pixels) equal to zero before performing the log2 (deriving from zero observed contacts) are shown in gray. Loci (pixels) where enrichment expected contact was zero before performing the ratio (observed / expected) are shown in black.

To plot and save the heatmap and histogram use the function ```plot_chromosome_enrich_data```:
```Python
plot_chromosome_enrich_data('HiCtool_chr6_40kb_normalized_enrich.txt', 
                            a_chr='6', 
                            bin_size=40000, 
                            full_matrix=False, 
                            start_coord=50000000, end_coord=54000000, 
                            species='hg38', 
                            my_dpi=1000, 
                            plot_histogram=True)
```
Instead of ``'HiCtool_chr6_40kb_normalized_enrich.txt'``, the object containing the contact matrix calculated above ``enrich_normalized_chr6`` can be passed as well.

Heatmap             |  Histogram
:-------------------------:|:-------------------------:
![](/figures/HiCtool_chr6_40kb_normalized_enrich.png)  |  ![](/figures/HiCtool_chr6_40kb_normalized_enrich_histogram.png)

**Additional example of the enrichment contact matrix for chromosome 6 at 1 Mb resolution**

In order to change the heatmap resolution, first data have to be calculated at the desired resolution set with the parameter ``bin_size`` of ``normalize_chromosome_enrich_data`` ([see section 2.2.](#22-normalized-enrichment-data)):
```Python
enrich_normalized_chr6 = normalize_chromosome_enrich_data(a_chr='6', 
                                                          bin_size=1000000, 
                                                          input_file='HiC_norm_binning.hdf5', 
                                                          species='hg38',
                                                          save_obs=False, 
                                                          save_expect=False)
```
Then, we plot the entire heatmap with a maximum and minimum cutoff for the log2 at 4 and -4 respectively:
```Python
plot_chromosome_enrich_data(enrich_normalized_chr6, 
                            a_chr='6', 
                            bin_size=1000000, 
                            full_matrix=True, 
                            species='hg38',
                            cutoff_max=4,
                            cutoff_min=-4,
                            plot_histogram=True)
```
Heatmap             |  Histogram
:-------------------------:|:-------------------------:
![](/figures/HiCtool_chr6_1mb_normalized_enrich.png)  |  ![](/figures/HiCtool_chr6_1mb_normalized_enrich_histogram.png)

# Data normalization with explicit-factor correction model of Yaffe and Tanay

This pipeline illustrates the procedure to normalize Hi-C **intra-chromosomal contact data only** following the explicit-factor model of [Yaffe and Tanay](http://www.nature.com/ng/journal/v43/n11/abs/ng.947.html).

If you are looking for a global normalization procedure for both intra- and inter-chromosomal maps, go to [this section where the data are normalizing using a matrix balancing approach](/tutorial/normalization-matrix-balancing.md).

## Table of contents

1. [Running HiFive functions]()
   - [Creating the Fend object]()
   - [Creating the HiCData object]()
   - [Creating the HiC project object]()
   - [Filtering HiC fends]()
   - [Estimating the HiC distance function]()
   - [Learning the correction model]()
2. [Normalizing the data]()
3. [Visualizing the data]()

## 1. Running HiFive functions

The script [HiCtool_hifive.py](/scripts/HiCtool_hifive.py) is used to run all the six steps needed in order to normalize the data, whose outputs are .hdf5 files. For more information about these functions, please see [HiFive’s API documentation](http://bxlab-hifive.readthedocs.org/en/latest/api.html). To run these steps execute the following code on the Python or iPython console:
```Python
execfile('HiCtool_hifive.py')
run_hifive('restrictionsites_gc_map_valid.bed','HiCfile_pair1.bam', 'HiCfile_pair2.bam','MboI','Yaffe-Tanay')
```
More details about all the steps performed here are illustrated in the following steps 1.1.-1.6. If not interested, go to section 2.

### 1.1. Creating the Fend object

A Fragment-end (Fend) object (hdf5 format) contains information about the fragments created by digestion of a genome by a specific restriction enzyme (RE). In our script, this information is supplied in the form of a BED-formatted file (```restrictionsites_gc_map_valid.bed```) containing information about the fragment ends like coordinates, GC content and mappability score (see preprocessing, step 5).

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

The **HiC project object** (hdf5 format) links the HiCData object with information about which fends to include in the analysis, model parameters and learned model values. This is the standard way of working with Hi-C data in HiFive and this object will be used for learning the correction model and downstream analysis.

To create a HiC project object use the function ```hifive.HiC```:
```Python
import hifive

hic = hifive.HiC('HiC_project_object.hdf5', 'w')
hic.load_data('HiC_data_object.hdf5')
hic.save()
```

where ```HiC_project_object.hdf5``` specifies the location to save the HiC object to and ```HiC_data_object.hdf5``` is the data object.

### 1.4. Filtering HiC fends

At this step we filter out fragments that do not have at least one interaction before learning correction parameters. Fragment ends within a distance of 500 kb are filtered out in step 6. This will allow to normalize data without confounding technical biases with features associated to biological-relevant structures (Yaffe and Tanay).

To filter out fends use the function ```hifive.filter_fends```:
```Python
import hifive

hic = hifive.HiC('HiC_project_object.hdf5')
hic.filter_fends(mininteractions=1, mindistance=0, maxdistance=0)
hic.save()
```

### 1.5. Estimating the HiC distance function

Estimation of the **distance-dependence relationship** from the data prior to normalization, in order to avoid biases that may result due to restriction site distribution characteristics or the influence of distance/signal relationship.

Restriction sites over the genome are unevenly distributed and this results in a large set of distances between fragments and their neighbors. Since the interaction frequency is strongly inversely-related to inter-fragment distance, this means that fragments surrounded by shorter ones will show higher nearby interactions than those with longer adjacent fragments, due to the uneven distribution of the restriction sites position.

To estimate the HiC distance function use ```hifive.find_distance_parameters```:
```Python
import hifive

hic = hifive.HiC('HiC_project_object.hdf5')
hic.find_distance_parameters(numbins=90, minsize=200, maxsize=0)
hic.save()
```
### 1.6. Learning the correction model

Algorithm to learn the correction model for Hi-C data. For the normalization, we take into account of fragments length, inter-fragment distance, GC content and mappability score biases, according to the information included in the Fend object. We also consider a minimum distance of 500 kb between fragments to take into account of the effect of biological biases (TSSs and CTCF bound sites) while learning the correction parameters.

To normalize the data using the binning algorithm ([Yaffe E. and Tanay A., 2011](http://www.ncbi.nlm.nih.gov/pubmed/22001755)) use ```hic.find_binning_fend_corrections```:
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




# Data normalization with explicit-factor correction model of Yaffe and Tanay

This pipeline presents the procedure to normalize and visualize **Hi-C intra-chromosomal contact data only** following the explicit-factor model of [Yaffe and Tanay](http://www.nature.com/ng/journal/v43/n11/abs/ng.947.html). This allows also to generate the O/E contact data used to calculate the **Pearson correlation matrix** (this can be used for A/B compartment analysis reported [here](/tutorial/compartment.md)).

## Table of contents

1. [Running HiFive functions](#1-running-hifive-functions)
2. [Normalizing the data](#2-normalizing-the-data)
   - [2.1. Normalizing fend data](#21-normalizing-fend-data)
   - [2.2. Normalizing enrichment O/E data and calculating the Pearson correlation matrix](#22-normalizing-enrichment-oe-data-and-calculating-the-pearson-correlation-matrix)
3. [Visualizing the data](#3-visualizing-the-data)
   - [3.1. Visualizing the contact data](#31-visualizing-the-contact-data)
   - [3.2. Visualizing the enrichment data](#32-visualizing-the-enrichment-data)
   - [3.3. Visualizing the Pearson correlation matrix](#33-visualizing-the-pearson-correlation-matrix)

## 1. Running HiFive functions

We resort to the HiFive package, and specifically the binning algorithm, to normalize the data using the approach of Yaffe and Tanay. HiFive allows to remove spurious ligation products, as well as PCR duplicates (even though deduplication was already performed during preprocessing) and non-informative reads before generating the contact matrix.

The Python script [HiCtool_hifive.py](/scripts/HiCtool_hifive.py) is used to run all the steps needed in order to generate the data used for normalization: 

- Creating the Fend object.
- Creating the HiCData object from a Fend object and mapped data in bam format. At this step spurious ligation products (paired-reads whose total distance to their respective restriction sites exceeds 500 bp) are removed. In addition, PCR duplicates are removed and reads with ends mapping to the same fragment and reads with ends mapping to adjacent fragments on opposite strands are also excluded, to consider the possibility of incomplete restriction enzyme digestion and fragment circularization.
- Creating the HiC project object, which stores the information about the observed contact data that we will use for the downstream analysis.
- Filtering fragments that do not have at least one interaction before learning correction parameters.
- Estimating the distance-dependence relationship from the data prior to normalization. This step allows to compute the expected contact matrix based on genomic distance which is used to output the O/E contact matrix.
- Learning the correction model for Hi-C data. For the normalization, we take into account of fragments length, inter-fragment distance, GC content and mappability score biases (see parameters ``--add_gc`` and ``--add_mappability`` in the code snippet below), according to the information included in the Fend object. We also consider a minimum distance of 500 kb between fragments to take into account of the effect of biological biases (TSSs and CTCF bound sites) while learning the correction parameters.

For more information about these functions, please see [HiFive’s API documentation](http://bxlab-hifive.readthedocs.org/en/latest/api.html). To run these steps execute the following command on the Unix console (update parameters properly):
```unix
python2.7 /HiCtool-master/scripts/HiCtool_hifive.py \
-f restrictionsites.bed \
--b1 HiCfile_pair1.bam \
--b2 HiCfile_pair2.bam \
-e MboI \
-m Yaffe-Tanay \
--add_gc 1 \
--add_mappability 1
```
where:

- ``-f`` is the FEND file in bed format from preprocessing.
- ``--b1`` is the first bam file from preprocessing.
- ``--b2`` is the second bam file from preprocessing.
- ``-e`` is the restriction enzyme or enzymes names between square brackets (example ``[MboI,Hinfl]``).
- ``-m`` is the normalization model used (``Yaffe-Tanay`` in this case).
- ``--add_gc`` set to 1 to consider also the fragment GC content in the normalization process. Set to 0 if you wish not to do so.
- ``--add_mappability`` set to 1 to consider also the fragment mappability score in the normalization process. Set to 0 if you wish not to do so.

**The following output files are generated:**

- ``fend_object.hdf5``
- ``HiC_data_object.hdf5``
- ``HiC_project_object.hdf5``
- ``HiC_project_object_with distance_parameters.hdf5``
- ``HiC_norm_binning.hdf5`` to be used in the following sections.

## 2. Normalizing the data

For the normalization, observed data and correction parameters to remove biases to obtain the corrected read counts are required. Therefore, the observed contact matrix and the fend expected contact matrix are calculated. In addition, the enrichment expected contact matrix (considering the genomic distance) is calculated to compute the O/E contact matrix and the Pearson correlation matrix of the O/E matrix (to be used in the A/B compartment analysis).

For each chromosome, the following matrices can be computed and saved to file. Data are compressed in a format to reduce storage occupation and improving saving and loading time ([see here for more details](/tutorial/HiCtool_compressed_format.md)).

- The **observed contact matrix** contains the observed reads count for each bin.
- The **fend expected contain matrix** contains the learned correction value to remove biases related to fends for each bin.
- The **enrichment expected contact matrix** contains the expected reads count for each bin, considering the linear distance between read pairs and the learned correction parameters.
- The **normalized fend contact matrix** contains the corrected read count for each bin.
- The **normalized enrichment (O/E) contact matrix** contains the observed/expected value for each bin.
- The **Pearson correlation matrix** calculated from the O/E contact matrix.

**Note!**

If you need only the normalized (corrected counts) contact matrices and do not need A/B compartment analysis, there is no need to calculate also the enrichment data. If you do not need the "expected" data, do not save it since they are the biggest files and the process may take longer time.

To normalize the data, the script [HiCtool_yaffe_tanay.py](/scripts/HiCtool_yaffe_tanay.py) is used.

### 2.1. Normalizing fend data

To calculate and save the **normalized intra-chromosomal contact matrix** for a chromosome ``--chr`` do as following:
```unix
python2.7 /HiCtool-master/scripts/HiCtool_yaffe_tanay.py \
--action normalize_fend \
-i HiC_norm_binning.hdf5 \
-c /HiCtool-master/scripts/chromSizes/ \
-b 40000 \
-s hg38 \
--chr 6 \
--save_obs 1 \
--save_expect 0
```
where:

- ``--action``: Action to perform (here ``normalize_fend``).
- ``-i``:  Norm binning file in ``.hdf5`` format obtained with ``HiCtool_hifive.py`` above.
- ``-c``: Path to the folder ``chromSizes`` with trailing slash at the end ``/``.
- ``-b``: The bin size (resolution) for the analysis.
- ``-s``: Species name.
- ``--chr``: The chromosome to normalize.
- ``--save_obs``: Set to 1 to save the observed contact matrix, 0 otherwise.
- ``--save_expect``: Set to 1 to save the fend expected data with correction values, 0 otherwise.

The parameter ``--chr`` can be used also to normalize multiple chromosomes (passed as a list between square brackets) at once and also multi-processing computation is provided if your machine supports it, using the parameters ``-p`` followed by the number of processors to be used in parallel (for example ``-p 2``):
```unix
chromosomes=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y]

python2.7 /HiCtool-master/scripts/HiCtool_yaffe_tanay.py \
--action normalize_fend \
-i HiC_norm_binning.hdf5 \
-c /HiCtool-master/scripts/chromSizes/ \
-b 40000 \
-s hg38 \
--chr $chromosomes \
--save_obs 1 \
--save_expect 0
-p 24
```

### 2.2. Normalizing enrichment O/E data and calculating the Pearson correlation matrix

This analysis is generally used to calculate A/B compartments, and it is performed at coarse resolution (1 Mb typically, or 500 kb).
To calculate and save the **O/E contact matrix** and the **Pearson correlation matrix** for a chromosome ``--chr`` do as following:
```unix
python2.7 /HiCtool-master/scripts/HiCtool_yaffe_tanay.py \
--action normalize_enrich \
-i HiC_norm_binning.hdf5 \
-c /HiCtool-master/scripts/chromSizes/ \
-b 1000000 \
-s hg38 \
--chr 6 \
--save_obs 1 \
--save_expect 0
```
where:

- ``--action``: action to perform (here ``normalize_enrich``).
- ``-i``:  Norm binning file in ``.hdf5`` format obtained with ``HiCtool_hifive.py`` above.
- ``-c``: Path to the folder ``chromSizes`` with trailing slash at the end ``/``.
- ``-b``: The bin size (resolution) for the analysis.
- ``-s``: Species name.
- ``--chr``: The chromosome to normalize.
- ``--save_obs``: Set to 1 to save the observed contact matrix, 0 otherwise.
- ``--save_expect``: Set to 1 to save the expected data, 0 otherwise.

The parameter ``--chr`` can be used also to normalize multiple chromosomes (passed as a list between square brackets) at once and also multi-processing computation is provided if your machine supports it, using the parameters ``-p`` followed by the number of processors to be used in parallel (for example ``-p 2``):
```unix
chromosomes=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y]

python2.7 /HiCtool-master/scripts/HiCtool_yaffe_tanay.py \
--action normalize_enrich \
-i HiC_norm_binning.hdf5 \
-c /HiCtool-master/scripts/chromSizes/ \
-b 1000000 \
-s hg38 \
--chr $chromosomes \
--save_obs 1 \
--save_expect 0 \
-p 24
```

## 3. Visualizing the data

This section allows to plot the contact data (observed, expected or normalized fend), the enrichment O/E contact matrix and the Pearson correlation matrix. For plotting functionalities, the script [HiCtool_yaffe_tanay.py](/scripts/HiCtool_yaffe_tanay.py) is used.

### 3.1. Visualizing the contact data

This part is to plot heatmaps and histograms of the contact data. You can use this to plot observed, fend normalized contact data or expected data, either fend expected (correction values) or enrichment expected (see above for definitions).

To plot and save the heatmap and histogram of the normalized data at 40 kb resolution run the following:
```unix
python2.7 /HiCtool-master/scripts/HiCtool_yaffe_tanay.py \
--action plot_map \
-i ./yaffe_tanay_40000/chr6_40000_normalized_fend.txt \
-c /HiCtool-master/scripts/chromSizes/ \
-b 40000 \
-s hg38 \
--chr 6 \
--coord [50000000,54000000] \
--data_type normalized \
--my_colormap [white,red] \
--cutoff_type percentile \
--cutoff 95 \
--max_color "#460000" \
--plot_histogram 1
```
where:

- ``--action``: Action to perform (here ``plot_map``).
- ``-i``:  Input contact matrix.
- ``-c``: Path to the folder ``chromSizes`` with trailing slash at the end ``/``.
- ``-b``: The bin size (resolution) for the analysis.
- ``-s``: Species name.
- ``--chr``: chromosome.
- ``--coord``: List of two integers with start and end coordinates for the chromosome to be plotted.
- ``--data_type``: Data type to label your data, example: observed, normalized, expected_fend etc.
- ``--my_colormap``: Colormap to be used to plot the data. You can choose among any colorbar at https://matplotlib.org/examples/color/colormaps_reference.html, or input a list of colors if you want a custom colorbar. Example: ``[white, red, black]``. Colors can be specified also HEX format. Default: ``[white,red]``.
- ``--cutoff_type``: To select a type of cutoff (``percentile`` or ``contact``) or plot the full range of the data (not declared). Default: ``percentile``.
- ``--cutoff``: To set a maximum cutoff on the number of contacts for the colorbar based on ``--cutoff_type``. Default: 95.
- ``--max_color``: To set the color of the bins with contact counts over ``--cutoff``. Default: "#460000".
- ``--plot_histogram``: Set to 1 to plot the histogram of the data distribution, 0 otherwise.

Heatmap             |  Histogram
:-------------------------:|:-------------------------:
![](/figures/HiCtool_chr6_40kb_normalized_fend_50000000_54000000.png)  |  ![](/figures/HiCtool_chr6_40kb_normalized_fend_histogram_50000000_54000000.png)


**Additional example of the contact matrix for chromosome 6 at 1 Mb resolution**

In order to change the heatmap resolution, first data have to be normalized at the desired resolution set with the parameter ``-b`` above ([see section 2.1.](#21-normalizing-fend-data)):

Then, we plot the entire heatmap (we also change here the colormap to white and blue):
```unix
python2.7 /HiCtool-master/scripts/HiCtool_yaffe_tanay.py \
--action plot_map \
-i ./yaffe_tanay_1000000/chr6_1000000_normalized_fend.txt \
-c /HiCtool-master/scripts/chromSizes/ \
-b 1000000 \
-s hg38 \
--chr 6 \
--data_type normalized \
--my_colormap [white,blue] \
--cutoff_type percentile \
--cutoff 95 \
--max_color "#460000"
```
![](/figures/HiCtool_chr6_1mb_normalized_fend.png)

### 3.2. Visualizing the enrichment data

This part is to plot the heatmap and histogram for the enrichment O/E data ("observed over expected"). The **log2 of the data** is plotted to quantify the positive enrichment (red) and the negative enrichment (blue). Loci (pixels) equal to zero before performing the log2 (deriving from zero observed contacts) are shown in gray. Loci (pixels) where enrichment expected contact was zero before performing the ratio (observed / expected) are shown in black.

To plot and save the heatmap and histogram use the following code (maximum and minimum cutoff for the log2 at 4 and -4 respectively):
```unix
python2.7 /HiCtool-master/scripts/HiCtool_yaffe_tanay.py \
--action plot_enrich \
-i ./yaffe_tanay_1000000/chr6_1000000_normalized_enrich.txt \
-c /HiCtool-master/scripts/chromSizes/ \
-b 1000000 \
-s hg38 \
--chr 6 \
--cutoff_max 4 \
--cutoff_min -4 \
--plot_histogram 1
```
Heatmap             |  Histogram
:-------------------------:|:-------------------------:
![](/figures/HiCtool_chr6_1mb_normalized_enrich.png)  |  ![](/figures/HiCtool_chr6_1mb_normalized_enrich_histogram.png)

**Additional example of the enrichment contact matrix for chromosome 6 at 40 kb resolution**

O/E data may be plotted even at finer resolution (note that it may require some time to compute the contact matrix). In order to change the heatmap resolution, first data have to be normalized at the desired resolution set with the parameter ``-b`` above ([see section 2.2.](#21-normalizing-enrichment-data)):

```unix
python2.7 /HiCtool-master/scripts/HiCtool_yaffe_tanay.py \
--action plot_enrich \
-i ./yaffe_tanay_40000/chr6_40000_normalized_enrich.txt \
-c /HiCtool-master/scripts/chromSizes/ \
-b 40000 \
-s hg38 \
--chr 6 \
--coord [50000000,100000000] \
--plot_histogram 1
```
Heatmap             |  Histogram
:-------------------------:|:-------------------------:
![](/figures/HiCtool_chr6_40kb_normalized_enrich_50000000_54000000.png)  |  ![](/figures/HiCtool_chr6_40kb_normalized_enrich_histogram_50000000_54000000.png)

### 3.3. Visualizing the Pearson correlation matrix

This part is to plot the heatmap of the Person correlation matrix derived from the O/E matrix at 1 Mb resolution. The input file ``chr6_1000000_correlation_matrix.txt`` was generated as output file [above](#22-normalizing-enrichment-oe-data-and-calculating-the-pearson-correlation-matrix).

```unix
python2.7 /HiCtool-master/scripts/HiCtool_yaffe_tanay.py \
--action plot_correlation \
-i ./yaffe_tanay_1000000/chr6_1000000_correlation_matrix.txt \
-c /HiCtool-master/scripts/chromSizes/ \
-b 1000000 \
-s hg38 \
--chr 6
```
![](/figures/HiCtool_chr6_1mb_correlation_matrix.png) 
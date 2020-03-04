# Data normalization with explicit-factor correction model of Yaffe and Tanay

This pipeline illustrates the procedure to normalize and visualize Hi-C **intra-chromosomal contact data only** following the explicit-factor model of [Yaffe and Tanay](http://www.nature.com/ng/journal/v43/n11/abs/ng.947.html).

## Table of contents

1. [Running HiFive functions](#1-running-hifive-functions)
2. [Normalizing the data](#2-normalizing-the-data)
   - [2.1. Normalizing fend data](#21-normalizing-fend-data)
   - [2.2. Normalizing enrichment data](#22-normalizing-enrichment-data)
3. [Visualizing the data](#3-visualizing-the-data)
   - [3.1. Visualizing the contact data](#31-visualizing-the-contact-data)
   - [3.2. Visualizing the enrichment data](#32-visualizing-the-enrichment-data)

## 1. Running HiFive functions

We resort to the HiFive package, and specifically the binning algorithm, to normalize the data using the approach of Yaffe and Tanay. HiFive allows to remove spurious ligation products, as well as PCR duplicates and non-informative reads before generating the contact matrix.

The Python script [HiCtool_hifive.py](/scripts/HiCtool_hifive.py) is used to run all the steps needed in order to generate the data used for normalization: 

- Creating the Fend object.
- Creating the HiCData object from a Fend object and mapped data in bam format. At this step spurious ligation products (paired-reads whose total distance to their respective restriction sites exceeds 500 bp) are removed. In addition, PCR duplicates are removed and reads with ends mapping to the same fragment and reads with ends mapping to adjacent fragments on opposite strands are also excluded, to consider the possibility of incomplete restriction enzyme digestion and fragment circularization.
- Creating the HiC project object, which stores the information about the observed contact data that we will use for the downstream analysis.
- Filtering fragments that do not have at least one interaction before learning correction parameters.
- Estimating the distance-dependence relationship from the data prior to normalization, in order to avoid biases that may result due to restriction site distribution characteristics or the influence of distance/signal relationship. Restriction sites over the genome are unevenly distributed and this results in a large set of distances between fragments and their neighbors. Since the interaction frequency is strongly inversely-related to inter-fragment distance, this means that fragments surrounded by shorter ones will show higher nearby interactions than those with longer adjacent fragments, due to the uneven distribution of the restriction sites position.
- Learning the correction model for Hi-C data. For the normalization, we take into account of fragments length, inter-fragment distance, GC content and mappability score biases, according to the information included in the Fend object. We also consider a minimum distance of 500 kb between fragments to take into account of the effect of biological biases (TSSs and CTCF bound sites) while learning the correction parameters.

For more information about these functions, please see [HiFive’s API documentation](http://bxlab-hifive.readthedocs.org/en/latest/api.html). To run these steps execute the following command on the Unix console (update parameters properly):
```unix
python /HiCtool-master/scripts/HiCtool_hifive.py \
-f restrictionsites.bed \
--b1 HiCfile_pair1.bam \
--b2 HiCfile_pair2.bam \
-e MboI \
-m Yaffe-Tanay
```
where:

- ``-f`` is the FEND file in bed format from preprocessing.
- ``--b1`` is the first bam file from preprocessing.
- ``--b2`` is the second bam file from preprocessing.
- ``-e`` is the restriction enzyme or enzymes names between square brackets (example [MboI,Hinfl]).
- ``-m`` is the normalization model used (Yaffe_Tanay in this case).

**The following output files are generated:**

- ``fend_object.hdf5``
- ``HiC_data_object.hdf5``
- ``HiC_project_object.hdf5``
- ``HiC_norm_binning.hdf5`` to be used in the following section.


## 2. Normalizing the data

For the normalization, observed data and correction parameters to remove biases to obtain the corrected read counts are required. Therefore, the observed contact matrix and the fend expected contact matrix are calculated. In addition, the enrichment expected contact matrix is calculated to compute the observed over expected enrichment values, considering also the distance between fends.

For each chromosome, the following five matrices can be computed and saved to file. Data are compressed in a format to reduce storage occupation and improving saving and loading time ([see here for more details](/tutorial/HiCtool_compressed_format.md)).

- The **observed data** contain the observed reads count for each bin.
- The **fend expected data** contain the learned correction value to remove biases related to fends for each bin.
- The **enrichment expected data** contain the expected reads count for each bin, considering the linear distance between read pairs and the learned correction parameters.
- The **normalized fend data** contain the corrected read count for each bin.
- The **normalized enrichment data** ("observed over expected" matrix) contain the enrichment value (O/E) for each bin.

**Note!**
If you need only the normalized contact matrices, there is no need to calculate also the enrichment data. If you do not need the expected data, do not save it since they are the biggest files and the process may take time.

To normalize the data, the script [HiCtool_yaffe_tanay.py](/scripts/HiCtool_yaffe_tanay.py) is used.

### 2.1. Normalizing fend data

To calculate and save the **normalized intra-chromosomal contact matrix** for a chromosome ``--chr`` do as following:
```unix
python /HiCtool-master/scripts/HiCtool_yaffe_tanay.py \
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

- ``--action``: action to perform (here ``normalize_fend``).
- ``-i``:  Norm binning file in ``.hdf5`` format obtained with ``HiCtool_hifive.py`` above.
- ``-c``: Path to the folder ``chromSizes`` with trailing slash at the end ``/``.
- ``-b``: The bin size (resolution) for the analysis.
- ``-s``: Species name.
- ``--chr``: The chromosome to normalize.
- ``--save_obs``: Set to 1 to save the observed contact matrix, 0 otherwise.
- ``--save_expect``: Set to 1 to save the fend expected data with correction values, 0 otherwise.

The parameter ``--chr`` can be used also to normalize multiple chromosomes (passed as a list between square brackets) at once and also multi-processing computation is provided if your machine supports it using the parameters ``-p``:
```unix
chromosomes=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y]

python /HiCtool-master/scripts/HiCtool_yaffe_tanay.py \
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

### 2.2. Normalizing enrichment data

To calculate and save the **"observed/expected" intra-chromosomal contact matrix** for a chromosome ``--chr`` do as following:
```unix
python /HiCtool-master/scripts/HiCtool_yaffe_tanay.py \
--action normalize_enrich \
-i HiC_norm_binning.hdf5 \
-c /HiCtool-master/scripts/chromSizes/ \
-b 40000 \
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

The parameter ``--chr`` can be used also to normalize multiple chromosomes (passed as a list between square brackets) at once and also multi-processing computation is provided if your machine supports it using the parameters ``-p``:
```unix
chromosomes=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y]

python /HiCtool-master/scripts/HiCtool_yaffe_tanay.py \
--action normalize_enrich \
-i HiC_norm_binning.hdf5 \
-c /HiCtool-master/scripts/chromSizes/ \
-b 40000 \
-s hg38 \
--chr $chromosomes \
--save_obs 1 \
--save_expect 0 \
-p 24
```

## 3. Visualizing the data

This section allows to plot the contact data (observed, expected or normalized fend) and the enrichment "observed over expected" data.

For plotting functionalities, the script [HiCtool_yaffe_tanay.py](/scripts/HiCtool_yaffe_tanay.py) is used.

### 3.1. Visualizing the contact data

This part is to plot heatmaps and histograms of the contact data. You can use this to plot observed, fend normalized contact data or expected data, either fend expected or enrichment expected (see above for definitions).

To plot and save the heatmap and histogram of the normalized data at 40 kb resolution run the following:
```unix
python /HiCtool-master/scripts/HiCtool_yaffe_tanay.py \
--action plot_map \
-i HiCtool_chr6_40kb_normalized_fend.txt \
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

- ``--action``: action to perform (here ``plot_map``).
- ``-i``:  Input contact matrix.
- ``-c``: Path to the folder ``chromSizes`` with trailing slash at the end ``/``.
- ``-b``: The bin size (resolution) for the analysis.
- ``-s``: Species name.
- ``--chr``: chromosome.
- ``--coord``: List of two integers with start and end coordinates for the chromosome to be plotted.
- ``--data_type``: Data type to label your data, example: observed, normalized, expected_fend etc.
- ``--my_colormap``: Colormap to be used to plot the data. You can choose among any colorbar at https://matplotlib.org/examples/color/colormaps_reference.html, or input a list of colors if you want a custom colorbar. Example: [white, red, black]. Colors can be specified also HEX format. Default: [white,red].
- ``--cutoff_type``: To select a type of cutoff (percentile or contact) or plot the full range of the data (not declared). Default: percentile.
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
python /HiCtool-master/scripts/HiCtool_yaffe_tanay.py \
--action plot_map \
-i HiCtool_chr6_1mb_normalized_fend.txt \
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

This part is to plot the heatmap and histogram for the enrichment normalized data ("observed over expected"). The **log2 of the data** is plotted to quantify the positive enrichment (red) and the negative enrichment (blue). Loci (pixels) equal to zero before performing the log2 (deriving from zero observed contacts) are shown in gray. Loci (pixels) where enrichment expected contact was zero before performing the ratio (observed / expected) are shown in black.

To plot and save the heatmap and histogram use the following code:
```unix
python /HiCtool-master/scripts/HiCtool_yaffe_tanay.py \
--action plot_enrich \
-i HiCtool_chr6_40kb_normalized_enrich.txt \
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

**Additional example of the enrichment contact matrix for chromosome 6 at 1 Mb resolution**

In order to change the heatmap resolution, first data have to be normalized at the desired resolution set with the parameter ``-b`` above ([see section 2.2.](#21-normalizing-enrichment-data)):

Then, we plot the entire heatmap with a maximum and minimum cutoff for the log2 at 4 and -4 respectively:
```unix
python /HiCtool-master/scripts/HiCtool_yaffe_tanay.py \
--action plot_enrich \
-i HiCtool_chr6_1mb_normalized_enrich.txt \
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

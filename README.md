# HiCtool

HiCtool is a Python library for processing and visualizing Hi-C data, including topological domain analysis.  
This is a short documentation to give an overview of the tool, the full documentation and code are available at [http://www.genomegitar.org](https:genomegitar.org).  

The cell line used here is from human (hg38) with GEO accession number [GSM1551550](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1551550). To run the example in this short tutorial you do not need to have installed all the following software/libraries but only:

- [numpy](http://scipy.org/)
- [scipy](http://scipy.org/)
- [matplotlib](http://matplotlib.org/)
- [math](https://docs.python.org/2/library/math.html)

To check if a module is installed, open your **python console** and try to import each module by typing:
```python
import my_module
```
If a module is not installed, go to your **unix console** and type the following to install it:
```unix
python install my_module
```

## Installation

HiCtool is in a pipeline format to allow extreme flexibility and easy usage. You do not need to install anything besides the following Python libraries, packages and software. Everything is open source.

**1. Python libraries [for python >2.7]:**

- [numpy](http://scipy.org/)
- [scipy](http://scipy.org/)
- [math](https://docs.python.org/2/library/math.html)
- [matplotlib](http://matplotlib.org/)
- [matplotlib.pyplot](http://matplotlib.org/api/pyplot_api.html#module-matplotlib.pyplot)
- [csv](https://docs.python.org/2/library/csv.html)
- [pybedtools](https://daler.github.io/pybedtools/)
- [pandas](https://pandas.pydata.org/)
- [multiprocessing](https://docs.python.org/2/library/multiprocessing.html)
- [biopython](http://biopython.org/)

**2. Python packages:**

- [hifive](http://bxlab-hifive.readthedocs.org/en/latest/introduction.html)
- [hmmlearn](https://github.com/hmmlearn/hmmlearn)

**3. Other software needed (for preprocessing only):**

- [BEDTools](http://bedtools.readthedocs.org/en/latest/)
- [Bowtie 2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
- [SAMTools](http://samtools.sourceforge.net/)
- [SRA Toolkit](http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc&f=fastq-dump)

## Data preprocessing

HiCtool provides a complete pipeline from the downloading of the raw data (SRA format) to the final BAM files that are used for the following analysis steps. In addition, instructions on how to generate a fragment end BED file to correct biases are provided.



Preprocessing steps:

1. Downloading the source data from GEO.
2. Pre-truncation of the reads that contain potential ligation junctions.
3. Mapping read pairs to the reference genome.
4. Filtering reads and selecting reads that are paired.
5. Creating the fragment-end (FEND) bed file.

These are the two BAM files and the FEND bed file if you want to proceed with the entire data normalization pipeline of the [full documentation](https://doc.genomegitar.org/data_analysis_and_visualization.html) (not reported here):

- [BAM file pair 1](https://sysbio.ucsd.edu/public/rcalandrelli/hictool_example/HiCfile_pair1.bam)
- [BAM file pair 2](https://sysbio.ucsd.edu/public/rcalandrelli/hictool_example/HiCfile_pair2.bam)
- [FEND bed file](https://sysbio.ucsd.edu/public/rcalandrelli/hictool_example/MboI_hg38_gc_map_valid.zip)

Download the data using ```wget``` followed by the link address of each file.

## Data analysis and visualization

The data analysis and visualization section provides the pipeline to normalize the data and plot the heatmaps. The normalization has been done using the Python package [HiFive](http://bxlab-hifive.readthedocs.org/en/latest/introduction.html) while for plotting Matplotlib is used, with the possibility also to add a histogram of the distribution of the data. Both observed and normalized counts can be plotted. In addition, we provide the possibility of plotting “observed over expected” (enrichment) contact heatmaps, where the expected counts are calculated considering both the learned correction parameters and the distance between read pairs, given the property that the average intrachromosomal contact probability for pairs of loci decreases monotonically with increasing of their linear genomic distance.

Data analysis and visualization steps:

1. Creating the Fend object.
2. Creating the HiCData object.
3. Creating the HiC project object.
4. Filtering HiC fends.
5. Estimating the HiC distance function.
6. Learning the correction model.
7. Normalizing the data.
8. Visualizing the data.

After the data are normalized, if both fend and enrichment data were calculated, these files will be produced (here for chromosome 6, at 1 mb and 40 kb resolution):

- [Observed data chr 6 (1 mb)](https://sysbio.ucsd.edu/public/rcalandrelli/hictool_example/HiCtool_chr6_1mb_observed.txt)
- [Normalized fend data chr 6 (1 mb)](https://sysbio.ucsd.edu/public/rcalandrelli/hictool_example/HiCtool_chr6_1mb_normalized_fend.txt)
- [Observed data chr 6 (40 kb)](https://sysbio.ucsd.edu/public/rcalandrelli/hictool_example/HiCtool_chr6_40kb_observed.txt)
- [Normalized fend data chr 6 (40 kb)](https://sysbio.ucsd.edu/public/rcalandrelli/hictool_example/HiCtool_chr6_40kb_normalized_fend.txt)
- [Normalized enrichment data chr 6 (40 kb)](https://sysbio.ucsd.edu/public/rcalandrelli/hictool_example/HiCtool_chr6_40kb_normalized_enrich.txt)

To download the data use the command ```wget``` on the **unix shell** and copy the url link from a file above, here an example for the normalized fend and enrichment data chr 6 (40 kb):

```unix
wget https://sysbio.ucsd.edu/public/rcalandrelli/hictool_example/HiCtool_chr6_40kb_normalized_fend.txt

wget https://sysbio.ucsd.edu/public/rcalandrelli/hictool_example/HiCtool_chr6_40kb_normalized_enrich.txt
```

Then download the Python script [HiCtool_normalization_visualization.py](https://sysbio.ucsd.edu/public/rcalandrelli/hictool_example/HiCtool_normalization_visualization.py):

```unix
wget https://sysbio.ucsd.edu/public/rcalandrelli/hictool_example/HiCtool_normalization_visualization.py
```

Open a Python or iPython console on the unix shell:
```unix
# Open a Python console
python
```
Execute the script:
```Python
execfile("HiCtool_normalization_visualization.py")
```
Plot the normalized fend data:
```Python
plot_chromosome_fend_data('HiCtool_chr6_40kb_normalized_fend.txt', a_chr='6', bin_size=40000, full_matrix=False, start_coord=50000000, end_coord=54000000, species='hg38', my_colormap=['white', 'red'], cutoff_type='percentile', cutoff=95, max_color='#460000', plot_histogram=True)
```
Here we plot data of chromosome 6, from 50 Mb to 54 Mb at a bin size of 40 kb, for species hg38. We use a colormap which goes from white (no contacts) to red (maximum contact) and we use a upper cut-off of the 95th percentile of the data to enhance higher order chromatin structure such as topological domains. We assign to the bins over the cut-off a spefic color (#460000) and also we choose to plot the distribution of the contact data as well.

The same can be done for the "observed over expected" data:
```Python
plot_chromosome_enrich_data('HiCtool_chr6_40kb_normalized_enrich.txt', a_chr='6', bin_size=40000, full_matrix=False, start_coord=50000000, end_coord=54000000, species='hg38', plot_histogram=True)
```
Or the normalized fend data at 1 mb resolution:
```Python
plot_chromosome_fend_data('HiCtool_chr6_1mb_normalized_fend.txt', a_chr='6', bin_size=1000000, full_matrix=True, species='hg38', my_colormap=['white', 'blue'], cutoff_type='percentile', cutoff=95, max_color='#460000', plot_histogram=True)
```
In this case we plot the entire contact matrix and we changed the color of the heatmap to blue using the parameter ```my_colormap```.


## Topological domain analysis



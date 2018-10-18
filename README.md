# HiCtool

HiCtool is a Python library for processing and visualizing Hi-C data, including topological domain analysis.  
This is a short documentation to give an overview of the tool, the full documentation and code are available at [http://www.genomegitar.org](https:genomegitar.org).

## Installation

HiCtool is in a pipeline format to allow extreme flexibility and easy usage. You do not need to install anything besides the following Python libraries, packages and software. Everything is open source.

**1. Python libraries [for python >2.7]:**

- [Numpy](http://scipy.org/)
- [Scipy](http://scipy.org/)
- [Math](https://docs.python.org/2/library/math.html)
- [Matplotlib](http://matplotlib.org/)
- [Matplotlib.pyplot](http://matplotlib.org/api/pyplot_api.html#module-matplotlib.pyplot)
- [csv](https://docs.python.org/2/library/csv.html)
- [pybedtools](https://daler.github.io/pybedtools/)
- [pandas](https://pandas.pydata.org/)
- [multiprocessing](https://docs.python.org/2/library/multiprocessing.html)
- [Biopython](http://biopython.org/)

**2. Python packages:**

- [HiFive](http://bxlab-hifive.readthedocs.org/en/latest/introduction.html)
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


## Data analysis and visualization

The data analysis and visualization section provides the pipeline to normalize the data and plot the heatmaps. The normalization has been done using the Python package [HiFive](http://bxlab-hifive.readthedocs.org/en/latest/introduction.html) while for plotting Matplotlib is used, with the possibility also to add a histogram of the distribution of the data. Both observed and normalized counts can be plotted. In addition, we provide the possibility of plotting “observed over expected” contact heatmaps, where the expected counts are calculated considering both the learned correction parameters and the distance between read pairs, given the property that the average intrachromosomal contact probability for pairs of loci decreases monotonically with increasing of their linear genomic distance.

Data analysis and visualization steps:

1. Creating the Fend object.
2. Creating the HiCData object.
3. Creating the HiC project object.
4. Filtering HiC fends.
5. Estimating the HiC distance function.
6. Learning the correction model.
7. Normalizing the data.
8. Visualizing the data.

After the data are normalized, if both fend and enrichment data were calculated, these files will be available (here for chromosome 6 only at 1 mb and 40 kb resolution for GEO accession number GMS1551550):

- [Observed data chr 6 (1 mb)](https://sysbio.ucsd.ed/public/rcalandrelli/hictool_example/HiCtool_chr6_1mb_observed.txt)
- [Normalized fend data chr 6 (1 mb)]()
- [Normalized enrichment data chr 6 (1 mb)]()

- [Observed data chr 6 (40 kb)]()
- [Normalized fend data chr 6 (40 kb)]()
- [Normalized enrichment data chr 6 (40 kb)]()

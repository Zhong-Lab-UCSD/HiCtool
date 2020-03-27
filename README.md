# HiCtool: a standardized pipeline to process and visualize Hi-C data (v2.2)

HiCtool is an open-source bioinformatic tool based on Python, which integrates several software to perform a standardized Hi-C data analysis, from the processing of raw data, to the visualization of heatmaps and the identification of topologically associated domains (TADs) and A/B compartments.

## Table of Contents

- [Overview](#overview)
- [Installation](#installation)
- [Tutorial](#tutorial)
- [Reference](#reference)
- [Version history](#version-history)
- [Support](#support)

## Overview

HiCtool leads the user step-by-step through a pipeline which consists of the following sections:

- Data preprocessing
- Data normalization and visualization
- A/B compartment analysis
- TAD analysis

HiCtool can implement contact data normalization following two approaches: 

- **The explicit-factor correction method reported by [Yaffe and Tanay](https://www.nature.com/articles/ng.947)** and performed by the library [HiFive](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0806-y). In this case, only intra-chromosomal analysis is performed, per each chromosome singularly and only single heatmaps can be plotted. In addition, there is the possibility to plot topological domains over the heatmap at a resolution of 40kb or lower. **This approach allows also to compute the O/E contact matrix, derive the Pearson correlation matrix and perform A/B compartment analysis**.
- **The matrix balancing approach performed by [Hi-Corrector](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4380031/).** In this case, a global analysis is performed including all the chromosomes and both intra- and inter-chromosomal maps. It is possible to visualize either single intra- and inter-chromosomal heatmap or the global all-by-all chromosomes heatmap (for the global heatmap visualization, resolution could be limited by your hardware). In addition, there is the possibility to plot topological domains over the intra-chromosomal heatmap (resolution of 40kb or lower) or plot the same maps from different samples on a side-by-side view for easy comparison.

## Installation

HiCtool is in a pipeline format based on single unix commands to allow easy usage. In order to use HiCtool, you need to install the following Python libraries, packages and software. Everything is open source. After that, you need the HiCtool source codes. **[Click here to download HiCtool](https://github.com/Zhong-Lab-UCSD/HiCtool/archive/master.zip)**, unzip the file, all the scripts are under the folder ``scripts``. [Hi-Corrector](http://zhoulab.usc.edu/Hi-Corrector/) source code is already inside this folder.

**1. Python libraries [for Python 2.7]:**

- [numpy](http://scipy.org/)
- [scipy](http://scipy.org/)
- [matplotlib](http://matplotlib.org/)
- [math](https://docs.python.org/2/library/math.html)
- [matplotlib.pyplot](http://matplotlib.org/api/pyplot_api.html#module-matplotlib.pyplot/)
- [csv](https://docs.python.org/2/library/csv.html)
- [pybedtools](https://daler.github.io/pybedtools/)
- [pandas](https://pandas.pydata.org/)
- [multiprocessing](https://docs.python.org/2/library/multiprocessing.html)
- [biopython](http://biopython.org/)
- [sklearn-decomposition](https://scikit-learn.org/)

**2. Python packages:**

- [hifive (v1.4)](http://bxlab-hifive.readthedocs.org/en/latest/introduction.html)
- [hmmlearn](https://github.com/hmmlearn/hmmlearn)

**3. Other software:**

- [BEDTools](http://bedtools.readthedocs.org/en/latest/)
- [Bowtie 2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
- [SAMTools](http://samtools.sourceforge.net/)
- [SRA Toolkit](http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc&f=fastq-dump)

## Tutorial

We have compiled a full tutorial to show the usage of the pipeline. Please check the [Tutorial Homepage](./tutorial/ReadMe.md).

## Reference

HiCtool was developed by Riccardo Calandrelli and Qiuyang Wu from Dr. Sheng Zhong's Lab at University of California San Diego. 

If you use HiCtool, please cite the paper: 

Calandrelli, R., Wu, Q., Guan, J., & Zhong, S. (2018). GITAR: An open source tool for analysis and visualization of Hi-C data. *Genomics, proteomics & bioinformatics.*

## Version history

### March 27, 2020

- Version 2.2 released:

   - A/B compartment analysis is now included.
   - Calculation and plotting of the Pearson correlation matrix from the O/E contact matrix (only Yaffe-Tanay method).
   - Calculation and plotting of the eigenvector for A/B compartment annotation.
   - Small bug fixes.

### July 25, 2019

- Version 2.1.1 released:

   - Optimized pre-processing and normalization to allow the analysis of very high depth sequenced files.
   - Added parameter to allow custom mapped reads filtering at the pre-processing step.
   - Bedpe file generation step added for pre-processed read pairs.
   - Deduplication step added in the pre-processing pipeline.
   - Optimized and faster FEND file generation for adding GC content and mappability score.
   - Small bug fixes.

### May 20, 2019

- Version 2.1 released:

   - HiCtool is now based only on unix commands. A Python script is given which contains utility functions for I/O of files generated with HiCtool in the Python console.
   - The entire pre-processing is now performed with a single-command.
   - Possibility of processing Hi-C data generated with a cocktail of restriction enzymes (i.e. as in the Arima Kit) or with restriction enzymes including "N" in their sequence.
   - HiCtool includes now additional species besides hg38 and mm10: hg19, mm9, dm6, susScr3, susScr11.
   - New function to visualize heatmaps from different samples or conditions on a side-by-side view for comparison.
   - Small bug fixes.

### March 28, 2019

- Version 2.0 released:

   - HiCtool code is now hosted on GitHub.
   - Added inter-chromosomal analysis and visualization.
   - Included an additional global normalization method based on a matrix balancing approach.
   - New function to plot the all-by-all chromosomes global contact matrix.
   - Possibility of saving contact matrices in tab separated format.
   - Possibility of plotting topological domains over the heatmaps.
   - Small bug fixes.

### December 2015 - October 2018

- The initial release of HiCtool v1.0 came out in late 2015. GITAR manuscript (including HiCtool) published in October 2018.

## Support

For issues related to the use of HiCtool or if you want to report a bug, please contact Riccardo Calandrelli at <rcalandrelli@eng.ucsd.edu>.

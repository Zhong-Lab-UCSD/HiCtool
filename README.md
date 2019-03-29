# HiCtool: a standardized pipeline to process and visualize Hi-C data (v2.0)

HiCtool is an open-source bioinformatic tool based on Python, which integrates several software to perform a standardized Hi-C data analysis, from the processing of raw data, to the visualization of heatmaps and the identification of topologically associated domains (TADs).

## Table of Contents

- [Overview](#overview)
- [Installation](#installation)
- [Tutorial](#tutorial)
- [API documentation](#api-documentation)
- [Version history](#version-history)
- [Reference](#reference)
- [Support](#support)

## Overview

We implemented a pipeline that is divided into three main sections:

- Data preprocessing
- Data normalization and visualization
- TAD analysis

HiCtool leads the user step-by-step through a pipeline, which goes from the raw Hi-C data to the computation, visualization, and optimized storage of contact matrices (intra- and inter-chromosomal) and TAD coordinates. 

HiCtool can implement contact data normalization following two approaches: 

- The explicit-factor correction method reported by [Yaffe and Tanay](https://www.nature.com/articles/ng.947) and performed by the library [HiFive](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0806-y). In this case, only intra-chromosomal analysis is performed, per each chromosome singularly.
- The matrix balancing approach performed by [Hi-Corrector](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4380031/). In this case, a global analysis is performed including all the chromosomes and both intra- and inter-chromosomal maps.

## Installation

HiCtool is in a pipeline format to allow extreme flexibility and easy usage. In order to use HiCtool, you need to install the following Python libraries, packages and software. Everything is open source.

**1. Python libraries [for python>2.7]:**

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

**2. Python packages:**

- [hifive](http://bxlab-hifive.readthedocs.org/en/latest/introduction.html)
- [hmmlearn](https://github.com/hmmlearn/hmmlearn)

**3. Other software:**

- [BEDTools](http://bedtools.readthedocs.org/en/latest/)
- [Bowtie 2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
- [SAMTools](http://samtools.sourceforge.net/)
- [SRA Toolkit](http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc&f=fastq-dump)
- [Hi-Corrector](http://zhoulab.usc.edu/Hi-Corrector/)

## Tutorial

We have compiled a full tutorial to show the usage of the pipeline. Please check the [Tutorial Homepage](./tutorial/ReadMe.md).

## API documentation

All the functions used in this documentation are reported with all the input parameters in the [API documentation](https://sysbio.ucsd.edu/public/rcalandrelli/HiCtool_API_documentation.pdf).

## Version history

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

- The initial release of HiCtool came out in late 2015. GITAR manuscript (including HiCtool) published in October 2018.

## Reference

HiCtool was developed by Riccardo Calandrelli and Qiuyang Wu from Dr. Sheng Zhong's Lab at University of California, San Diego. 

If you use HiCtool, please cite the paper: 

Calandrelli, R., Wu, Q., Guan, J., & Zhong, S. (2018). GITAR: An open source tool for analysis and visualization of Hi-C data. *Genomics, proteomics & bioinformatics.*

## Support

For issues related to the use of HiCtool or if you want to report a bug, please contact Riccardo Calandrelli at <rcalandrelli@eng.ucsd.edu>.

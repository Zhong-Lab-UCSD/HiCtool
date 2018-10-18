# HiCtool

HiCtool is a Python library for processing and visualizing Hi-C data, including topological domain analysis.  
The full documentation is available at [http://www.genomegitar.org](https:genomegitar.org).

## Installation

HiCtool is in a pipeline format to allow extreme flexibility and easy usage. You do not need to install anything besides the following Python libraries, packages and software. Everything is open source.

1. Python libraries [for python >2.7]:

- Numpy
- Scipy
Math
Matplotlib
Matplotlib.pyplot
csv
pybedtools
pandas
multiprocessing
Biopython
2. Python packages:

HiFive
hmmlearn
3. Other software needed:

BEDTools
Bowtie 2
SAMTools
SRA Toolkit

## Data preprocessing

1. Downloading the source data from GEO.
2. Pre-truncation of the reads that contain potential ligation junctions.
3. Mapping read pairs to the reference genome.
4. Filtering reads and selecting reads that are paired.
5. Creating the fragment-end (FEND) bed file.

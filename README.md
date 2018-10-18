# HiCtool

HiCtool is a Python library for processing and visualizing Hi-C data, including topological domain analysis.  
The full documentation is available at [http://www.genomegitar.org](https:genomegitar.org).

## Data preprocessing

1. Downloading the source data from GEO.
2. Pre-truncation of the reads that contain potential ligation junctions.
3. Mapping read pairs to the reference genome.
4. Filtering reads and selecting reads that are paired.
5. Creating the fragment-end (FEND) bed file.

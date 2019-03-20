# Data preprocessing

## Table of Contents

1. [Downloading the raw data from GEO](#1-downloading-the-raw-data-from-geo)
2. [Pre-truncation of the reads that contain potential ligation junctions](#2-pre-trunction-of-the-reads-that-contain-potential-ligation-junction)
3. Mapping read pairs to the reference genome.
4. Filtering reads and selecting reads that are paired.
5. Creating the fragment-end (FEND) bed file.

## 1. Downloading the raw data from GEO

The source data in sra format are downloaded via GEO accession number using the command ```fastq-dump``` of [SRA Toolkit](http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc&f=fastq-dump).

Before proceeding, you may need to [setup the output directory](https://github.com/ncbi/sra-tools/wiki/Toolkit-Configuration) where the sra files will be saved. After having installed SRA Toolkit, go to the path where the software has been installed, under the subfolder “bin”, and run the following command line:
```unix
./vdb-config -i
```
This will open an interface that will allow you to setup/change your output directory.

To download the data related to a GEO accession number, go to the bottom of that page and click on the SRA number under the section “Relations”. After that, under the section “Runs” you will find the SRR files, then run the following:
```unix
fastq-dump SRRXXXXXXX --split-3
```
where ```SRRXXXXXXX``` has to be replaced with the specific number of the run you want to download (**[SRR1658570](https://www.ncbi.nlm.nih.gov/sra?term=SRX764936) in this documentation**).

To be more specific, this code will either download the SRA file under the output directory you have set with the interface above, but it will also convert the SRA file into fastq files and dump each read into separate files in the current working directory, to produce:

- SRRXXXXXXX_1.fastq
- SRRXXXXXXX_2.fastq
- SRRXXXXXXX.fastq

where paired-end reads in SRRXXXXXXX.sra are split and stored into **SRRXXXXXXX_1.fastq** and **SRRXXXXXXX_2.fastq**, and SRRXXXXXXX.fastq (if present) contains reads with no mates.

**Note!**
To produce our final results, use this GEO accession number: **[GSM1551550](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1551550)**.

## 2. Pre-trunction of the reads that contain potential ligation junction

After the fastq files are obtained, pre-truncation is performed on the reads that contain potential ligation junctions to keep the longest piece without a junction sequence ([Ay et al., 2015](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0745-7)). To do so, use the code in [pre_truncation.py](/scripts/pre_truncation.py) and run the following code on your Python or iPython console:
```Python
execfile('pre_truncation.py')
pre_truncation('SRR1658570_1.fastq', 'MboI')
pre_truncation('SRR1658570_2.fastq', 'MboI')
```
The output files will have the same filename with extension ```.trunc.fastq```. If a different restriction enzyme (RE) than HindIII, MboI, NcoI and DpnII has been used in the Hi-C experiment, then run:
```Python
execfile('pre_truncation.py')
pre_truncation('a_file.fastq', 'custom_RE', 'custom_ligation_junction')
```
where ```custom_ligation_junction``` is the ligation junction sequence of nucleotides associated to the restriction enzyme you are using (see [Ay et al., 2015](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0745-7) for more details).

A log file named **pre_truncation_log.txt** is generated with the information about the percentage of reads that have been truncated. This is also printed on the console:
```Python
SRR1658570_1.fastq
202095066 reads (length = 101 bp); of these:
29851195 (14.78%) contained a potential ligation junction and have been truncated.
SRR1658570_2.fastq
202095066 reads (length = 101 bp); of these:
28681691 (14.2%) contained a potential ligation junction and have been truncated.
```
The length distribution of the truncated reads is also plotted and saved to file.

![](/figures/SRR1658570_1.fastq_truncated_reads.png =250x)


# Data preprocessing

## Table of Contents

1. Downloading the raw data from GEO.
2. Pre-truncation of the reads that contain potential ligation junctions.
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

After the fastq files are obtained, pre-truncation is performed on the reads that contain potential ligation junctions to keep the longest piece without a junction sequence ([Ay et al., 2015](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0745-7)). To do so, use the code in [pre_truncation.py](/scripts/pre_truncation.py) (see API Documentation) and run the following code on your Python or iPython console:


# Data preprocessing

## Table of Contents

1. [Downloading the raw data from GEO](#1-downloading-the-raw-data-from-geo)
2. [Pre-truncation of the reads that contain potential ligation junctions](#2-pre-trunction-of-the-reads-that-contain-potential-ligation-junction)
3. [Mapping read pairs to the reference genome](#3-mapping-read-pairs-to-reference-genome)
4. [Filtering reads and selecting reads that are paired](#4-filtering-reads-and-selecting-reads-that-are-paired)
5. [Creating the fragment-end (FEND) bed file](#5-creating-the-fragment-end-fend-bed-file)

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

![](/figures/SRR1658570_1.fastq_truncated_reads.png)

![](/figures/SRR1658570_2.fastq_truncated_reads.png)

## 3. Mapping read pairs to reference genome

[Bowtie 2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) is used for mapping the read pairs, and reads are mapped **independently** to avoid any proximity constraint. To align the reads, first build a **corresponding index** for the reference genome (execute this step only one time for reference genome):
```unix
bowtie2-build hg38.fa index
```
```hg38.fa``` is the reference sequence in FASTA format, the output files in bt2 format are named with the prefix 'index'.

Now align the reads to the reference sequence:
```unix
(bowtie2 -p 32 -x index SRR1658570_1.trunc.fastq -S HiCfile1.sam) 2>HiCfile1_log.txt
(bowtie2 -p 32 -x index SRR1658570_2.trunc.fastq -S HiCfile2.sam) 2>HiCfile2_log.txt
```
where:

- The ```-p``` argument refers to a specified number of parallel search threads (update it accordingly to your number of available cores). This can be useful to decrease the processing time in aligning the reads.
- The ```-x``` argument specifies the basename of the index for the reference genome. The basename is the name of any of the index files up to but not including the final ‘/.1.bt2’, ‘/.2.bt2’, etc.
- The ```-S``` argument specifies the output file in sam format.
- ```HiCfile1_log.txt``` and ```HiCfile2_log.txt``` are log files containing the statistics of the alignment:

```unix
HiCfile1_log.txt
202095066 reads; of these:
202095066 (100.00%) were unpaired; of these:
5770798 (2.86%) aligned 0 times
156759009 (77.57%) aligned exactly 1 time
39565259 (19.58%) aligned >1 times
97.14% overall alignment rate

HiCfile2_log.txt
202095066 reads; of these:
202095066 (100.00%) were unpaired; of these:
13381441 (6.62%) aligned 0 timess
149852422 (74.15%) aligned exactly 1 time
38861203 (19.23%) aligned >1 times
93.38% overall alignment rate
```

## 4. Filtering reads and selecting reads that are paired

[SAMtools](http://samtools.sourceforge.net/) is used to extract the headers and perform filtering on the mapped reads ([Yaffe and Tanay, 2011](http://www.nature.com/ng/journal/v43/n11/abs/ng.947.html)):
```unix
samtools view -H HiCfile1.sam > header1.txt
samtools view -H HiCfile2.sam > header2.txt

samtools view -F 4 -q 30 HiCfile1.sam > HiCfile1_hq.sam
samtools view -F 4 -q 30 HiCfile2.sam > HiCfile2_hq.sam
```
- ```-F 4``` is used to discard unmapped reads.
- ```-q 30``` is used to select reads that were uniquely mapped with a MAPQ >= 30, i.e. the estimated probability of mapping error is <= 0.1%.

**Unused files are removed** using the command ```rm``` to reduce memory occupance since they may be quite big. If you want to keep your files, do not run those commands.

To update the log files from mapping with information about reads mapped with MAPQ >= 30, use the following code:
```unix
n1=`wc -l HiCfile1_hq.sam | awk '{print $1}'`
n2=`wc -l HiCfile2_hq.sam | awk '{print $1}'`

nt1=`wc -l HiCfile1.sam | awk '{print $1}'`
h1=`wc -l header1.txt | awk '{print $1}'`
ntot1=`expr $nt1 - $h1`

nt2=`wc -l HiCfile2.sam | awk '{print $1}'`
h2=`wc -l header2.txt | awk '{print $1}'`
ntot2=`expr $nt2 - $h2`

perc1=$(awk -v n1=$n1 -v ntot1=$ntot1 'BEGIN { print 100*n1/ntot1 }' | cut -c1-5)
perc2=$(awk -v n2=$n2 -v ntot2=$ntot2 'BEGIN { print 100*n2/ntot2 }' | cut -c1-5)

printf "\n----------\n"$ntot1" reads; of these:\n  "$n1" ("$perc1"%%) aligned with MAPQ>=30" >> HiCfile1_log.txt
printf "\n----------\n"$ntot2" reads; of these:\n  "$n2" ("$perc2"%%) aligned with MAPQ>=30" >> HiCfile2_log.txt

rm HiCfile1.sam
rm HiCfile2.sam
```
After filtering, reads that are not paired from the two SAM files (HiCfile1_hq.sam and HiCfile2_hq.sam) are screened out and both files are converted to BAM format:
```unix
awk '{print $1}' HiCfile1_hq.sam | sort > readnames1.txt
awk '{print $1}' HiCfile2_hq.sam | sort > readnames2.txt
comm -12 readnames1.txt readnames2.txt > paired_reads.txt

grep -Fwf paired_reads.txt HiCfile1_hq.sam | \
cat header1.txt - | \
samtools view -b -@ 32 - > HiCfile_pair1.bam

grep -Fwf paired_reads.txt HiCfile2_hq.sam | \
cat header2.txt - | \
samtools view -b -@ 32 - > HiCfile_pair2.bam

rm HiCfile1_hq.sam
rm HiCfile2_hq.sam
```
- ```readnames1.txt``` contains the names of the high quality mapped reads in HiCfile1_hq.sam.
- ```readnames2.txt``` contains the names of the high quality mapped reads in HiCfile2_hq.sam.
- ```paired_reads.txt``` contains the names of the reads that have both high quality sides.
- ```HiCfile_pair1.bam``` is the bam file associated to the first mate.
- ```HiCfile_pair2.bam``` is the bam file associated to the second mate.

The -@ parameter is used to set the number of BAM compression threads to use in addition to the main thread (update it accordingly to your number of available cores). This can be useful to decrease the processing time.

**Both the bam files will serve as inputs for the normalization pipeline.**

To update the log files with pairing statistics, use the following code:
```unix
n=`wc -l paired_reads.txt | awk '{print $1}'`

ntot1=`wc -l readnames1.txt | awk '{print $1}'`
ntot2=`wc -l readnames2.txt | awk '{print $1}'`

perc1=$(awk -v n1=$n -v ntot1=$ntot1 'BEGIN { print 100*n1/ntot1 }' | cut -c1-5)
perc2=$(awk -v n2=$n -v ntot2=$ntot2 'BEGIN { print 100*n2/ntot2 }' | cut -c1-5)

printf "; of these:\n    "$n" ("$perc1"%%) were paired and saved into HiCfile_pair1.bam" >> HiCfile1_log.txt
printf "; of these:\n    "$n" ("$perc2"%%) were paired and saved into HiCfile_pair2.bam" >> HiCfile2_log.txt
```
The final and updated log files look like these:
```unix
HiCfile1_log.txt
202095066 reads; of these:
202095066 (100.00%) were unpaired; of these:
5770798 (2.86%) aligned 0 times
156759009 (77.57%) aligned exactly 1 time
39565259 (19.58%) aligned >1 times
97.14% overall alignment rate

----------
202095066 reads; of these:
172973813 (85.59%) aligned with MAPQ>=30; of these:
143415284 (82.91%) were paired and saved into HiCfile_pair1.bam

HiCfile2_log.txt
202095066 reads; of these:
202095066 (100.00%) were unpaired; of these:
13381441 (6.62%) aligned 0 times
149852422 (74.15%) aligned exactly 1 time
38861203 (19.23%) aligned >1 times
93.38% overall alignment rate

----------
202095066 reads; of these:
161438783 (79.88%) aligned with MAPQ>=30; of these:
143415284 (88.83%) were paired and saved into HiCfile_pair2.bam
```

## 5. Creating the fragment-end (FEND) bed file


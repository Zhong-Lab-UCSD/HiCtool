# Data preprocessing

This is the first section of the pipeline and it allows to pre-process the raw Hi-C data (fastq files), in order to generate input files for the normalization step. For more information about the Python functions used here check the [API documentation](https://sysbio.ucsd.edu/public/rcalandrelli/HiCtool_API_documentation.pdf).

## Table of Contents

1. [Preprocessing the data](#1-preprocessing-the-data)
   - [1.1. Downloading the raw data from GEO](#11-downloading-the-raw-data-from-geo)
2. [Creating the fragment-end (FEND) bed file](#2-creating-the-fragment-end-fend-bed-file)

## 1. Preprocessing the data

In order to start the Hi-C data analysis and preprocess your data you have to provide: 

- Two fastq files, respectively of the first and the second reads of the pairs. If you wish to use public datasets on GEO and you need instructions to download the data, see [section 1.1.](#11-downloading-the-raw-data-from-geo). 
- The Bowtie2 indexes of your reference genome.
- The restriction enzyme used in the Hi-C experiment.

**Note!** To reproduce the results of this tutorial, use the public dataset with this GEO accession number: **[GSM1551550](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1551550)**.

HiCtool allows to process Hi-C data generated for one of the following species:

- hg19
- hg38
- mm9
- mm10
- dm6
- susScr3
- susScr11

HiCtool allows to process Hi-C data generated using one or more of the following restriction enzymes:

- HindIII
- MboI
- DpnII
- Sau3AI
- BglII
- NcoI
- Hinfl

The **Arima Kit** uses a cocktail of restriction enzymes which includes MboI and Hinfl.
If your experiment was performed using a different species or restriction enzyme(s), please contact Riccardo Calandrelli at <rcalandrelli@eng.ucsd.edu>.

If you do not have the Bowtie2 genome index, please run the following in order to generate it:
```unix
bowtie2-build hg38.fa index
```
```hg38.fa``` is the reference genome sequence in FASTA format (in this case of hg38), the output files in ``.bt2`` format are named with the prefix ``index``.

**The data preprocessing is performed with a single unix command line (replace parameters in the code below properly) and comprises the following steps:**

1. Pre-truncation of the reads that contain potential ligation junctions to keep the longest piece without a junction sequence ([Ay et al., 2015](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0745-7)).
2. Independent mapping of the read pairs to the reference genome to avoid any proximity constraint.
3. Removing the unmapped reads and selecting reads that were uniquely mapped with a MAPQ >= 30, i.e. the estimated probability of mapping error is <= 0.1%.

```unix
# Make the bash script executable
chmod u+x /HiCtool-master/scripts/HiCtool_run_preprocessing.sh

# Run the script
/HiCtool-master/scripts/HiCtool_run_preprocessing.sh \
-h /HiCtool-master/scripts/ \
-o /your_output_directory/ \
-1 /myfastq_path/file1.fastq \
-2 /myfastq_path/file2.fastq \
-e MboI \
-g /path_to_the_genome_indexes/index \
-p 32
```
where:

- ``-h``: path to the HiCtool scripts with the final trailing slash ``/``.
- ``-o``: path to save the output files. If the folder does not exist, it is created automatically.
- ``-1``: the fastq file with the first reads of the pairs.
- ``-2``: the fastq file with the second reads of the pairs.
- ``-e``: the restriction enzyme or enzymes passed between square brackets (example: [MboI,Hinfl] for the cocktail of the Arima Kit).
- ``-g``: Bowtie2 genome indexes. Only the filename should be passed here without extension, in this case ``index``.
- ``-p``: the number of parallel threads (processors) to use for alignment and preprocessing. The more the fastest the process.

The structure of the output directory is the following:
```unix
~/your_output_directory
	|___ file1.trunc.fastq
	|___ file2.trunc.fastq
	|___ pre_truncation_log.txt
	|___ file1.fastq_truncated_reads.pdf
	|___ file2.fastq_truncated_reads.pdf
	|___ HiCfile1_pair1.bam
	|___ HiCfile2_pair1.bam
	|___ HiCfile1_log.txt
	|___ HiCfile2_log.txt
```

**The following output files are generated:**

- ``file1.trunc.fastq`` and ``file2.trunc.fastq`` are the two fastq file with pre-truncated reads.
- ``pre_truncation_log.txt`` with the information about the percentage of reads that have been truncated. This is also printed on the console during pre-processing:
```unix
SRR1658570_1.fastq
202095066 reads (length = 101 bp); of these:
29851195 (14.78%) contained a potential ligation junction and have been truncated.
SRR1658570_2.fastq
202095066 reads (length = 101 bp); of these:
28681691 (14.2%) contained a potential ligation junction and have been truncated.
```
- ``file1.fastq_truncated_reads.pdf`` and ``file2.fastq_truncated_reads.pdf`` are plots of the length distribution of the truncated reads, for both the fastq files:

![](/figures/SRR1658570_1.fastq_truncated_reads.png)

![](/figures/SRR1658570_2.fastq_truncated_reads.png)

- ``HiCfile_pair1.bam`` and ``HiCfile_pair2.bam`` that are the bam files of the pre-truncated first and second reads in the pairs respectively, generated after alignment and filtering.
- ``HiCfile1_log.txt`` and ``HiCfile2_log.txt`` are the log files with alignment and filtering statistics for the first and second reads in the pairs respectively.
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

### 1.1. Downloading the raw data from GEO

The source data in sra format are downloaded via GEO accession number using the command ``fastq-dump`` of [SRA Toolkit](http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc&f=fastq-dump).

Before proceeding, you may need to [setup the output directory](https://github.com/ncbi/sra-tools/wiki/Toolkit-Configuration) where the sra files will be saved. After having installed SRA Toolkit, go to the path where the software has been installed, under the subfolder ``/bin``, and run the following command line:
```unix
./vdb-config -i
```
This will open an interface that will allow you to setup/change your output directory.

To download the data related to a GEO accession number, go to the bottom of that page and click on the SRA number under the section “Relations”. After that, under the section “Runs” you will find the SRR files, then run the following:
```unix
fastq-dump SRRXXXXXXX --split-3
```
where ``SRRXXXXXXX`` has to be replaced with the specific number of the run you want to download (**[SRR1658570](https://www.ncbi.nlm.nih.gov/sra?term=SRX764936) in this documentation**).

To be more specific, this code will either download the SRA file under the output directory you have set with the interface above, but it will also convert the SRA file into fastq files and dump each read into separate files in the current working directory, to produce:

- ``SRRXXXXXXX_1.fastq``
- ``SRRXXXXXXX_2.fastq``
- ``SRRXXXXXXX.fastq``

where paired-end reads in ``SRRXXXXXXX.sra`` are split and stored into **``SRRXXXXXXX_1.fastq``** and **``SRRXXXXXXX_2.fastq``**, and ``SRRXXXXXXX.fastq`` (if present) contains reads with no mates.

**Note!** To produce our final results, use this GEO accession number: **[GSM1551550](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1551550)**.


## 2. Creating the fragment-end (FEND) bed file

The fragment-end (FEND) bed file is used to normalize the data and it contains restriction site coordinates and optional additional information related to fragment properties, such as GC content and mappability score (needed only if you wish to normalize your data using the [explicit-factor approach from Yaffe and Tanay](/tutorial/normalization-yaffe-tanay.md). Specifically, for each fragment the GC content of 200 bp upstream and downstream to the restriction site is computed. For the mappability score, the entire genome sequence is split into artificial reads (50 bp reads, starting every 10 bp) and then mapped back to the genome. For each fragment end the mappability score is then defined to be the portion of artificial reads mapped uniquely to the genome (MAPQ > 30) within a 500-bp window upstream and downstream to the fragment. Fragment ends with a mappability score less than 0.5 are then discarded (Yaffe and Tanay, 2011).

Since generating the FEND bed file may be time consuming, we provide the most common FEND files available for download (DpnII is the same restriction site than MboI):

- [HindIII-hg38](http://data.genomegitar.org/HindIII_hg38_gc_map_valid.zip)
- [MboI-hg38 (or DpnII-hg38)](http://data.genomegitar.org/MboI_hg38_gc_map_valid.zip) (file used in this documentation)
- [Arima_Kit-hg38](http://data.genomegitar.org/arima_gc_map.zip) (restriction ezymes cocktail with MboI and Hinfl)
- [NcoI-hg38](http://data.genomegitar.org/NcoI_hg38_gc_map_valid.zip)
- [HindIII-mm10](http://data.genomegitar.org/HindIII_mm10_gc_map_valid.zip)
- [MboI-mm10 (or DpnII-mm10)](http://data.genomegitar.org/MboI_mm10_gc_map_valid.zip)

**Perform the following steps ONLY if you need to generate a new fragment end bed file (because you are using another species or a different restriction enzyme than those provided above). Otherwise, download the file of interest above and go to the [data normalization section](https://github.com/Zhong-Lab-UCSD/HiCtool/tree/master/tutorial#2-data-normalization-and-visualization).**

***

**The generation of the FEND file is performed with a single unix command line (replace parameters in the code below accordingly) and comprises the following steps:**

1. Generation of the restriction enzyme sequence.
2. Multiple alignment of the restriction site using Bowtie 2, to locate all the coordinates of the restriction enzyme sites.
3. Conversion of the output sam file to bed format.

```unix
# Make the bash script executable
chmod u+x /HiCtool-master/scripts/HiCtool_generate_fend_file.sh

# Run the script
/HiCtool-master/scripts/HiCtool_generate_fend_file.sh \
-h /HiCtool-master/scripts/ \
-o your_output_directory \
-e MboI \
-g /path_to_the_genome_indexes/index \
-s hg38 \
-p 32
```
where:

- ``-h``: path to the HiCtool scripts with the final trailing slash ``/``.
- ``-o``: path to save the output files. If the folder does not exist, it is created automatically.
- ``-e``: the restriction enzyme or enzymes passed between square brackets (example: [MboI,Hinfl] for the cocktail of the Arima Kit).
- ``-g``: Bowtie2 genome indexes. Only the filename should be passed here without extension, in this case ``index``.
- ``-s``: species name.
- ``-p``: the number of parallel threads (processors) to use for alignment and preprocessing. The more the fastest the process.

The structure of the output directory is the following:
```unix
~/your_output_directory
	|___ restrictionsites.bed
	|___ restrictionsites_log.txt
```

**The following output files are generated:**

- ``restrictionsites.bed`` is the FEND bed file.
- ``restrictionsites_log.txt`` is the log file with the restriction site multiple alignment statistics.

- **If you will be using the [Hi-Corrector normalization approach](https://github.com/Zhong-Lab-UCSD/HiCtool/blob/master/tutorial/normalization-matrix-balancing.md), you do NOT need to perform the following steps, because additional information such as GC content or mappability score are not needed. You can use ``restrictionsites.bed`` as the FEND file.**
- **Proceed with the following section ONLY if you will be using [Yaffe and Tanay's normalization approach](https://github.com/Zhong-Lab-UCSD/HiCtool/blob/master/tutorial/normalization-yaffe-tanay.md).**

***

You can decide if adding both the GC content and mappability score information or only one of them. Our suggestion is to add both at this point, so you will have a complete FEND file, then you may decide later on if not considering one of these features for normalization. This process is very time consuming, so we suggest to use the highest number of processors available. Using 24 threads, we took around 18 hours to add the GC content and mappability score information for hg38-MboI. This time may exponentially increase if you are using a mixture of restriction sites in your experiment: for the ARIMA kit cocktail (over 23 millions restriction sites for hg38) we took around ...

**The update of the FEND file to add GC content and/or mappability score is performed with a single unix command line (replace parameters in the code below accordingly) and comprises the following steps:**

1. Splitting ``restrictionsites.bed`` into single bed files, one per each chromosome.
2. Adding GC content information per each single chromosome (if GC content selected).
3. Adding mappability score information per each single chromosome (if mappability selected).
4. Merging all the files, removing fragments with mappability score < 0.5, and parsing the file format to generate the final FEND bed file.

```unix
# Make the bash script executable
chmod u+x /HiCtool-master/scripts/HiCtool_fend_file_add_gc_map.sh

# Run the script
/HiCtool-master/scripts/HiCtool_fend_file_add_gc_map.sh \
-h /HiCtool-master/scripts/ \
-o your_output_directory \
-s hg38 \
-p 32 \
-b /a_path/gc5Base.bw \
-m /a_path/hg38.fa
```
where:

- ``-h``: path to the HiCtool scripts with the final trailing slash ``/``.
- ``-o``: path to save the output files. If the folder does not exist, it is created automatically.
- ``-s``: species name.
- ``-p``: the number of parallel threads (processors) to use for alignment and preprocessing. The more the fastest the process.
- ``-b``: gc5Base.bw file with the GC content information (download it from this website: hgdownload.cse.ucsc.edu/gbdb/your_species/bbi/). **If the GC content information has been already calculated for the reference genome** (i.e. you have run already this code using the same reference genome), you will have a ``GC_info`` directory with several ``txt`` files (one per each chromosome); therefore you can input this directory path to avoid this redundant step. If ``-b`` is not declared, the GC content information is not added.
- ``-m``: the reference genome in fasta format to add the mappability score to the FEND file. **If the mappability score has been already calculated for the reference genome** (i.e. you have run already this code using the same reference genome), you will have a ``mappability_info`` directory with several ``txt`` files (one per each chromosome); therefore you can input this directory path to avoid this redundant step. If ``-m`` is not declared, the mappability score information is not added.

If both ``-b`` and ``-m`` are declared, the structure of the output directory is the following (``i`` stands for a general chromosome, and there are gonna be one of those files per each chromosome):
```unix
~/your_output_directory
	|___ restrictionsites.bed
	|___ chri.bed
	|___ chri_gc.bed
	|___ chri_map.bed
	|___ restrictionsites_gc_map.bed
	|___ GC_info
		|___ chri.txt
	|___ mappability_info
		|___ chri.txt
		|___ artificial_reads.log
```
**The following output files are generated:**

- ``restrictionsites.bed`` is the FEND bed file without information of GC content or mappability score.
- ``chri.bed`` are single bed files, one per each chromosome, derived by the splitting of ``restrictionsites.bed``.
- ``chri_gc.bed`` are the correspondent of ``chri.bed`` with GC content information added.
- ``chri_map.bed`` are the correspondent of ``chri.bed`` with mappability score information added.
- ``restrictionsites_gc_map.bed`` is the final FEND file with GC content and mappability score.
- ``GC_info`` is a folder which contains the GC content information, one file per chromosome, specific to the reference genome. If you wish to compute another FEND file for the same reference genome but different restriction enzyme, you may keep this folder and input it later in ``-b``.
- ``mappability_info`` is a folder which contains the mappability score information, one file per chromosome, specific to the reference genome. If you wish to compute another FEND file for the same reference genome but different restriction enzyme, you may keep this folder and input it later in ``-m``. This folder contains also ``artificial_reads.log`` which contains the alignment statistics of the artificial reads used to compute the mappability score of the fragments.

 **``restrictionsites_gc_map.bed`` is the FEND file that you will use in [Yaffe and Tanay's normalization approach](https://github.com/Zhong-Lab-UCSD/HiCtool/blob/master/tutorial/normalization-yaffe-tanay.md).**

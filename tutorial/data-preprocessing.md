# Data preprocessing

This is the first section of the pipeline and it allows to pre-process the raw Hi-C data (fastq files), in order to generate input files for the normalization step.

## Table of Contents

1. [Downloading the raw data from GEO](#1-downloading-the-raw-data-from-geo)
2. [Pre-truncation of the reads that contain potential ligation junctions](#2-pre-trunction-of-the-reads-that-contain-potential-ligation-junction)
3. [Mapping read pairs to the reference genome](#3-mapping-read-pairs-to-reference-genome)
4. [Filtering reads and selecting reads that are paired](#4-filtering-reads-and-selecting-reads-that-are-paired)
5. [Creating the fragment-end (FEND) bed file](#5-creating-the-fragment-end-fend-bed-file)

## 1. Downloading the raw data from GEO

**Note!** If you have your fastq files generated from your custom experiment, you can skip this first step and go to [step 2](#2-pre-trunction-of-the-reads-that-contain-potential-ligation-junction).

The source data in sra format are downloaded via GEO accession number using the command ``fastq-dump`` of [SRA Toolkit](http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc&f=fastq-dump).

Before proceeding, you may need to [setup the output directory](https://github.com/ncbi/sra-tools/wiki/Toolkit-Configuration) where the sra files will be saved. After having installed SRA Toolkit, go to the path where the software has been installed, under the subfolder “bin”, and run the following command line:
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

- SRRXXXXXXX_1.fastq
- SRRXXXXXXX_2.fastq
- SRRXXXXXXX.fastq

where paired-end reads in SRRXXXXXXX.sra are split and stored into **SRRXXXXXXX_1.fastq** and **SRRXXXXXXX_2.fastq**, and SRRXXXXXXX.fastq (if present) contains reads with no mates.

**Note!** To produce our final results, use this GEO accession number: **[GSM1551550](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1551550)**.

## 2. Pre-trunction of the reads that contain potential ligation junction

After the fastq files are obtained, pre-truncation is performed on the reads that contain potential ligation junctions to keep the longest piece without a junction sequence ([Ay et al., 2015](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0745-7)). To do so, use the code in [HiCtool_pre_truncation.py](/scripts/HiCtool_pre_truncation.py) and run the following code on your Python or iPython console:
```Python
execfile('HiCtool_pre_truncation.py')
pre_truncation('SRR1658570_1.fastq', 'MboI')
pre_truncation('SRR1658570_2.fastq', 'MboI')
```
where the first argument is the fastq file, the second argument is the restriction enzyme (see API documentation).

The output files will have the same filename with extension ```.trunc.fastq```. If a different restriction enzyme (RE) than HindIII, MboI, NcoI and DpnII has been used in the Hi-C experiment, then run the following for each fastq file:
```Python
execfile('HiCtool_pre_truncation.py')
pre_truncation('a_file.fastq', 'custom_RE', 'custom_ligation_junction')
```
where ``custom_ligation_junction`` is the ligation junction sequence of nucleotides associated to the restriction enzyme you are using (see [Ay et al., 2015](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0745-7) for more details), ``custom_RE`` is your restriction enzyme name.

A log file named **pre_truncation_log.txt** is generated with the information about the percentage of reads that have been truncated. This is also printed on the console:
```unix
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

[Bowtie 2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) is used for mapping the read pairs, and reads are **mapped independently** to avoid any proximity constraint. To align the reads, first build the index for the reference genome (execute this step only once per reference genome):
```unix
bowtie2-build hg38.fa index
```
```hg38.fa``` is the reference sequence in FASTA format, the output files in ``bt2`` format are named with the prefix ``index``.

Now align the reads to the reference sequence:
```unix
(bowtie2 -p 32 -x index SRR1658570_1.trunc.fastq -S HiCfile1.sam) 2>HiCfile1_log.txt
(bowtie2 -p 32 -x index SRR1658570_2.trunc.fastq -S HiCfile2.sam) 2>HiCfile2_log.txt
```
where:

- ``-p`` refers to a specified number of parallel search threads (update it accordingly to your number of available cores). This can be useful to decrease the processing time in aligning the reads.
- ``-x`` specifies the basename of the index for the reference genome. The basename is the name of any of the index files up to but not including the final ``/.1.bt2``, ``/.2.bt2``, etc.
- ``-S`` specifies the output file in ``sam`` format.
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

[SAMtools](http://samtools.sourceforge.net/) is used to perform filtering on the mapped reads ([Yaffe and Tanay, 2011](http://www.nature.com/ng/journal/v43/n11/abs/ng.947.html)):
```unix
# extracting the headers
samtools view -H HiCfile1.sam > header1.txt
samtools view -H HiCfile2.sam > header2.txt

# read filtering
samtools view -F 4 -q 30 HiCfile1.sam > HiCfile1_hq.sam
samtools view -F 4 -q 30 HiCfile2.sam > HiCfile2_hq.sam
```
- ``-F 4`` is used to discard unmapped reads.
- ``-q 30`` is used to select reads that were uniquely mapped with a MAPQ >= 30, i.e. the estimated probability of mapping error is <= 0.1%.

**Unused files are removed** using the command ```rm``` to reduce memory occupance since they may be quite big. If you want to keep your files, do not run those command lines.

To **update the log files** from mapping with information about reads mapped with MAPQ >= 30, use the following code:
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
After filtering, reads that are not paired from the two SAM files (``HiCfile1_hq.sam`` and ``HiCfile2_hq.sam``) are screened out and both files are converted to BAM format:
```unix
awk '{print $1}' HiCfile1_hq.sam | sort > readnames1.txt
awk '{print $1}' HiCfile2_hq.sam | sort > readnames2.txt
comm -12 readnames1.txt readnames2.txt > paired_reads.txt

# Select reads that are paires with the second sam file
grep -Fwf paired_reads.txt HiCfile1_hq.sam | \
cat header1.txt - | \
samtools view -b -@ 32 - > HiCfile_pair1.bam

# Select reads that are paired with the first sam file
grep -Fwf paired_reads.txt HiCfile2_hq.sam | \
cat header2.txt - | \
samtools view -b -@ 32 - > HiCfile_pair2.bam

rm HiCfile1_hq.sam
rm HiCfile2_hq.sam
```
- ``readnames1.txt`` contains the names of the high quality mapped reads in ``HiCfile1_hq.sam``.
- ``readnames2.txt`` contains the names of the high quality mapped reads in ``HiCfile2_hq.sam``.
- ``paired_reads.txt`` contains the names of the reads that have both high quality sides.
- ``HiCfile_pair1.bam`` is the BAM file associated to the first mate.
- ``HiCfile_pair2.bam`` is the BAM file associated to the second mate.
- ``-@`` is used to set the number of BAM compression threads to use in addition to the main thread (update it accordingly to your number of available cores). This can be useful to decrease the processing time.

**Both the bam files will serve as inputs for the normalization pipeline.**

To **update the log files** with pairing statistics, use the following code:
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

The fragment-end (FEND) bed file is used to normalize the data and it contains restriction site coordinates and additional information related to fragment properties (GC content and mappability score). Specifically, for each fragment the GC content of 200 bp upstream and downstream to the restriction site is computed. For the mappability score, the entire genome sequence is split into artificial reads (50 bp reads, starting every 10 bp) and then mapped back to the genome. For each fragment end the mappability score is then defined to be the portion of artificial reads mapped uniquely to the genome (MAPQ > 30) within a 500-bp window upstream and downstream to the fragment. Fragment ends with a mappability score less than 0.5 are then discarded (Yaffe and Tanay, 2011).

Since the following steps may be time consuming, we provide the most common FEND files available for download (DpnII is the same restriction site than MboI):

- [HindIII-hg38](http://data.genomegitar.org/HindIII_hg38_gc_map_valid.zip)
- [MboI-hg38](http://data.genomegitar.org/MboI_hg38_gc_map_valid.zip) (file used in this documentation)
- [NcoI-hg38](http://data.genomegitar.org/NcoI_hg38_gc_map_valid.zip)
- [HindIII-mm10](http://data.genomegitar.org/HindIII_mm10_gc_map_valid.zip)
- [MboI-mm10](http://data.genomegitar.org/MboI_mm10_gc_map_valid.zip)

**Perform the following steps ONLY if you need to generate a new fragment end bed file (because you are using another species or a different restriction enzyme than those provided above). Otherwise, download the file of interest and go to the [data normalization section](https://github.com/Zhong-Lab-UCSD/HiCtool/tree/master/tutorial#2-data-normalization-and-visualization).**

***

In order to align all the restriction sites for a certain cutting enzyme, a ``fastq`` file related to the enzyme cutting site has to be provided. For the quality score of the restriction enzyme sequence, we can simply add a default average score ``I``:
```unix
echo -e "@HindIII\nAAGCTT\n+\nIIIIII" > HindIII.fastq
echo -e "@MboI\nGATC\n+\nIIII" > MboI.fastq
echo -e "@NcoI\nCCATGG\n+\nIIIIII" > NcoI.fastq
```
After this, implement the multiple alignment command in Bowtie 2 to locate all the coordinates of the restriction enzyme sites:
```unix
(bowtie2 -p 32 -k 8000000 -x your_genome_index -U MboI.fastq -S restrictionsites.sam) 2>restrictionsites_log.txt
```
where the ``-k`` argument changes Bowtie 2 research behavior. By default, ``bowtie2`` searches for distinct, valid alignments for each read. When it finds a valid alignment, it continues looking for alignments that are nearly as good or better and the best alignment found is reported. When ``-k <int>`` is specified, ``bowtie2`` searches for at most ``<int>`` distinct, valid alignments for each read. The search terminates when it can not find more distinct valid alignments, or when it finds ``<int>``, which happens first.

In order to use the restriction sites file as an input for the following analysis, we have to convert the ``sam`` file to ``bed`` file via [SAMtools](http://samtools.sourceforge.net/) and [bedtools](http://bedtools.readthedocs.org/en/latest/):
```unix
samtools view -b restrictionsites.sam | bedtools bamtobed -i > restrictionsites.bed
```
- **If you will be using the [Hi-Corrector normalization approach](https://github.com/Zhong-Lab-UCSD/HiCtool/blob/master/tutorial/normalization-matrix-balancing.md), you do NOT need to perform the following steps, because additional information such as GC content or mappability score are not needed. You can use ``restrictionsites.bed`` as the FEND file.**
- **Proceed with the following steps if you will be using [Yaffe and Tanay's normalization approach](https://github.com/Zhong-Lab-UCSD/HiCtool/blob/master/tutorial/normalization-yaffe-tanay.md).**

***

First ``restrictionsites.bed`` is split into separate files, one per each chromosome (**update the chromosomes list according to your species**):
```unix
chromosomes=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrX" "chrY")

for i in "${chromosomes[@]}"; do
awk -v var="$i" '(NR>1) && ($1==var)' restrictionsites.bed > $i.bed
done
```

The GC content and mappability score information must be in comma-separated format. First, each feature is computed for each separate chromosome, finally the files are merged together and parsed to generate a unique bed file. For this part, **Python multiprocessing** is used to consistently reduce the computation time. It is recommended to use the highest number of threads available in your processor (the highest up to the total number of chromosomes of your species).

- Download the GC content information for the species of interest from the UCSC website at http://hgdownload.cse.ucsc.edu/gbdb/hg38/bbi/ (replace the species name in the link if needed) and use the file ``gc5Base.bw``. This file is in BigWig format, before running step 2) we need to convert it to BedGraph format (tab separated file with 4 columns: chromosome, start, end, score; in this case "score" is the GC content). After this, the file has to be splitted into separate txt files, one per each chromosome, named as **chr1.txt, chr2.txt, … , chrX.txt, chrY.txt**. To do so, run the following Unix script (note that given the dimension of the files, this process may require sometime):
```unix
wget link_to_gc5Base.bw
bigWigToBedGraph gc5Base.bw gc5Base.bedGraph

chromosomes=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrX" "chrY")
for i in "${chromosomes[@]}"; do
awk -v var="$i" '(NR>1) && ($1==var)' gc5Base.bedGraph | awk -v OFS='\t' '{print $1, $2, $3, $4}' > $i.txt
done
```
- Add the GC content information using the Python script [HiCtool_add_fend_gc_content.py](/scripts/HiCtool_add_fend_gc_content.py). Open the script, **update the parameters on the top and save**. Then just execute the script to add the gc content information (using 24 threads, we took around 9 hours for all the chromosomes of hg38-MboI):
```Python
execfile('HiCtool_add_fend_gc_content.py')
```
- Generate artificial reads using the Python function inside [HiCtool_artificial_reads.py](/scripts/HiCtool_artificial_reads.py) (this step is required only once per reference genome):
```Python
execfile('HiCtool_artificial_reads.py')
generate_artificial_reads(genome_file, output_reads_file)
```
where ``genome_file`` is the reference genome sequence in ``fasta`` format, ``output_reads_file`` is the file where to save reads in ``fastq`` format.
- Align the artificial reads to the reference genome, remove unmapped reads and generate a bed formatted file where the last field is the MAPQ score (this step is required only once per reference genome). This is the same format of the GC content files downloaded from UCSC, where the last field was the GC percentage instead of MAPQ:
```unix
(bowtie2 -p 32 -x your_genome_index artificial_reads.fastq -S artificial_reads.sam) 2>artificial_reads.log
samtools view -F 4 artificial_reads.sam > artificial_reads_mapped.sam
awk -v OFS='\t' '{print $3, $4-1, $4-1+50, $5}' artificial_reads_mapped.sam > artificial_reads_mapped.txt
```
- Split the mapped reads into separate files, one per chromosome (**update the list of chromosome names if a different species is used**):
```unix
chromosomes=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrX" "chrY")
for i in "${chromosomes[@]}"; do
awk -v var="$i" '(NR>1) && ($1==var)' artificial_reads_mapped.txt | awk -v OFS='\t' '{print $1, $2, $3, $4}' > $i.txt
done
```
- Add the mappability score information using the Python script [HiCtool_add_fend_mappability.py](/scripts/HiCtool_add_fend_mappability.py). Open the script, **update the parameters on the top and save**. Then just execute the script to add the mappability information (using 24 threads, we took around 9 hours for all the chromosomes of hg38-MboI):
```Python
execfile('HiCtool_add_fend_mappability.py')
```
- Sort by coordinate, merge the files together, remove fragment ends with a mappability score < 0.5, parse GC content and mappability score in comma separated format and add the header:
```unix
chromosomes=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrX" "chrY")

for i in "${chromosomes[@]}"; do
sort -k 2,2n "$i"_restrictionsites_gc_map.bed | cat >> restrictionsites_gc_map.bed
done

awk '(NR>1) && ($9 >= 0.5) && ($10 >= 0.5)' restrictionsites_gc_map.bed | awk -v OFS='\t' '{print $1, $2, $3, $4, $5, $6, $7 "," $8, $9 "," $10}' > restrictionsites_gc_map_valid_noHeader.bed
echo -e 'chr\tstart\tstop\tname\tscore\tstrand\tgc\tmappability' | cat - restrictionsites_gc_map_valid_noHeader.bed > restrictionsites_gc_map_valid.bed

rm restrictionsites_gc_map_valid_noHeader.bed
```
**restrictionsites_gc_map_valid.bed** is the final FEND bed file that will be used in the normalization pipeline to remove biases.

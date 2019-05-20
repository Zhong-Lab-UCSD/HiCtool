while getopts h:o:1:2:e:g:p: option
do
case "${option}"
in
h) hictoolPath=${OPTARG};; # The path where are the HiCtool scripts with the final trailing slash.
o) outputPath=${OPTARG};; # The path where to save the output files.
1) fastq1=${OPTARG};; # The fastq file with the first reads of the pairs.
2) fastq2=${OPTARG};; # The fastq file with the second reads of the pairs.
e) restrictionEnzyme=${OPTARG};; # The restriction enzyme or enzymes passed between square brackets (example: [enzyme1,enzyme2]).
g) genomeIndex=${OPTARG};; # The Bowtie2 genome indexes of the reference genome (only filename without extension).
p) threads=${OPTARG};; # The number of parallel threads to use for alignment and pre-processing. The more the fastest the process.
esac
done

echo "Start data preprocessing: $(date)"

checkMakeDirectory(){
if [ ! -e "$1" ]; then
mkdir -p "$1"
fi
}
checkMakeDirectory $outputPath
cd $outputPath

python $hictoolPath"HiCtool_pre_truncation.py" -i [$fastq1,$fastq2] -e $restrictionEnzyme

fastq1_trunc="${fastq1%%.*}.trunc.fastq"
fastq2_trunc="${fastq2%%.*}.trunc.fastq"

echo -n "Aligning "$fastq1_trunc" ... "
(bowtie2 -p $threads -x $genomeIndex $fastq1_trunc -S HiCfile1.sam) 2>HiCfile1_log.txt
echo "Done!"
echo -n "Aligning "$fastq2_trunc" ... "
(bowtie2 -p $threads -x $genomeIndex $fastq2_trunc -S HiCfile2.sam) 2>HiCfile2_log.txt
echo "Done!"

# extracting the headers and read filtering
echo -n "Filtering HiCfile1.sam ... "
samtools view -H HiCfile1.sam > header1.txt
samtools view -F 4 -q 30 HiCfile1.sam > HiCfile1_hq.sam
echo "Done!"
echo -n "Filtering HiCfile2.sam ... "
samtools view -H HiCfile2.sam > header2.txt
samtools view -F 4 -q 30 HiCfile2.sam > HiCfile2_hq.sam
echo "Done!"

echo -n "Building log files ... "
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
echo "Done!"

rm HiCfile1.sam
rm HiCfile2.sam

echo -n "Selecting reads that are paired ... "
awk '{print $1}' HiCfile1_hq.sam | sort > readnames1.txt
awk '{print $1}' HiCfile2_hq.sam | sort > readnames2.txt
comm -12 readnames1.txt readnames2.txt > paired_reads.txt

# Select reads that are paires with the second sam file
grep -Fwf paired_reads.txt HiCfile1_hq.sam | \
cat header1.txt - | \
samtools view -b -@ $threads - > HiCfile_pair1.bam
rm HiCfile1_hq.sam

# Select reads that are paired with the first sam file
grep -Fwf paired_reads.txt HiCfile2_hq.sam | \
cat header2.txt - | \
samtools view -b -@ $threads - > HiCfile_pair2.bam
rm HiCfile2_hq.sam

echo "Done!"

echo -n "Updating log files ... "
n=`wc -l paired_reads.txt | awk '{print $1}'`

ntot1=`wc -l readnames1.txt | awk '{print $1}'`
ntot2=`wc -l readnames2.txt | awk '{print $1}'`

perc1=$(awk -v n1=$n -v ntot1=$ntot1 'BEGIN { print 100*n1/ntot1 }' | cut -c1-5)
perc2=$(awk -v n2=$n -v ntot2=$ntot2 'BEGIN { print 100*n2/ntot2 }' | cut -c1-5)

printf "; of these:\n    "$n" ("$perc1"%%) were paired and saved into HiCfile_pair1.bam" >> HiCfile1_log.txt
printf "; of these:\n    "$n" ("$perc2"%%) were paired and saved into HiCfile_pair2.bam" >> HiCfile2_log.txt
echo "Done!"

rm header1.txt
rm header2.txt
rm readnames1.txt
rm readnames2.txt
rm paired_reads.txt

echo "End data preprocessing: $(date)"

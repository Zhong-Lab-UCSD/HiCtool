while getopts h:o:e:g:s:p:b:m: option
do
case "${option}"
in
h) hictoolPath=${OPTARG};; # The path of the HiCtool scripts
o) directory=${OPTARG};; # The output directory to work on and save the fend file
e) restrictionEnzyme=${OPTARG};; # The restriction enzyme (or enzymes passed in the form of a Python list)
g) genomeIndex=${OPTARG};; # The Bowtie2 genome indexes of the reference genome
s) species=${OPTARG};; # Species under analysis: hg19, hg38, mm9, mm10, dm6, susScr2 are available
p) threads=${OPTARG};; # The number of threads to use
esac
done

if [ $species = 'hg38' -o $species = 'hg19' ]
then
	chromosomes=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrX" "chrY")
elif [ $species = 'mm10' -o $species = 'mm9' ]
then
	chromosomes=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chrX" "chrY")
elif [ $species = 'dm6' ]
then
	chromosomes=("chr2L" "chr2R" "chr3L" "chr3R" "chr4" "chrX" "chrY")
elif [ $species = 'susScr3' -o $species = 'susScr11' ]
then
	chromosomes=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chrX" "chrY")
else
	echo "ERROR! Wrong species inserted! Check the species spelling or insert an available species: hg38, hg19, mm10, mm9, dm6, susScr3, susScr11. If your species is not listed, please contact Riccardo Calandrelli at <rcalandrelli@eng.ucsd.edu>."
	exit
fi

### Create output directory
checkMakeDirectory(){
	if [ ! -e "$1" ]; then
		mkdir -p "$1"
	fi
}
checkMakeDirectory $directory

echo "Start fend file generation: $(date)"
cd $directory

python2.7 $hictoolPath"HiCtool_build_enzyme_fastq.py" -e $restrictionEnzyme

echo -n "Aligning restriction sites for "$restrictionEnzyme" ... "
(bowtie2 -p $threads -k 10000000 -x $genomeIndex -U restriction_enzyme.fastq -S restrictionsites.sam) 2>restrictionsites_log.txt
echo "Done!"

echo -n "Converting sam file to bed file ... "
samtools view -b -@ $threads restrictionsites.sam | bedtools bamtobed -i > restrictionsites.bed
rm restrictionsites.sam
echo "Done!"

echo "End fend file generation: $(date)"
echo "Your fend file is restrictionsites.bed"

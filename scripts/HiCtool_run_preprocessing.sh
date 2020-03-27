while getopts h:o:1:2:e:q:g:p:c: option
do
case "${option}"
in
h) hictoolPath=${OPTARG};; # The path where are the HiCtool scripts with the final trailing slash.
o) outputPath=${OPTARG};; # The path where to save the output files.
1) fastq1=${OPTARG};; # The fastq file with the first reads of the pairs.
2) fastq2=${OPTARG};; # The fastq file with the second reads of the pairs.
e) restrictionEnzyme=${OPTARG};; # The restriction enzyme or enzymes passed between square brackets (example: [enzyme1,enzyme2]).
q) quality=${OPTARG};; # The MAPQ value to filter mapped reads from SAMTOOLS
g) genomeIndex=${OPTARG};; # The Bowtie2 genome indexes of the reference genome (only filename without extension).
p) threads=${OPTARG};; # The number of parallel threads to use for alignment and pre-processing. The more the fastest the process.
c) chunk_size=${OPTARG};; # The number of lines per each temporary fastq file in order to avoid memory errors and multiprocessing to speed up the process. Each temporary file is processed by a separate processor if multiple threads are used.
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

echo -n "Calculating total lines of the fastq files ... "
fastq_lines=`wc -l $fastq1 | awk '{print $1}'`
echo "Done!"

if [ -z $chunk_size ]
then
	echo "chunk_size not declared."
	python2.7 $hictoolPath"HiCtool_pre_truncation.py" -i [$fastq1,$fastq2] -e $restrictionEnzyme -p $threads

	tot_reads_1=$(awk -F'\t' '{sum+=$1;} END{print sum;}' "${fastq1%%.*}_log.txt")
	tot_reads_2=$(awk -F'\t' '{sum+=$1;} END{print sum;}' "${fastq2%%.*}_log.txt")

	tot_trunc_reads_1=$(awk -F'\t' '{sum+=$3;} END{print sum;}' "${fastq1%%.*}_log.txt")
	tot_trunc_reads_2=$(awk -F'\t' '{sum+=$3;} END{print sum;}' "${fastq2%%.*}_log.txt")

	perc1=$(awk -v n1=$tot_trunc_reads_1 -v ntot1=$tot_reads_1 'BEGIN { print 100*n1/ntot1 }' | cut -c1-5)
	perc2=$(awk -v n1=$tot_trunc_reads_2 -v ntot1=$tot_reads_2 'BEGIN { print 100*n1/ntot1 }' | cut -c1-5)

	read_length=$(awk 'FNR == 1 {print $2}' "${fastq1%%.*}_log.txt")

	printf $fastq1"\n"$tot_reads_1" reads (length = "$read_length" bp); of these:\n  "$tot_trunc_reads_1" ("$perc1"%%)\n\n" > pre_truncation_log.txt
	printf $fastq2"\n"$tot_reads_2" reads (length = "$read_length" bp); of these:\n  "$tot_trunc_reads_2" ("$perc2"%%)" >> pre_truncation_log.txt

	rm "${fastq1%%.*}_log.txt"
	rm "${fastq2%%.*}_log.txt"

elif ! [ -z $chunk_size ] && [ $chunk_size -ge $fastq_lines ]
then
	echo "chunk_size not consider because greater that the total lines of the fastq file."
	python2.7 $hictoolPath"HiCtool_pre_truncation.py" -i [$fastq1,$fastq2] -e $restrictionEnzyme -p $threads

	tot_reads_1=$(awk -F'\t' '{sum+=$1;} END{print sum;}' "${fastq1%%.*}_log.txt")
	tot_reads_2=$(awk -F'\t' '{sum+=$1;} END{print sum;}' "${fastq2%%.*}_log.txt")

	tot_trunc_reads_1=$(awk -F'\t' '{sum+=$3;} END{print sum;}' "${fastq1%%.*}_log.txt")
	tot_trunc_reads_2=$(awk -F'\t' '{sum+=$3;} END{print sum;}' "${fastq2%%.*}_log.txt")

	perc1=$(awk -v n1=$tot_trunc_reads_1 -v ntot1=$tot_reads_1 'BEGIN { print 100*n1/ntot1 }' | cut -c1-5)
	perc2=$(awk -v n1=$tot_trunc_reads_2 -v ntot1=$tot_reads_2 'BEGIN { print 100*n1/ntot1 }' | cut -c1-5)

	read_length=$(awk 'FNR == 1 {print $2}' "${fastq1%%.*}_log.txt")

	printf $fastq1"\n"$tot_reads_1" reads (length = "$read_length" bp); of these:\n  "$tot_trunc_reads_1" ("$perc1"%%)\n\n" > pre_truncation_log.txt
	printf $fastq2"\n"$tot_reads_2" reads (length = "$read_length" bp); of these:\n  "$tot_trunc_reads_2" ("$perc2"%%)" >> pre_truncation_log.txt

	rm "${fastq1%%.*}_log.txt"
	rm "${fastq2%%.*}_log.txt"

elif ! [ -z $chunk_size ] && [ $chunk_size -lt $fastq_lines ]
then
	echo -n "Using chunk_size to split the fastq files ... "
	if (( $chunk_size % 4 )) ; then
		chunk_size=`expr $chunk_size - $(($chunk_size % 4))`
	fi
	# Splitting the first fastq file
	k=$chunk_size
	count=1
	while [ $k -lt $fastq_lines ]
	do
		start=`expr $k - $chunk_size + 1`
		quit=`expr $k + 1`
		sed -n "$start,"$k"p;"$quit"q" $fastq1 > "${fastq1%%.*}_temp_"$count".fastq"
		count=`expr $count + 1`
		k=`expr $k + $chunk_size`
	done
	start=`expr $k - $chunk_size + 1`
	sed -n "$start,"$fastq_lines"p" $fastq1 > "${fastq1%%.*}_temp_"$count".fastq"

	# Splitting the second fastq file
	k=$chunk_size
	count=1
	while [ $k -lt $fastq_lines ]
	do
		start=`expr $k - $chunk_size + 1`
		quit=`expr $k + 1`
		sed -n "$start,"$k"p;"$quit"q" $fastq2 > "${fastq2%%.*}_temp_"$count".fastq"
		count=`expr $count + 1`
		k=`expr $k + $chunk_size`
	done
	start=`expr $k - $chunk_size + 1`
	sed -n "$start,"$fastq_lines"p" $fastq2 > "${fastq2%%.*}_temp_"$count".fastq"
	echo "Done!"

	# Generate list of temporary files to pass to pre-truncation
	for i in *"_temp_"*".fastq"; do
		temp_fastq=`echo $temp_fastq","$i` 
	done
	temp_fastq="[""${temp_fastq:1}""]"

	python2.7 $hictoolPath"HiCtool_pre_truncation.py" -i $temp_fastq -e $restrictionEnzyme -p $threads
	
	cat "${fastq1%%.*}_temp_"*".trunc.fastq" > "${fastq1%%.*}.trunc.fastq"
	cat "${fastq2%%.*}_temp_"*".trunc.fastq" > "${fastq2%%.*}.trunc.fastq"

	rm *"_temp_"*".fastq"

	cat "${fastq1%%.*}"*"log"* > "${fastq1%%.*}_log.txt"
	cat "${fastq2%%.*}"*"log"* > "${fastq2%%.*}_log.txt"

	tot_reads_1=$(awk -F'\t' '{sum+=$1;} END{print sum;}' "${fastq1%%.*}_log.txt")
	tot_reads_2=$(awk -F'\t' '{sum+=$1;} END{print sum;}' "${fastq2%%.*}_log.txt")

	tot_trunc_reads_1=$(awk -F'\t' '{sum+=$3;} END{print sum;}' "${fastq1%%.*}_log.txt")
	tot_trunc_reads_2=$(awk -F'\t' '{sum+=$3;} END{print sum;}' "${fastq2%%.*}_log.txt")

	perc1=$(awk -v n1=$tot_trunc_reads_1 -v ntot1=$tot_reads_1 'BEGIN { print 100*n1/ntot1 }' | cut -c1-5)
	perc2=$(awk -v n1=$tot_trunc_reads_2 -v ntot1=$tot_reads_2 'BEGIN { print 100*n1/ntot1 }' | cut -c1-5)

	read_length=$(awk 'FNR == 1 {print $2}' "${fastq1%%.*}_log.txt")

	printf $fastq1"\n"$tot_reads_1" reads (length = "$read_length" bp); of these:\n  "$tot_trunc_reads_1" ("$perc1"%%)\n\n" > pre_truncation_log.txt
	printf $fastq2"\n"$tot_reads_2" reads (length = "$read_length" bp); of these:\n  "$tot_trunc_reads_2" ("$perc2"%%)" >> pre_truncation_log.txt

	rm *"_temp_"*
	rm "${fastq1%%.*}_log.txt"
	rm "${fastq2%%.*}_log.txt"
	
fi

fastq1_trunc="${fastq1%%.*}.trunc.fastq"
fastq2_trunc="${fastq2%%.*}.trunc.fastq"

echo -n "Aligning "$fastq1_trunc" ... "
(bowtie2 -p $threads -x $genomeIndex $fastq1_trunc -S HiCfile1.sam) 2>HiCfile1_log.txt
echo "Done!"
echo -n "Aligning "$fastq2_trunc" ... "
(bowtie2 -p $threads -x $genomeIndex $fastq2_trunc -S HiCfile2.sam) 2>HiCfile2_log.txt
echo "Done!"

# extracting the headers and read filtering
echo "Filtering and deduplicating HiCfile1.sam ... "
samtools view -H HiCfile1.sam > header1.txt
samtools view -u -h -F 4 -q $quality HiCfile1.sam | \
samtools sort -@ $threads -n - -o - | \
samtools fixmate -m - - | \
samtools sort -@ $threads - -o - | \
samtools markdup -r - HiCfile1_hq_nodup.bam
echo "Done!"
echo "Filtering and deduplicating HiCfile2.sam ... "
samtools view -H HiCfile2.sam > header2.txt
samtools view -u -h -F 4 -q $quality HiCfile2.sam | \
samtools sort -@ $threads -n - -o - | \
samtools fixmate -m - - | \
samtools sort -@ $threads - -o - | \
samtools markdup -r - HiCfile2_hq_nodup.bam
echo "Done!"

echo -n "Building log files ... "
n1=`samtools view HiCfile1_hq_nodup.bam | wc -l | awk '{print $1}'`
n2=`samtools view HiCfile2_hq_nodup.bam | wc -l | awk '{print $1}'`

nt1=`wc -l HiCfile1.sam | awk '{print $1}'`
h1=`wc -l header1.txt | awk '{print $1}'`
ntot1=`expr $nt1 - $h1`

nt2=`wc -l HiCfile2.sam | awk '{print $1}'`
h2=`wc -l header2.txt | awk '{print $1}'`
ntot2=`expr $nt2 - $h2`

perc1=$(awk -v n1=$n1 -v ntot1=$ntot1 'BEGIN { print 100*n1/ntot1 }' | cut -c1-5)
perc2=$(awk -v n2=$n2 -v ntot2=$ntot2 'BEGIN { print 100*n2/ntot2 }' | cut -c1-5)

printf "\n----------\n"$ntot1" reads; of these:\n  "$n1" ("$perc1"%%) aligned with MAPQ>="$quality" and are deduplicated" >> HiCfile1_log.txt
printf "\n----------\n"$ntot2" reads; of these:\n  "$n2" ("$perc2"%%) aligned with MAPQ>="$quality" and are deduplicated" >> HiCfile2_log.txt
echo "Done!"

rm HiCfile1.sam
rm HiCfile2.sam

samtools view HiCfile1_hq_nodup.bam > HiCfile1_hq_nodup.sam
samtools view HiCfile2_hq_nodup.bam > HiCfile2_hq_nodup.sam
rm HiCfile1_hq_nodup.bam
rm HiCfile2_hq_nodup.bam

echo "Selecting paired reads ... "
awk '{print $1}' HiCfile1_hq_nodup.sam | sort > readnames1.txt
awk '{print $1}' HiCfile2_hq_nodup.sam | sort > readnames2.txt
comm -12 readnames1.txt readnames2.txt > paired_reads.txt
echo "Done!"

if [ -z $chunk_size ]
then
	# Select reads of the first sam file that are paires with the second sam file
	echo -n "Extracting paired reads from the first sam file and convert it to bam format ..."
	grep -Fwf paired_reads.txt HiCfile1_hq_nodup.sam | \
	cat header1.txt - | \
	samtools view -b -@ $threads - > HiCfile_pair1.bam
	rm HiCfile1_hq_nodup.sam
	echo "Done!"

	# Select reads of the second sam file that are paired with the first sam file
	echo -n "Extracting paired reads from the second sam file and convert it to bam format ..."
	grep -Fwf paired_reads.txt HiCfile2_hq_nodup.sam | \
	cat header2.txt - | \
	samtools view -b -@ $threads - > HiCfile_pair2.bam
	rm HiCfile2_hq_nodup.sam
	echo "Done!"

elif ! [ -z $chunk_size ] && [ $chunk_size -ge $fastq_lines ]
then
	# Select reads of the first sam file that are paires with the second sam file
	echo -n "Extracting paired reads from the first sam file and convert it to bam format ..."
	grep -Fwf paired_reads.txt HiCfile1_hq_nodup.sam | \
	cat header1.txt - | \
	samtools view -b -@ $threads - > HiCfile_pair1.bam
	rm HiCfile1_hq_nodup.sam
	echo "Done!"

	# Select reads of the second sam file that are paired with the first sam file
	echo -n "Extracting paired reads from the second sam file and convert it to bam format ..."
	grep -Fwf paired_reads.txt HiCfile2_hq_nodup.sam | \
	cat header2.txt - | \
	samtools view -b -@ $threads - > HiCfile_pair2.bam
	rm HiCfile2_hq_nodup.sam
	echo "Done!"

elif ! [ -z $chunk_size ] && [ $chunk_size -lt $fastq_lines ]
then

	# Splitting the paired reads file
	echo -n "Using chunk_size to split the paired read IDs file ..."
	paired_reads=paired_reads.txt
	paired_reads_lines=`wc -l $paired_reads | awk '{print $1}'`
	chunk_size_paired_reads=`expr $chunk_size / 5`
	k=$chunk_size_paired_reads
	count=1
	while [ $k -lt $paired_reads_lines ]
	do
		start=`expr $k - $chunk_size_paired_reads + 1`
		quit=`expr $k + 1`
		sed -n "$start,"$k"p;"$quit"q" $paired_reads > "paired_reads_temp_"$count".txt"
		count=`expr $count + 1`
		k=`expr $k + $chunk_size_paired_reads`
	done
	start=`expr $k - $chunk_size_paired_reads + 1`
	sed -n "$start,"$paired_reads_lines"p" $paired_reads > "paired_reads_temp_"$count".txt"
	echo "Done!"
	
	# Search for paired reads from each temporary file
	echo "Extracting paired reads into temporary sam files and merging ..."
	for i in "paired_reads_temp"*; do
		grep -Fwf $i HiCfile1_hq_nodup.sam > "HiCfile1_${i%%.*}.sam"
		grep -Fwf $i HiCfile2_hq_nodup.sam > "HiCfile2_${i%%.*}.sam"
	done

	cat "HiCfile1_paired_reads_temp"* > HiCfile1_paired.sam
	cat "HiCfile2_paired_reads_temp"* > HiCfile2_paired.sam

	rm *"temp"*
	echo "Done!"
	
	echo -n "Converting the first sam file to bam format ..."
	cat header1.txt HiCfile1_paired.sam | \
	samtools view -b -@ $threads - > HiCfile_pair1.bam
	rm HiCfile1_paired.sam
	echo "Done!"
	
	echo -n "Converting the second sam file to bam format ..."
	cat header2.txt HiCfile2_paired.sam | \
	samtools view -b -@ $threads - > HiCfile_pair2.bam
	rm HiCfile2_paired.sam
	echo "Done!"

	rm HiCfile1_hq_nodup.sam
	rm HiCfile2_hq_nodup.sam

fi

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

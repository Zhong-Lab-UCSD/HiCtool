while getopts h:o:s:p:b:m: option
do
case "${option}"
in
h) hictoolPath=${OPTARG};; # The path of the HiCtool scripts
o) directory=${OPTARG};; # The output directory to work on and save the fend file
s) species=${OPTARG};; # Species under analysis: hg19, hg38, mm9, mm10, dm6, susScr2 are available
p) threads=${OPTARG};; # The number of threads to use
b) gc5Base=${OPTARG};; # gc5Base.bw file with the GC content information (download it from this website: hgdownload.cse.ucsc.edu/gbdb/your_species/bbi/) or GC_info directory path if the GC content information has been already calculated for the reference genome.
m) genomeFasta=${OPTARG};; # genome in fasta format to add the mappability score to the fend file or mappability_info directory path if the mappability score information has been already calculated for the reference genome.
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

echo "Start fend file update: $(date)"
cd $directory

## Splitting restrictionsites.bed to generate separate files to be used either to add GC content information or mappability information or both.
if ! [ -z "$gc5base" ] || ! [ -z "$genomeFasta" ]
then
	echo -n "Splitting restrictionsites.bed, one bed file per each chromosome ... "
	for i in "${chromosomes[@]}"; do
		awk -v var="$i" '(NR>1) && ($1==var)' restrictionsites.bed > $i.bed
	done
	echo "Done!"
fi

### Adding GC content information
if ! [ -z "$gc5Base" ]
then
	if [ $threads -gt "${#chromosomes[@]}" ]
	then
		python_threads="${#chromosomes[@]}"
	else
		python_threads=$threads
	fi
	
	if ! [ -d "$gc5Base" ]
	then
		echo -n "Splitting GC content information, one file per each chromosome ... "
		bigWigToBedGraph $gc5Base gc5Base.bedGraph
		checkMakeDirectory GC_info
		for i in "${chromosomes[@]}"; do
			awk -v var="$i" '(NR>1) && ($1==var)' gc5Base.bedGraph | awk -v OFS='\t' '{print $1, $2, $3, $4}' > GC_info/$i.txt
		done
		echo "Done!"

		python2.7 $hictoolPath"HiCtool_add_fend_gc_content.py" -c $hictoolPath"chromSizes/" -s $species -r $directory -g $directory"GC_info/" -p $python_threads
    else
		python2.7 $hictoolPath"HiCtool_add_fend_gc_content.py" -c $hictoolPath"chromSizes/" -s $species -r $directory -g $gc5Base -p $python_threads
	fi
fi

### Adding mappability score information
if ! [ -z "$genomeFasta" ]
then
	if [ $threads -gt "${#chromosomes[@]}" ]
	then
		python_threads="${#chromosomes[@]}"
	else
		python_threads=$threads
	fi

	if ! [ -d "$genomeFasta" ]
	then
		checkMakeDirectory mappability_info
		python2.7 $hictoolPath"HiCtool_artificial_reads.py" -g $genomeFasta -o $directory"mappability_info/artificial_reads.fastq"
		
		echo -n "Mapping artificial reads ... "
		(bowtie2 -p $threads -x $genomeIndex mappability_info/artificial_reads.fastq -S mappability_info/artificial_reads.sam) 2>mappability_info/artificial_reads.log
		rm mappability_info/artificial_reads.fastq
		echo "Done!"
		
		echo -n "Selecting mapped reads only and generating separate files per each chromosome with mappability score information ... "
		samtools view -F 4 mappability_info/artificial_reads.sam | \
		awk -v OFS='\t' '{print $3, $4-1, $4-1+50, $5}' > mappability_info/artificial_reads_mapped.txt
		
		rm mappability_info/artificial_reads.sam

		for i in "${chromosomes[@]}"; do
			awk -v var="$i" '(NR>1) && ($1==var)' mappability_info/artificial_reads_mapped.txt | awk -v OFS='\t' '{print $1, $2, $3, $4}' > mappability_info/$i.txt
		done
		rm mappability_info/artificial_reads_mapped.txt
		echo "Done!"
		
		python2.7 $hictoolPath"HiCtool_add_fend_mappability.py" -c $hictoolPath"chromSizes/" -s $species -r $directory -a $directory"mappability_info/" -p $python_threads
	else
		python2.7 $hictoolPath"HiCtool_add_fend_mappability.py" -c $hictoolPath"chromSizes/" -s $species -r $directory -a $genomeFasta -p $python_threads
	fi
fi


### Merge files together and generate the final fend bed file
if ! [ -z "$gc5Base" ] && [ -z "$genomeFasta" ]
then
	echo -e 'chr\tstart\tstop\tname\tscore\tstrand\tgc' > header.txt
	echo -n "Merging bed files together and parsing ... "
	for i in "${chromosomes[@]}"; do
		sort -k 2,2n "$i"_gc.bed | cat >> restrictionsites_temp.bed
	done
	
	awk -v OFS='\t' '{print $1, $2, $3, $4, $5, $6, $7 "," $8}' restrictionsites_temp.bed | \
	cat header.txt - > restrictionsites_gc.bed
	
	rm restrictionsites_temp.bed
	rm header.txt

	echo "Done!"
	echo "End fend file update: $(date)"
	echo "Your fend file is restrictionsites_gc.bed"

elif [ -z "$gc5Base" ] && ! [ -z "$genomeFasta" ]
then
	echo -e 'chr\tstart\tstop\tname\tscore\tstrand\tmappability' > header.txt
	echo -n "Merging bed files together, selecting reads with MAP_score >= 0.5 and parsing ... "
	for i in "${chromosomes[@]}"; do
		sort -k 2,2n "$i"_map.bed | cat >> restrictionsites_temp.bed
	done

	awk '(NR>1) && ($7 >= 0.5) && ($8 >= 0.5)' restrictionsites_temp.bed | \
	awk -v OFS='\t' '{print $1, $2, $3, $4, $5, $6, $7 "," $8}' | \
	cat header.txt - > restrictionsites_map.bed

	rm restrictionsites_temp.bed
	rm header.txt

	echo "Done!"
	echo "End fend file update: $(date)"
	echo "Your fend file is restrictionsites_map.bed"

elif ! [ -z "$gc5Base" ] && ! [ -z "$genomeFasta" ]
then
	echo -e 'chr\tstart\tstop\tname\tscore\tstrand\tgc\tmappability' > header.txt
	echo -n "Merging bed files together, selecting reads with MAP_score >= 0.5 and parsing ... "
	for i in "${chromosomes[@]}"; do
		paste -d'\t' "$i"_gc.bed <(awk -F'\t' '{print $7"\t"$8}' "$i"_map.bed) | sort -k 2,2n | cat >> restrictionsites_temp.bed
	done

	awk '(NR>1) && ($9 >= 0.5) && ($10 >= 0.5)' restrictionsites_temp.bed | \
	awk -v OFS='\t' '{print $1, $2, $3, $4, $5, $6, $7 "," $8, $9 "," $10}' | \
	cat header.txt - > restrictionsites_gc_map.bed

	rm restrictionsites_temp.bed
	rm header.txt

	echo "Done!"
	echo "End fend file update: $(date)"
	echo "Your fend file is restrictionsites_gc_map.bed"
fi

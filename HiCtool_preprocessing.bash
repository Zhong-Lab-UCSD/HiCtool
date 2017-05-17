# 2. Converting data from sra format to fastq format

for i in $(ls *.sra); do
    fastq-dump $i --split-3
done

rm *.sra


# 3. Mapping paired-end reads to the reference sequence

num1=0

for i in $(ls *_1.fastq); do
    num1=$(( $num1 + 1))
    string1+=$i,
done

length1=${#string1}-1
mate1s=${string1:0:$length1}


num2=0

for i in $(ls *_2.fastq); do
    num2=$(( $num2 + 1))
    string2+=$i,
done

length2=${#string2}-1
mate2s=${string2:0:$length2}

bowtie2 -p 8 -x index -1 $mate1s -2 $mate2s -S HiCfile.sam

rm *.fastq
samtools view -bS HiCfile.sam > HiCfile.bam
rm *.sam
samtools sort -m 5000000000 HiCfile.bam HiCfile.sort
rm HiCfile.bam


# 4. Removing PCR duplicates from the bam file

samtools sort -m 5000000000 -n HiCfile.sort.bam HiCfile.namesort
rm HiCfile.sort.bam

samtools fixmate HiCfile.namesort.bam HiCfile.fixmate_namesort.bam
rm HiCfile.namesort.bam

samtools sort -m 5000000000 HiCfile.fixmate_namesort.bam HiCfile.fixmate_sort
rm HiCfile.fixmate_namesort.bam

samtools rmdup HiCfile.fixmate_sort.bam HiCfile_noDup.sort.bam
rm HiCfile.fixmate_sort.bam


# 5. Splitting the bam file to separate the two reads in a pair

samtools view -h -f 0x40 HiCfile_noDup.sort.bam > HiCfile_pair1.bam
samtools view -h -f 0x80 HiCfile_noDup.sort.bam > HiCfile_pair2.bam
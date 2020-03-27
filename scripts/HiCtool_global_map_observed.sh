while getopts h:i:b:s:p: option
do
case "${option}"
in
h) hictool=${OPTARG};; # the path to the HiCtool scripts.
i) input_file=${OPTARG};; # Project object file in .hdf5 format.
b) bin_size=${OPTARG};; # Bin size.
s) species=${OPTARG};; # Species.
p) threads=${OPTARG};; # Number of parallel threads to be used.
esac
done

dir=$(dirname "${input_file}")

### Create output directory
checkMakeDirectory(){
echo -e "checking directory: $1"
if [ ! -e "$1" ]; then
	echo -e "\tmakedir $1"
	mkdir -p "$1"
fi
}

start_time=$(date)

output_dir=$dir"/observed_"$bin_size"/" # output directory
checkMakeDirectory $output_dir

python2.7 $hictool"HiCtool_global_map_rows.py" \
-i $input_file \
-o $output_dir \
-b $bin_size \
-s $species \
-c $hictool"chromSizes/" \
-p $threads

while read chromosome size; do
	cat $output_dir"matrix_full_line_observed_"$chromosome".txt" >> $output_dir"HiCtool_observed_global_"$bin_size".txt"
done < $hictool"chromSizes/"$species".chrom.sizes"

echo "Start generating global observed matrix: "$start_time
echo "End generating global observed matrix: $(date)"

### Split each matrix line by single contact matrices
echo -n "Splitting global matrix into single tab separated matrices ..."
while read chromosome_row size_row; do
matrix_line=$output_dir"matrix_full_line_observed_"$chromosome_row".txt"
k=0
while read chromosome_col size_col; do
start=`expr $k + 1`
end=`expr $k + $size_col / $bin_size`
cut -d$'\t' -f"$start"-"$end" $matrix_line > $output_dir"chr"$chromosome_row"_chr"$chromosome_col"_"$bin_size".txt"
k=`expr $k + $size_col / $bin_size`
done < $hictool"chromSizes/"$species".chrom.sizes"
done < $hictool"chromSizes/"$species".chrom.sizes"
echo "Done!"

rm $output_dir"matrix_full_line"*

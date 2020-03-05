while getopts h:i:b:m:q:u:s: option
do
case "${option}"
in
h) hictool=${OPTARG};; # the path to the HiCtool scripts.
i) input_mat_file=${OPTARG};; # the observed global contact matrix in tab delimited format.
b) bin_size=${OPTARG};; # bin size.
m) total_mem=${OPTARG};; # The memory size. Its unit is Megabytes (MB).
q) max_iter=${OPTARG};; # maximum number of iterations performed in the algorithm.
u) row_sum_after_norm=${OPTARG};; # Row sum after the normalization (if not declared automatically calculated as avg_matrix * number of rows).
s) species=${OPTARG};; # species.
esac
done

dir=$(pwd)
total_rows=`wc -l $input_mat_file | awk '{print $1}'`

if [ -z $row_sum_after_norm ]
then
	echo "Rowsum after normalization not declared. Using avg_matrix * number of rows."
	a=$(awk '{for(i=1;i<=NF;i++)x+=$i;print x}' $input_mat_file | tail -n 1)
	row_sum_after_norm=$(( $a / $total_rows ))
else
	echo "Using rowsum after normalization equal to "$row_sum_after_norm"."
fi

checkMakeDirectory(){
echo -e "checking directory: $1"
if [ ! -e "$1" ]; then
	echo -e "\tmakedir $1"
	mkdir -p "$1"
fi
}

start_time=$(date)

output_dir=$dir"/output_ic_mes" # output directory
checkMakeDirectory $output_dir

has_header_line=0  # input file doesn't have header line
has_header_column=0 # input file doesn't have header column

### Calculate biases
cmd=$hictool"Hi-Corrector1.2/bin/ic_mes"
bias_factor_file="$output_dir/output.bias" # output file consists of a vector of bias factors
log_file="$output_dir/output.log" # log file recording the verbose console output of the ic command

echo "$cmd $input_mat_file $total_mem $total_rows $max_iter $has_header_line $has_header_column $bias_factor_file > $log_file"
$cmd $input_mat_file $total_mem $total_rows $max_iter $has_header_line $has_header_column $bias_factor_file > $log_file

### Generate normalized contact matrices
cmd=$hictool"Hi-Corrector1.2/bin/export_norm_data"
normalized_matrix="$output_dir/output_normalized.txt" # output file consists of a vector of bias factors

echo "$cmd $input_mat_file $total_rows $has_header_line $has_header_column $total_mem $bias_factor_file $row_sum_after_norm $normalized_matrix"
$cmd $input_mat_file $total_rows $has_header_line $has_header_column $total_mem $bias_factor_file $row_sum_after_norm $normalized_matrix

echo "Start data normalization time: "$start_time
echo "End data normalization time: $(date)"

output_directory=$dir"/normalized_"$bin_size
checkMakeDirectory $output_directory
output_mat_file="$(echo "$(basename "$input_mat_file")" | sed s/observed/normalized/)"
mv $normalized_matrix $output_directory"/$output_mat_file"


### Split normalized contact matrix in lines by chromosome
echo -n "Splitting global matrix into single tab separated matrices ..."
k=0
while read chromosome size; do
	start=`expr $k + 1`
	end=`expr $k + $size / $bin_size`
	quit=`expr $end + 1`
	sed -n "$start,"$end"p;"$quit"q" $output_directory"/$output_mat_file" > $output_directory"/matrix_full_line_normalized_"$chromosome".txt"
	k=`expr $k + $size / $bin_size`
done < $hictool"chromSizes/"$species".chrom.sizes"

### Split each matrix line by single contact matrices
while read chromosome_row size_row; do
	matrix_line=$output_directory"/matrix_full_line_normalized_"$chromosome_row".txt"
	k=0
	while read chromosome_col size_col; do
		start=`expr $k + 1`
		end=`expr $k + $size_col / $bin_size`
		cut -d$'\t' -f"$start"-"$end" $matrix_line > $output_directory"/chr"$chromosome_row"_chr"$chromosome_col"_"$bin_size".txt"
		k=`expr $k + $size_col / $bin_size`
	done < $hictool"chromSizes/"$species".chrom.sizes"
done < $hictool"chromSizes/"$species".chrom.sizes"
echo "Done!"

rm $output_directory"/matrix_full_line"*

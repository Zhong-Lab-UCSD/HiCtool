while getopts h:i:m:q:r:s: option
do
case "${option}"
in
h) hiCorrector=${OPTARG};; # the path to the Hi-Corrector software.
i) input_mat_file=${OPTARG};; # the observed global contact matrix in tab delimited format.
m) total_mem=${OPTARG};; # The memory size. Its unit is Megabytes (MB).
q) max_iter=${OPTARG};;
r) total_rows=${OPTARG};; # The number of rows or columns of the input chromatin contact frequency matrix to be normalized.
s) row_sum_after_norm=${OPTARG};; # Row sum after the normalization.
esac
done

dir=$(dirname "${input_mat_file}")

### Create output directory
checkMakeDirectory(){
echo -e "checking directory: $1"
if [ ! -e "$1" ]; then
echo -e "\tmakedir $1"
mkdir -p "$1"
fi
}

echo "Start data normalization with Hi-Corrector: $(date)"
start_time=$(date)

output_dir=$dir"/output_ic_mes" # output directory
checkMakeDirectory $output_dir

has_header_line=0  # input file doesn't have header line
has_header_column=0 # input file doesn't have header column


### Calculate biases
cmd=$hiCorrector"bin/ic_mes"
bias_factor_file="$output_dir/output.bias" # output file consists of a vector of bias factors
log_file="$output_dir/output.log" # log file recording the verbose console output of the ic command

echo "$cmd $input_mat_file $total_mem $total_rows $max_iter $has_header_line $has_header_column $bias_factor_file > $log_file"
$cmd $input_mat_file $total_mem $total_rows $max_iter $has_header_line $has_header_column $bias_factor_file > $log_file


### Generate normalized contact matrices
cmd=$hiCorrector"/bin/export_norm_data"
normalized_matrix="$output_dir/output_normalized.txt" # output file consists of a vector of bias factors

echo "$cmd $input_mat_file $total_rows $has_header_line $has_header_column $total_mem $bias_factor_file $row_sum_after_norm $normalized_matrix"
$cmd $input_mat_file $total_rows $has_header_line $has_header_column $total_mem $bias_factor_file $row_sum_after_norm $normalized_matrix

output_mat_file="$(echo "$(basename "$input_mat_file")" | sed s/observed/normalized/)"
mv $normalized_matrix $dir"/$output_mat_file"

echo "Start data normalization time: "$start_time
echo "End data normalization time: $(date)"

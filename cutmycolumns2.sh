# Syntax: bash cutmycolumns.sh <index_file> <big_file> <output_file>

index_file=$1
big_file=$2
output_file=$3

touch temp_result.txt
for i in `cat ${index_file}`; do awk -v i=$i '{print $i}' ${big_file} |paste temp_result.txt - >temp_result2.txt; cp temp_result2.txt temp_result.txt; done
cut -f2- temp_result.txt >$output_file
rm temp_result.txt temp_result2.txt temp

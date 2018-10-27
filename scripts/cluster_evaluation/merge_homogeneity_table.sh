
evalue="1e-50"
inflation='*'
coverage="20"
filter=only_poly

# evalue=$1
# inflation=$2
# coverage=$3
# coverage="20"
cluster_homogeneity_dir=data/clustering_result/Viruses/clustering_parameter_variation_homogeneity_evaluation
output_dir=results/clustering_evaluation/homogeneity_results_merge
mkdir -p $output_dir

merge_file_wt="$output_dir/merge_Viruses_evalue${evalue}_I${inflation}_coverage${coverage}_homogeneity${filter}.csv"

merge_file=`echo "$merge_file_wt" | sed -e 's/\*/_all/g'` # replace the wildcard if the variable by all

echo merge file is $merge_file
if [ ! -f ${merge_file} ] || [ "$force" == true ]; then

  head -1  $cluster_homogeneity_dir/Viruses_evalue_1e-20coverage50_I2_homogeneity.csv > $merge_file

  file_of_interest=$cluster_homogeneity_dir/Viruses_evalue_${evalue}coverage${coverage}_I${inflation}_homogeneity.csv
  for file in $file_of_interest;
  do
    if [ "$filter" == only_poly ]; then
      cat $file | grep -v "cluster_id" | awk '$6 != "0"' >> $merge_file # select only polyprotein
    else
      cat $file | grep -v "cluster_id"  >> $merge_file
    fi
  done

else
  echo merge file already exist
fi

  echo $merge_file
  # head $merge_file
#
# csv_file = args[1]
# evalue= args[2]
# inflation=args[3]
# coverage=args[4]
# output_dir = args[5]
echo RSCRIPT
# Rscript scripts/R_scripts/merge_homogeneity_density_plot.r "$merge_file" "$output_dir/homogeneity_density" "$evalue" "$inflation" "$coverage"




taxonomy_file="data/taxonomy/taxonomy_virus.txt"
alignment_dir="data/alignment/Viruses_1e-5_coverage90_I2/"
window=20

rm ${alignment_dir}stat_of_all_cluster.csv
touch ${alignment_dir}stat_of_all_cluster.csv
for aln in $alignment_dir*aln;
do

  output_file=${aln%.*}.csv
  echo $output_file
  python3 scripts/multiple_alignment_analysis.py $aln $output_file $taxonomy_file $window
  cat $output_file >> ${alignment_dir}stat_of_all_cluster.csv
done

mkdir -p /tmp/$USER/
grep ^'cluster' ${alignment_dir}stat_of_all_cluster.csv | uniq > /tmp/$USER/tmp_file
grep -v ^'cluster' ${alignment_dir}stat_of_all_cluster.csv >> /tmp/$USER/tmp_file
cat /tmp/$USER/tmp_file > ${alignment_dir}stat_of_all_cluster.csv

rm -rf /tmp/$USER/



#Make graph from the csv file result
prefix_plot_file=$(basename $alignment_dir)
Rscript scripts/R_scripts/alignment_analysis_plot.r ${alignment_dir}stat_of_all_cluster.csv $prefix_plot_file

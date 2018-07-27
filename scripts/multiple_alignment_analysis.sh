
set -e # exit if command fail

taxonomy_file="data/taxonomy/taxonomy_virus.txt"
# alignment_dir="data/alignment/Viruses_1e-5_coverage90_I2/"
alignment_dir='data/alignment/Viruses_evalue_1e-40coverage40_I2/'
alignment_dir="data/alignment/RefSeq_download_date_2018-07-21/Viruses_evalue_1e-30coverage20_I2/"
windows="5 10 20 30"
echo window $windows
touch ${alignment_dir}stat_of_all_cluster.csv
rm ${alignment_dir}stat_of_all_cluster.csv

for aln in $alignment_dir*aln;
do

  output_file=${aln%.*}.csv
  echo $output_file

  python3 scripts/multiple_alignment_analysis.py $aln $output_file $taxonomy_file "$windows"
  cat $output_file >> ${alignment_dir}stat_of_all_cluster.csv
done

mkdir -p /tmp/$USER/
grep ^'cluster' ${alignment_dir}stat_of_all_cluster.csv | uniq > tmp
grep -v ^'cluster' ${alignment_dir}stat_of_all_cluster.csv >> tmp
cat tmp > ${alignment_dir}stat_of_all_cluster.csv

rm -rf /tmp/$USER/



#Make graph from the csv file result
prefix_plot_file=$(basename $alignment_dir)
graph_output_dir="results/figures/alignment_analysis_plot/$prefix_plot_file/"

mkdir -p $graph_output_dir

for window in $windows:
do
  echo $window
  Rscript scripts/R_scripts/alignment_analysis_plot.r ${alignment_dir}stat_of_all_cluster.csv $graph_output_dir $window
done

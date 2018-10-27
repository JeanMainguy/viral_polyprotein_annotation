

taxonomy_file="data/taxonomy/taxonomy_virus.txt"
gff_file='data/interpro_results/interproscan-5.30-69.0/domains_viral_sequences.gff3'
# alignment_dir="data/alignment/Viruses_1e-5_coverage90_I2/"

windows="30"


global_alignment_dir='data/alignment/Viruses/RefSeq_download_date_2018-08-13'

for specific_alignment_dir in $global_alignment_dir/*;
do

  # specific_alignment_dir="data/alignment/Viruses/RefSeq_download_date_2018-08-13/Viruses_evalue_1e-40coverage70_I2"
  mkdir -p ${specific_alignment_dir}/stat/
  echo $specific_alignment_dir

  stat_group_file=${specific_alignment_dir}/stat/stat_cleavage_site_groups.csv
  alignement_stat_file=${specific_alignment_dir}/stat/stat_alignments.csv # one line per cluster

  if [ ! -f $alignement_stat_file ]  ; then
    time python3 scripts/multiple_alignment_analysis.py $specific_alignment_dir "$windows" $stat_group_file $alignement_stat_file $taxonomy_file $gff_file
  else
    echo file exist already $positive_negative_file
  fi

  # cross_validation_stat_file=${specific_alignment_dir}/stat/cross_validation_stat_alignments.csv # one line per cluster
  # positive_negative_file=${specific_alignment_dir}/stat/cross_validation_positive_negative.txt # one line per cluster

  #   time python3 scripts/cross_validation_alignment_analysis.py $specific_alignment_dir "$windows" $cross_validation_stat_file $positive_negative_file $taxonomy_file $gff_file
  # else
  #   echo file exist already $positive_negative_file
  # fi

  #Make graph from the csv file result
  # prefix_plot_file=$(basename $specific_alignment_dir)
  # graph_output_dir="results/figures/alignment_analysis_plot/$prefix_plot_file"
  # Rscript scripts/R_scripts/distribution_cluster_category.r $alignement_stat_file $graph_output_dir
  #
  # echo DONE...


done

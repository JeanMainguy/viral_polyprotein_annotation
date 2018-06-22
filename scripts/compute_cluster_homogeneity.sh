



set -e # exit if command fail

force=false # if true even if file exist we recompte them
taxon='Viruses'
# taxon='Alphavirus'
# taxon="Retro-transcribing viruses"
force=false
coverage_min=0

coverage_max=90

coverage_intervalle=10


taxon=${taxon// /_} #replace space by underscore
taxon=${taxon//,/} # replace coma by nothing

###HOMOGENEITY MEASUREMENTs
taxonomy_file="data/taxonomy/taxonomy_virus.txt"
alternative_taxon_id_file="data/taxonomy/heterogeneous_taxon_id_taxonomy_virus.txt"
unclassified_term_file="data/taxonomy/unclassified_terms.txt"


tmpdir=/tmp/$USER/

cluster_dir=data/clustering_result/${taxon}/inflation_test/
homogeneity_dir=data/clustering_result/${taxon}/homogeneity_evaluation/

mkdir -p $homogeneity_dir
mkdir -p $tmpdir$homogeneity_dir
# rm -rf $homogeneity_dir*

for file in ${cluster_dir}*.out;
do
  echo $file
  base=$(basename $file)
  output_file=${homogeneity_dir}${base%.*}_tax_homogeneity.csv
  echo $output_file
  if [ ! -f $output_file ] || [ "$force" == true ]; then

    python3 scripts/taxonomic_homogeneity_checking.py $file ${tmpdir}$output_file $taxonomy_file $alternative_taxon_id_file $unclassified_term_file
    mv ${tmpdir}$output_file ${homogeneity_dir}
  else
    echo homogeneity file $output_file exist already
  fi

done


for c in $(seq $coverage_min $coverage_intervalle $coverage_max);
do
  echo coverage $c
  for file in $homogeneity_dir*coverage${c}*.csv;
  do
    cat $file >> ${homogeneity_dir}merge_coverage${c}.tmp
  done

  grep ^'cluster_id'  ${homogeneity_dir}merge_coverage${c}.tmp | uniq > ${homogeneity_dir}merge_coverage${c}.csv
  grep -v ^'cluster_id' ${homogeneity_dir}merge_coverage${c}.tmp >> ${homogeneity_dir}merge_coverage${c}.csv
  rm ${homogeneity_dir}merge_coverage${c}.tmp

  # #Make graph from the csv file result
  # prefix_plot_file=coverage${c}
  # Rscript scripts/R_scripts/homogeneity_evaluation_plot.r ${homogeneity_dir}coverage80.csv $prefix_plot_file
done


##SCRIPT PLOT HOMOGENEITY
output_dir=results/figures/${taxon}_clustering_homogeneity/
mkdir -p $output_dir

for file_csv in ${homogeneity_dir}merge*csv;
do
  prefix_plot_file=$(basename $file_csv)
  prefix_plot_file=${prefix_plot_file%.*}
  echo plot of  ${output_dir}$prefix_plot_file
  count=`ls -1 ${output_dir}$prefix_plot_file* 2>/dev/null | wc -l`
  if [ $count == 0 ] || [ "$force" == true ];
  then
    Rscript scripts/R_scripts/homogeneity_evaluation_plot.r $file_csv ${output_dir}$prefix_plot_file
  else
    echo plot exist already
  fi
done

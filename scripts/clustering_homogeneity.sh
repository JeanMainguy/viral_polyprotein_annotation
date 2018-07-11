
set -e # exit if command fail

###HOMOGENEITY MEASUREMENTs

taxonomy_file="data/taxonomy/${taxonomy_version}taxonomy_virus.txt"
alternative_taxon_id_file="data/taxonomy/${taxonomy_version}heterogeneous_taxon_id_taxonomy_virus.txt"
unclassified_term_file="data/taxonomy/unclassified_terms.txt"

cluster_dir="data/clustering_result/Viruses/clustering_parameter_variation"

output_dir="data/clustering_result/Viruses/clustering_parameter_variation_homogeneity_evaluation"
mkdir -p $output_dir
summary_file="${output_dir}summary_stat.csv"
stat_protein_file='results/stat_viral_protein/stat_proteins_Viruses.csv'


python3 scripts/taxonomic_homogeneity_checking.py $cluster_dir $output_dir $taxonomy_file $alternative_taxon_id_file $unclassified_term_file $stat_protein_file $summary_file

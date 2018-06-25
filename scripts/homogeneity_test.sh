
set -e # exit if command fail

taxonomy_file="data/taxonomy/taxonomy_virus.txt"
alternative_taxon_id_file="data/taxonomy/heterogeneous_taxon_id_taxonomy_virus.txt"
unclassified_term_file="data/taxonomy/unclassified_terms.txt"

# cluster_dir='data/clustering_result/Retro-transcribing_viruses/inflation_test/'
# homogeneity_dir='data/clustering_result/Retro-transcribing_viruses/homogeneity_evaluation/'
#
# 
# file=data/clustering_result/Viruses/inflation_test/Viruses_1e-5_coverage80_I1_4.out
# file="data/clustering_result/Alphavirus/inflation_test/Alphavirus_1e-5_coverage50_I1_8.out"
# file="data/clustering_result/Alphavirus/inflation_test/Alphavirus_1e-5_coverage90_I1_2.out"
# file='test/inflation_test/Viruses_1e-5_coverage50_I8.out'

file='test/Viruses_1e-5_coverage50_I1_4.out'

#
python3 scripts/visualisation_taxonomic_tree_cluster.py $file

# prefix_plot_file=$(basename $homogeneity_dir)
# Rscript scripts/R_scripts/homogeneity_evaluation_plot.r ${homogeneity_dir}Viruses_1e-5_coverage80_I1_4_tax_homogeneity.csv test_plot

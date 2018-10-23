
set -e # exit if command fail
taxonomy_file="results_db_viral_2018-10-19/genomes_index/taxonomy_virus.txt"
sp_treshold='90'
gff_domain_file='results_db_viral_2018-10-19/intermediate_files/interproscan_results/domains_viral_sequences.gff3'

stat_protein_file="results_db_viral_2018-10-19/viral_protein_stat/stat_proteins_Picornavirales.csv"
# stat_protein_file="results/stat_viral_protein/stat_proteins_Picornavirales.csv"

taxon='Viruses'
taxon='Picornavirales'
stat_output_dir='test/'


python3 scripts/domains_annotation_stat.py $taxon \
                                            ${stat_output_dir}/$RefSeq_download_date \
                                            $taxonomy_file \
                                            $sp_treshold \
                                            $gff_domain_file \
                                            $stat_protein_file

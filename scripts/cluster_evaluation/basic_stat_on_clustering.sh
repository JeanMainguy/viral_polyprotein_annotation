
# cluster_file='data/clustering_result/Viruses/inflation_test/Viruses_evalue_1e-40coverage40_I2.out'
cluster_file=$1
nb_cluster=$(wc -l ${cluster_file})

nb_of_cluster_with_polyprotein=$(grep Peptide ${cluster_file} -c)

nb_of_Polyprotein=$(grep Peptide ${cluster_file} | grep Peptide -o | wc -l)

nb_genome_withannotatedP=$(cat ${cluster_file} | grep Peptide | sed 's/\t\t*/\n/g' | grep 'Peptide' | cut -d'|' -f1 | sort | uniq | wc -l)
nb_of_unannotated_genomecluster_withPoly=$(cat ${cluster_file} | grep Peptide | sed 's/\t\t*/\n/g' | grep 'None' | cut -d'|' -f1 | sort | uniq | wc -l)

nb_protein_cluster_withP=$(grep Peptide ${cluster_file} | grep None -o | wc -l)
singleton=$(cat ${cluster_file} | grep "Peptide" | grep $'\t' -v | wc -l)

echo nb_cluster $nb_cluster
echo  nb_of_cluster_with_polyprotein $nb_of_cluster_with_polyprotein

echo nb of annotated polyprotein $nb_of_Polyprotein
echo nb_protein_clustered_withP $nb_protein_cluster_withP

echo nb_of_unannotated_genomecluster_withPoly $nb_of_unannotated_genomecluster_withPoly
echo singleton $singleton

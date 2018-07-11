#!/bin/bash

set -e # exit if command fail
force=false
TMPDIR=/tmp/$USER/
mkdir -p $TMPDIR

#Cluster Input file
cluster_file='data/clustering_result/Viruses/clustering_parameter_variation/Viruses_evalue_1e-20coverage50_I1_8.out'
faa_db="data/viral_proteins/Viruses_protein_db.faa"

name_dir=$(basename $cluster_file)
name_dir=${name_dir%.*}
alignement_dir=data/alignment/$name_dir/
mkdir -p $alignement_dir

#Stat file to identify protein that have been considered as annotated polyprotein
stat_protein_file="results/stat_viral_protein/stat_proteins_Viruses.csv"

#Python Output file
clusters_with_polyprotein=${alignement_dir}clusters_with_identified_polyprotein.out

if [ ! -f ${clusters_with_polyprotein} ] || [ "$force" == true ]; then # we run the python script only if its output file does not exist...
  python3 scripts/cluster_of_interest_identification.py $cluster_file $clusters_with_polyprotein $stat_protein_file
else
  echo the clusters_with_polyprotein file ${clusters_with_polyprotein} already exist
fi


## -- ALIGNEMENT --

nb_aln_file=$(ls -a ${alignement_dir}*.aln | wc -l)
# if there is already aln file in alignement_dir we don't recompute alignment step
if [ ${nb_aln_file} == '0' ] || [ "$force" == true ]; then
  echo LETS ALIGNE
   sbatch --export=cluster_file=$clusters_with_polyprotein,alignement_dir=$alignement_dir,faa_db=$faa_db scripts/multiple_alignment_proteins.sh
 else
   echo the alignement_dir already contain aln file $clusters_with_polyprotein. no need to recompute this step
 fi


## -- INTERPRO DOMAIN SEARCH --
interproscan_version=$(ls /mirror/interpro/ | grep interproscan*.* -o)
interpro_dir="data/interpro_results/${interproscan_version}/"
mkdir -p $interpro_dir
already_computed_id=${interpro_dir}seq_header_already_processed.txt
final_interpro_result=${interpro_dir}domains_viral_sequences.gff3

if [ ! -f ${already_computed_id} ]; then
  echo creation of a empty already_computed_id file in $interpro_dir because it does not already exist
  echo it means that all id seq found in ${clusters_with_polyprotein} will be process
  touch $already_computed_id
fi

if [ ! -f ${final_interpro_result} ]; then
  echo creation of a empty final_interpro_result file in $interpro_dir because it does not already exist
  # echo it means that all id seq found in ${clusters_with_polyprotein} will be process
  touch $final_interpro_result
fi

cat  ${clusters_with_polyprotein} | sed -e 'y/\t/\n/' > ${interpro_dir}tmp_new_id_list.txt  #replace tab by newline

id_list_to_process=${interpro_dir}complement_new_id_list.txt
comm -23 <(sort ${interpro_dir}tmp_new_id_list.txt) <(sort $already_computed_id) > ${id_list_to_process} # id found  tmp_new_id_list.txt  and not in seq_header_list.txt

nb_new_seq_to_process=$(cat ${id_list_to_process} | wc -l)
echo $nb_new_seq_to_process
if [ ${nb_new_seq_to_process} != '0' ]; then # if new seq is not empty we search interpro domains annotations
    echo Interproscan will process ${nb_new_seq_to_process} new sequences
    echo  sbatch --export=interpro_dir=${interpro_dir},faa_db=$faa_db scripts/interpro_domain_search.sh
   sbatch --export=interpro_dir=${interpro_dir},faa_db=$faa_db scripts/interpro_domain_search.sh
 else
   echo There is no new sequence to process...
 fi
# interpro sbatch file should add complement id list to the main list at the end
# to not recompute the same sequences over and over
rm ${interpro_dir}tmp_new_id_list.txt

echo end
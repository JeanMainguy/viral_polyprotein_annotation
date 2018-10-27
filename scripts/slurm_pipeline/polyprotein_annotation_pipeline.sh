#!/bin/bash
module load pyfasta
module load ncbiblastplus


################################################################################
## FUNCTIONS
################################################################################
exit_if_fail () {
  if [[ $? -ne 0 ]];
  then
    echo FAIL
    exit 1
  fi
}

wait_until_jobs_nb_decrease () {

  nb_jobs_max=$1
  current_nb_jobs=`squeue -u $USER -h -o "%A" | wc -l`
  while (( current_nb_jobs > nb_jobs_max )) ; do
    sleep 5
    squeue -u $USER
    current_nb_jobs=`squeue -u $USER -h -o "%A" | wc -l`
    echo $current_nb_jobs
  done
}

wait_until_jobs_finish () {
  echo lets wait until the sbatch jobs end

  squeue -u $USER
  for job_id in $@ ; do
    cnt=0
    echo $job_id
    while [[ ! $? -ne 0 ]] ; do
      ((cnt++))
      squeue -u $USER
      echo $cnt
      echo "l($cnt) * 2" | bc -l
      sleep `echo "l($cnt) * 2" | bc -l`
      squeue -u $USER -h -o "%A" | grep $job_id
    done

  done
  echo latence sleep
  sleep 10
}

get_RefSeq_date () {
  file=$1
  real=`realpath ${file}` # /proj/viral_polyprotein_annotation/data/taxonomy/RefSeq_download_date_2018-07-21/taxonomy_virus.txt
  real_dir=`dirname $real` #/proj/viral_polyprotein_annotation/data/taxonomy/RefSeq_download_date_2018-07-21/
  RefSeq_download_date=`basename $real_dir` # RefSeq_download_date_2018-07-21
  echo $RefSeq_download_date
}


create_sub_cluster_files () {
  echo create_sub_cluster_files
  echo $1
  echo $2
  cluster_file=$1
  seq_threshold=$2

  nb_seq_in_sub_file=0
  cluster_cmpt=0
  sub_file=${cluster_file}_sub_file$cluster_cmpt
  rm -f $sub_file

  while read l; do
    ((cluster_cmpt++))
    echo $l >> $sub_file
    echo $sub_file
    nb_seq_in_cluster=$(echo $l | grep -o ' ' | wc -l)
    ((nb_seq_in_sub_file+=$nb_seq_in_cluster))
    echo NB SEQ IN SUB FILE $nb_seq_in_sub_file
    if (( nb_seq_in_sub_file > seq_threshold )); then
      echo CREATION OF A NEW SUB FILE
      sub_file=${cluster_file}_sub_file$cluster_cmpt

      nb_seq_in_sub_file=0
    fi

  done <$cluster_file

}

################################################################################
## INPUT AND OUTPUT FILES
################################################################################


seq_fasta_dir='data/viral_proteins/'
stat_output_dir='results/stat_viral_protein/'
ncbi_db_path="/mirror/ncbi/current/"
interpro_path="/mirror/interpro/"


################################################################################
## VARIABLES
################################################################################
## extraction and basic stat
tresholdSP="90"
taxon='Viruses'
# taxon='ssRNA viruses'
# taxon='Retro-transcribing viruses'
# taxon='Picornavirales'

## blast
blast_evalue='1e-5'

## filtering and mcl clustering
coverages='20 30 40 50 60 70'
evalues_filtering='1e-20 1e-30 1e-40 1e-50 1e-60' # 1e-50 1e-20' #'1e-140 1e-160'
inflations='1.4 1.8 2 3'

force=false
# nb jobs launch simultaneously on the cluster
nb_jobs_max=10


echo $taxon
echo Parameters C $coverages E $evalues_filtering I $inflations
sleep 1

taxon_name_for_path=${taxon// /_} #replace space by underscore
taxon_name_for_path=${taxon_name_for_path//,/} # replace coma by nothing

RefSeq_download_date="RefSeq_download_date_`stat -c %y ${ncbi_db_path}genomes/refseq/viral/ | cut -d' ' -f1`"

################################################################################
echo '## TAXONOMY INDEX FILE CREATION'
################################################################################

#output
taxonomy_index_dir="data/taxonomy/$RefSeq_download_date"
taxonomy_file="$taxonomy_index_dir/taxonomy_virus.txt"

bash scripts/create_taxonomy_file.sh $ncbi_db_path $taxonomy_index_dir
exit_if_fail


################################################################################
echo "\n## EXTRACTION VIRAL PROTEINS AND CREATION OF BASIC STAT_FILE\n"
################################################################################

#Output...
sequence_dir="data/viral_proteins/${RefSeq_download_date}"
faa_db="$sequence_dir/${taxon_name_for_path}_protein_db.faa"

stat_output_dir="results/stat_viral_protein/${RefSeq_download_date}/"
stat_protein_file="${stat_output_dir}/stat_proteins_${taxon_name_for_path}.csv"

if [ ! -f $faa_db ] || [ ! -f $stat_protein_file ] ; then
  echo $taxonomy_file
  echo $seq_output_dir
  echo $sequence_dir
  bash scripts/extraction_viral_protein.sh "$taxon" $tresholdSP $taxonomy_file $sequence_dir $stat_output_dir
  exit_if_fail
else
  echo protein sequence database exists already : $faa_db
fi

################################################################################
echo '\n## BLAST ALL VS ALL\n'
################################################################################

#output
blast_result_dir="data/blast_result/${taxon_name_for_path}/${RefSeq_download_date}"

blast_result="${blast_result_dir}/${taxon_name_for_path}_blast_evalue${blast_evalue}.out"

if [ ! -f ${blast_result} ]; then
  bash scripts/blast_allvsall_preparation_array.sh "$taxon" $blast_evalue $sequence_dir $blast_result_dir
  exit_if_fail
else
  echo blast result exist already $blast_result
fi

ls $blast_result
exit_if_fail

################################################################################
echo '## BLAST FILTERING AND MCL CLUSTERING'
################################################################################

clustering_dir="data/clustering_result/${taxon_name_for_path}/${RefSeq_download_date}"

blast_filter_dir="$blast_result_dir/blast_result_filter"

mkdir -p $blast_filter_dir
rm -f $blast_filter_dir/*

python3 scripts/filter_blast_result.py $blast_result $blast_filter_dir "$coverages" "$evalues_filtering"
exit_if_fail

echo $blast_filter_dir

JOB_IDS=()
cluster_files=()
for file in $blast_filter_dir/*;
  do
  echo $file
  name=${taxon_name_for_path}_$(basename $file)
  name=${name%.*}

  for inflation in $inflations; #2 #1.2 1.4 1.6 1.8 2 3 5 8; #$(seq 2 2 8);
  do
    clustering_result=${clustering_dir}/${name}_I${inflation//./_}.out
    cluster_files+=("${clustering_result}")
    # if the file exist we don't recompute the clustering
    if [ ! -f $clustering_result ] || [ "$force" == true ]; then
      echo $clustering_result
      job_id_message=`sbatch --export=filtered_blast_result_file=$file,inflations="$inflation",clustering_dir=$clustering_dir,taxon=$taxon_name_for_path scripts/homologous_clustering.sh`
      echo $job_id_message
      job_id=`echo $job_id_message | cut -d' ' -f4`
      JOB_IDS+=($job_id)
      wait_until_jobs_nb_decrease $nb_jobs_max
    else
      echo the file $clustering_result exist already. We dont recompute the clustering
    fi
  done

done

wait_until_jobs_finish ${JOB_IDS[@]}


################################################################################
echo "## IDENTIFY CLUSTER WITH POLYPROTEINS"
################################################################################
cluster_to_aln=()
alignement_dir_general=data/alignment/$taxon_name_for_path/$RefSeq_download_date
mkdir -p $alignement_dir_general
for cluster_file in "${cluster_files[@]}";
do
  echo $cluster_file
  tail $cluster_file
  sleep 1
  #Python Output file
  name_dir=$(basename $cluster_file)
  name_dir=${name_dir%.*}

  mkdir -p $alignement_dir_general/$name_dir/
  clusters_with_polyprotein=${alignement_dir_general}/$name_dir/clusters_with_identified_polyprotein.out


  if [ ! -f ${clusters_with_polyprotein} ] || [ "$force" == true ]; then # we run the python script only if its output file does not exist...
    echo creation of a file containing only cluster with polyproteins
    python3 scripts/cluster_of_interest_identification.py $cluster_file $clusters_with_polyprotein $stat_protein_file
    exit_if_fail
  else
    echo the clusters_with_polyprotein file ${clusters_with_polyprotein} already exist
  fi

  # Split cluster that have framshiffted proteins
  #Python Output file
  mkdir -p $alignement_dir_general/${name_dir}_splitted/
  clusters_with_polyprotein_splitted=${alignement_dir_general}/${name_dir}_splitted/clusters_with_identified_polyprotein_splitted.out

  if [ ! -f ${clusters_with_polyprotein_splitted} ] || [ "$force" == true ]; then # we run the python script only if its output file does not exist...
    echo  Split cluster that have framshiffted protein
    python3 scripts/split_cluster_with_framshiffted_protein.py $clusters_with_polyprotein $clusters_with_polyprotein_splitted
    exit_if_fail
  else
    echo the clusters_with_polyprotein file ${clusters_with_polyprotein_splitted} already exist
  fi

  cluster_to_aln+=($clusters_with_polyprotein)
  cluster_to_aln+=($clusters_with_polyprotein_splitted)

done


################################################################################
echo "## ALIGNMENTS"
################################################################################
JOB_IDS=()
seq_threshold_by_job=60



for cluster_file in ${cluster_to_aln[@]};
do
  echo $cluster_file
  name_dir=$(basename $cluster_file)
  name_dir=${name_dir%.*}
  alignement_dir=$(dirname $cluster_file)
  nb_aln_file=$(ls -a ${alignement_dir}/*.aln | wc -l)

  # if there is already aln file in alignement_dir we don't recompute alignment step
  if [ ${nb_aln_file} == '0' ] || [ "$force" == true ]; then
    echo SUB FILE CREATION
    create_sub_cluster_files $cluster_file $seq_threshold_by_job
    sleep 5
    for sub_cluster_file in ${cluster_file}_sub_file*;
      do
        echo alignment of $sub_cluster_file
        extension="${sub_cluster_file##*.}"
        var=`echo ${extension//[^0-9]/}`

        job_id_message=`sbatch --export=cluster_file=$sub_cluster_file,alignement_dir=$alignement_dir,faa_db=$faa_db,var=$var scripts/multiple_alignment_proteins.sh`
        job_id=`echo $job_id_message | cut -d' ' -f4`
        JOB_IDS+=($job_id)

        wait_until_jobs_nb_decrease $nb_jobs_max
      done

   else
     echo the alignement_dir already contain aln file $clusters_with_polyprotein. no need to recompute this step
   fi

done

wait_until_jobs_finish ${JOB_IDS[@]}

echo DONE


################################################################################
echo "## INTERPROSCAN: DOMAIN ANNOTATIONS"
################################################################################


interproscan_version=$(ls $interpro_path | grep interproscan*.* -o)
interpro_dir="data/interpro_results/${interproscan_version}"
seq_id_list=${interpro_dir}/tmp_new_id_list.txt
# ####### WARNING
# interpro_dir="test/interpro_results/${interproscan_version}"
mkdir -p $interpro_dir
## MERGE ALL CLUSTER FILE INTO ONE
rm -f ${interpro_dir}/merge_all_tmp_new_id_list.txt
for cluster_file in ${cluster_to_aln[@]};
do
  cat $cluster_file | sed -e 'y/\t/\n/'  >> ${interpro_dir}/merge_all_tmp_new_id_list.txt
done


cat ${interpro_dir}/merge_all_tmp_new_id_list.txt | sort | uniq > $seq_id_list

rm ${interpro_dir}/merge_all_tmp_new_id_list.txt

bash scripts/interproscan_preparation.sh $interpro_dir $clustering_dir $faa_db $seq_id_list
exit_if_fail


echo INTERPRO DONE
exit

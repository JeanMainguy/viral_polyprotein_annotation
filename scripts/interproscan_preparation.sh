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


interpro_dir=$1
cluster_dir=$2
faa_db=$3
seq_id_list=$4

mkdir -p $interpro_dir

already_computed_id=${interpro_dir}/seq_header_already_processed.txt
final_interpro_result=${interpro_dir}/domains_viral_sequences.gff3

if [ ! -f ${already_computed_id} ]; then
  echo creation of a empty already_computed_id file in $interpro_dir/ because it does not already exist
  echo it means that all id seq found in ${cluster_dir} will be process
  touch $already_computed_id
fi

if [ ! -f ${final_interpro_result} ]; then
  echo creation of a empty final_interpro_result file in $interpro_dir/ because it does not already exist
  # echo it means that all id seq found in ${clusters_with_polyprotein} will be process
  touch $final_interpro_result
fi

echo seq_id_list tail
tail $seq_id_list

id_list_to_process=${interpro_dir}/complement_new_id_list.txt
comm -23 <(sort $seq_id_list) <(sort $already_computed_id) > ${id_list_to_process} # id found  tmp_new_id_list.txt  and not in seq_header_list.txt
echo new id_list_to_process tail
tail ${id_list_to_process}
nb_new_seq_to_process=$(cat ${id_list_to_process} | wc -l)
echo nb new seq to process by interproscan : $nb_new_seq_to_process

if [ ${nb_new_seq_to_process} != '0' ]; then # if new seq is not empty we search interpro domains annotations

  echo Interproscan will process ${nb_new_seq_to_process} new sequences
  echo  sbatch --export=interpro_dir=${interpro_dir},faa_db=$faa_db scripts/interpro_domain_search.sh

  job_id_message=`sbatch --export=interpro_dir=${interpro_dir},faa_db=$faa_db scripts/interpro_domain_search.sh`
  job_id=`echo $job_id_message | cut -d' ' -f4`
  wait_until_jobs_finish "$job_id"

else
  echo There is no new sequence to process in interproscan
fi
# interpro sbatch file should add complement id list to the main list at the end
# to not recompute the same sequences over and over
echo rm seq_id_list
rm $seq_id_list

echo end

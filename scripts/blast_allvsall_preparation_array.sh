#!/bin/bash

module load pyfasta
module load ncbiblastplus

# USAGE: bash scripts/blast_allvsall_preparation_array.sh <taxon name> <blast_evalue> <sequence_dir> <blast_result_dir>
# This script takes 4 arguments and performs the following steps:
    ## * split protein database
    ## * launch blast sbatch job in array on the splitted sequence files
    ## * wait until all jobs have ended
    ##   and gather all the blast result into one result file


################################################################################
## FUNCTIONS
################################################################################

wait_until_jobs_finish () {
  echo lets wait until the sbatch jobs end
  squeue -u $USER
  for job_id in $@ ; do
    echo $job_id
    while [[ ! $? -ne 0 ]] ; do
      squeue -u $USER
      sleep 10
      squeue -u $USER -h -o "%A" | grep $job_id
    done

  done
  echo ALL JOB ENDED $@
  squeue -u $USER
  sleep 5
}

################################################################################
## PARAMETERS AND VARIABLES
################################################################################

if [ -z ${1+x} ];
then
  echo "taxon is not provided." # $taxon is used by default";
  exit 1
else
  taxon=$1;
  echo "taxon provided is $taxon";
fi

if [ -z ${2+x} ];
then
  evalue='1e-5'
  echo "Blast evalue threshold is not provided" #. $tresholdSP is used by default";
  exit 1
else
  evalue=$2;
  echo "evalue provided is $evalue";
fi


################################################################################
## INPUT AND OUTPUT FILES
################################################################################

if [ -z ${3+x} ];
then
  echo "sequence_dir is not provided" #. $tresholdSP is used by default";
  exit 1
else
  sequence_dir=$3;
  echo "sequence_dir is $sequence_dir";
fi

if [ -z ${4+x} ];
then
  echo "blast_result_dir is not provided" #. $tresholdSP is used by default";
  exit 1
else
  blast_result_dir=$4;
  echo "blast_result_dir is $blast_result_dir";
fi


taxon=${taxon// /_} #replace space by underscore
taxon=${taxon//,/} # replace coma by nothing

################################################################################
### SPLIT PROTEIN DATABASE
################################################################################
protein_db_faa="${sequence_dir}/${taxon}_protein_db.faa"

splitted_fasta_dir="${sequence_dir}/${taxon}_splitted_fasta_files"

# splitted_fasta_dir="${seq_dir}"
mkdir -p $splitted_fasta_dir
# rm $splitted_fasta_dir*.faa

nb_total_seq=`grep -c ^'>' $protein_db_faa`
nb_seq_per_file=500

nb_of_pieces=$((1 + nb_total_seq / nb_seq_per_file))

echo nb seq in $protein_db_faa is $nb_total_seq
echo nb seq desired per splitted filed is $nb_seq_per_file
echo the db  is then splited into $nb_of_pieces pieces

pyfasta split -n $nb_of_pieces $protein_db_faa

# rename the sub files to make it match with slurm array expectation
for file in ${sequence_dir}/${taxon}_protein_db.*.faa;
 do
   file_name=$(basename $file)
   file_path=$(dirname $file)

   new_file_name=$(echo $file_name | sed -r "s/([^.]*).0*([^0]*)/\1.\2/") # remove leading 0

   echo file_name $file_name
   echo new file_name $new_file_name

   mv $file $splitted_fasta_dir/$new_file_name
done

echo "$splitted_fasta_dir/${taxon}_protein_db.${nb_of_pieces}.faa"

if [ -f "${splitted_fasta_dir}/${taxon}_protein_db..faa" ];
 then
  mv "${splitted_fasta_dir}/${taxon}_protein_db..faa" "$splitted_fasta_dir/${taxon}_protein_db.${nb_of_pieces}.faa"
else
  mv "${splitted_fasta_dir}/${taxon}_protein_db.split.faa" "$splitted_fasta_dir/${taxon}_protein_db.${nb_of_pieces}.faa"
fi

################################################################################
##LAUNCH BLAST JOBS
################################################################################

echo construct blast db
makeblastdb -in ${protein_db_faa} -dbtype prot

echo launch array slurm
submit_message="$(sbatch -a 1-${nb_of_pieces} --export=input_sequence_dir=$sequence_dir,output_blast_dir=$blast_result_dir,taxon=$taxon,evalue=$evalue scripts/blast_allvsall_slurm_array.sh)"

echo $submit_message

jid1=`echo $submit_message | cut -d' ' -f4`
echo $jid1
echo lauched slurm concat that run after the blast array job has ended

wait_until_jobs_finish $jid1

################################################################################
#CONCAT BLAST RESULT FILES
################################################################################

ls $blast_result_dir/

final_result=$blast_result_dir/${taxon}_blast_evalue${evalue}.out
echo creation of final result file $final_result
rm -f $final_result

for file in $blast_result_dir/slurm*${taxon}_blast_evalue${evalue}_*.out;
do
  echo concat $file into final result
  cat $file >> $final_result
done

echo remove intermediary file in dir $blast_result_dir
rm $blast_result_dir/slurm*${taxon}_blast_evalue${evalue}_*.out

echo remove intermediary file in dir data/viral_proteins/${taxon}_splitted_fasta_files
rm -rf data/viral_proteins/${taxon}_splitted_fasta_files


################################################################################
#CREATE SYMBOLIC LINK
################################################################################

echo remove symb link of Blast result if exist
if [ -L data/blast_result/${taxon}/${taxon}_blast_evalue${evalue}.out ];
then
  echo link exist exist.. we remove it
  find  data/blast_result/${taxon}/${taxon}_blast_evalue${evalue}.out -type l -delete
else
  echo link does not exist
fi
echo creation of symblink
real_path_outputdir=`realpath ${blast_result_dir}`
ln -s ${real_path_outputdir}/${taxon}_blast_evalue${evalue}.out data/blast_result/${taxon}/

echo END of blast all vs all

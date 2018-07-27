#!/bin/bash

module load pyfasta
module load ncbiblastplus
set -e # exit if command fail


seq_output_dir='data/viral_proteins/'

taxon='ssRNA viruses'
taxon='Alphavirus'
taxon='Viruses'
evalue='1e-5'
# taxon="Retro-transcribing viruses"
# taxon='ssRNA viruses'

taxon=${taxon// /_} #replace space by underscore
taxon=${taxon//,/} # replace coma by nothing

protein_db_faa="${seq_output_dir}${taxon}_protein_db.faa"

# Extract RefSeq_download_date
real_p=`realpath ${protein_db_faa}` # /proj/viral_poly[...]data/viral_proteins//RefSeq_download_date_2018-07-21/protein_file...
real_dir=`dirname $real_p` #/proj/viral_polyprotein_annotation/data/taxonomy/RefSeq_download_date_2018-07-21/
RefSeq_download_date=`basename $real_dir` # RefSeq_download_date_2018-07-21

echo The RefSeq_download_date is $RefSeq_download_date


splitted_fasta_dir="${seq_output_dir}${taxon}_splitted_fasta_files/"
# splitted_fasta_dir="${seq_output_dir}"
mkdir -p $splitted_fasta_dir
# rm $splitted_fasta_dir*.faa

nb_total_seq=`grep -c ^'>' $protein_db_faa`
nb_seq_per_file=2000

nb_of_pieces=$((1 + nb_total_seq / nb_seq_per_file))

echo nb seq in $protein_db_faa is $nb_total_seq
echo nb seq desired per splitted filed is $nb_seq_per_file
echo the db  is then splited into $nb_of_pieces pieces

pyfasta split -n $nb_of_pieces $protein_db_faa

# rename the sub files to make it match with slurm array expectation
for file in ${seq_output_dir}${taxon}_protein_db.*.faa;
 do
  base=$(basename $file)
  echo base $base
  new_base=$(echo $base | sed "s/${taxon}_protein_db.0*/${taxon}_protein_db./") # remove leading 0
  echo new base $new_base
  mv $file $splitted_fasta_dir$new_base

done
echo "$splitted_fasta_dir${taxon}_protein_db.${nb_of_pieces}.faa"
mv "${splitted_fasta_dir}${taxon}_protein_db..faa" "$splitted_fasta_dir${taxon}_protein_db.${nb_of_pieces}.faa"

echo construct blast db
makeblastdb -in ${protein_db_faa} -dbtype prot


echo launch array slurm
submit_message="$(sbatch -a 1-${nb_of_pieces} --export=taxon=$taxon,fasta_dir=$fasta_dir,evalue=$evalue,RefSeq_download_date=$RefSeq_download_date scripts/blast_allvsall_slurm_array.sh)"

echo $submit_message
echo $submit_message | cut -d' ' -f4
jid1=`echo $submit_message | cut -d' ' -f4`
echo $jid1
echo lauched slurm concat that run after the blast array job has ended
sbatch  --dependency=afterok:$jid1 --export=taxon=$taxon,evalue=$evalue,RefSeq_download_date=$RefSeq_download_date scripts/blast_allvsall_concat_array.sh

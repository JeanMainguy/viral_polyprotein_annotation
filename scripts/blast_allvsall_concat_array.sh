#!/bin/bash
#
#SBATCH --job-name=concat_blast_result
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL
#SBATCH --mem=300M
#SBATCH --output=log/%x-%A_%a.out
#SBATCH --error=log/%x-%A_%a.out

set -e # exit if command fail

#PARAMETERs FOR PROTEIN EXTRACTION
taxonomy_file="data/taxonomy/taxonomy_virus.txt"

if [ -z "$evalue" ];
then
  echo evalue not found default value used
  evalue="1e-5"
fi
echo evalue : $evalue

if [ -z "$taxon" ];
then
  echo taxon not found default value used
  taxon='Viruses'
fi
echo taxon $taxon

if [ -z "$RefSeq_download_date" ];
then
  echo RefSeq_download_date variable is not there we have to quit..
  sleep 20
  exit
fi
echo RefSeq_download_date $RefSeq_download_date

taxon=${taxon// /_} #replace space by underscore
taxon=${taxon//,/} # replace coma by nothing


blast_result_dir="data/blast_result/${taxon}/$RefSeq_download_date"
ls $blast_result_dir

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

echo ----end----

#!/bin/bash
#
#SBATCH --job-name=concat_blast_result
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL
#SBATCH --mem=300M
#SBATCH --output=log/%x-%j.out
#SBATCH --error=log/%x-%j.out

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

taxon=${taxon// /_} #replace space by underscore
taxon=${taxon//,/} # replace coma by nothing


blast_result_dir="data/blast_result/${taxon}_${evalue}"
ls $blast_result_dir

final_result=$blast_result_dir/${taxon}_blast_evalue${evalue}.out
echo creation of final result file $final_result
touch $final_result

for file in $blast_result_dir/slurm*${taxon}_blast_evalue${evalue}_*.out;
do
  echo concat $file into final result
  cat $file >> $final_result
done

echo remove intermediary file in dir $blast_result_dir
rm $blast_result_dir/slurm*${taxon}_blast_evalue${evalue}_*.out

echo remove intermediary file in dir data/viral_proteins/${taxon}_splitted_fasta_files
rm -rf data/viral_proteins/${taxon}_splitted_fasta_files


echo ----end----

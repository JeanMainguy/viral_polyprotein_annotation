#!/bin/bash
#
#SBATCH --job-name=blast_all_vs_all
#SBATCH --cpus-per-task=8
#SBATCH --mail-type=ALL
#SBATCH --mem=300M
#SBATCH --constraint=array-8core
#SBATCH --output=log/%x-%A_%a.out
#SBATCH --error=log/%x-%A_%a.out
#SBATCH --nice=5000

module load ncbiblastplus
set -e # exit if command fail

#PARAMETERs FOR PROTEIN EXTRACTION
taxonomy_file="data/taxonomy/taxonomy_virus.txt"

if [ -z "$fasta_dir" ];
then
  echo seq_output_dir not found default value used
  fasta_dir='data/viral_proteins/'
fi
echo fasta_dir : $fasta_dir
echo fasta splitted directory
fasta_splitted_dir=${fasta_dir}${taxon}_splitted_fasta_files/

if [ -z "$taxon" ];
then
  echo taxon not found default value used
  taxon='Viruses'
fi
echo taxon $taxon

if [ -z "$evalue" ];
then
  echo evalue not found default value used
  evalue="1e-5"
fi
echo evalue : $evalue

if [ -z "$RefSeq_download_date" ];
then
  echo RefSeq_download_date variable is not define we have to quit..
  sleep 20
  exit
fi
echo RefSeq_download_date $RefSeq_download_date

taxon=${taxon// /_} #replace space by underscore
taxon=${taxon//,/} # replace coma by nothing

blast_result_dir="data/blast_result/${taxon}/$RefSeq_download_date"
mkdir -p $blast_result_dir

result_name="${taxon}_blast_evalue${evalue}_${SLURM_ARRAY_TASK_ID}.out"

query="${fasta_splitted_dir}${taxon}_protein_db.${SLURM_ARRAY_TASK_ID}.faa"

protein_db_faa="${fasta_dir}${taxon}_protein_db.faa"

echo Blast $query against all

# /usr/bin/time blastp -db ${protein_db_faa} -query ${protein_db_faa} -evalue "${evalue}" -out $result_all_vs_all -num_threads 16 -outfmt "7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore length positive qcovs qcovhsp qcovus"
/usr/bin/time blastp  -db ${protein_db_faa} -query ${query} -evalue "${evalue}" -out $TMPDIR$result_name -num_threads 16 -outfmt "6 qseqid sseqid qcovs qcovhsp evalue bitscore"


echo mv file from tmp to final dir
mv ${TMPDIR}${result_name}  ${blast_result_dir}/

echo ----end----

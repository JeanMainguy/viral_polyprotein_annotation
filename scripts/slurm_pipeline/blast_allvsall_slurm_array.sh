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

################################################################################
### VARIABLES CHECKING
################################################################################

# TAXON
if [ -z "$taxon" ];
then
  echo taxon not found
  exit 1
  echo default value used
  taxon='Viruses'
fi
echo taxon $taxon

if [ -z "$evalue" ];
then
  echo evalue not found
  exit 1
  echo default value used
  evalue="1e-5"
fi
echo evalue : $evalue

# INPUT
if [ -z "$input_sequence_dir" ];
then
  echo seq_output_dir not found
  exit 1
  echo default value used
  input_sequence_dir='data/viral_proteins'
fi
echo input_sequence_dir : $input_sequence_dir
echo fasta splitted directory
fasta_splitted_dir=${input_sequence_dir}/${taxon}_splitted_fasta_files/

#OUTPUT
if [ -z "$output_blast_dir" ];
then
  echo output_blast_dir not found default value used
  sleep 5
  exit 1
fi

taxon=${taxon// /_} #replace space by underscore
taxon=${taxon//,/} # replace coma by nothing

mkdir -p $output_blast_dir

result_name="${taxon}_blast_evalue${evalue}_${SLURM_ARRAY_TASK_ID}.out"

query="${fasta_splitted_dir}${taxon}_protein_db.${SLURM_ARRAY_TASK_ID}.faa"

protein_db_faa="${input_sequence_dir}/${taxon}_protein_db.faa"


################################################################################
### BLAST: QUERY AGAINST ALL SEQUENCES
################################################################################

echo Blast $query against all

/usr/bin/time blastp  -db ${protein_db_faa} -query ${query} -evalue "${evalue}" -out $TMPDIR$result_name -num_threads 16 -outfmt "6 qseqid sseqid qcovs qcovhsp evalue bitscore"

## echo mv file from tmp to final dir
mv ${TMPDIR}${result_name}  ${output_blast_dir}/

echo ----end----

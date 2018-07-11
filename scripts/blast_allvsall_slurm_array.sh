#!/bin/bash
#
#SBATCH --job-name=viral_proteins_blast_all_vs_all
#SBATCH --cpus-per-task=16
#SBATCH --mail-type=ALL
#SBATCH --mem=300M
#SBATCH --constraint=array-20core
#SBATCH --output=log/%x-%j.out
#SBATCH --error=log/%x-%j.out

module load ncbiblastplus
set -e # exit if command fail

#PARAMETERs FOR PROTEIN EXTRACTION
taxonomy_file="data/taxonomy/taxonomy_virus.txt"
seq_output_dir='data/viral_proteins/'

taxon='ssRNA viruses'
taxon='Alphavirus'
taxon='Viruses'

taxon=${taxon// /_} #replace space by underscore
taxon=${taxon//,/} # replace coma by nothing

# blast parameters
evalue="1e-5"

blast_result_dir="data/blast_result/${taxon}_${evalue}"
mkdir -p $blast_result_dir

result_name="${taxon}_blast_evalue${evalue}_${SLURM_ARRAY_TASK_ID}.out"

query="${seq_output_dir}${taxon}_protein_db.${SLURM_ARRAY_TASK_ID}.faa"

protein_db_faa="${seq_output_dir}${taxon}_protein_db.faa"

echo Blast $query against all

# /usr/bin/time blastp -db ${protein_db_faa} -query ${protein_db_faa} -evalue "${evalue}" -out $result_all_vs_all -num_threads 16 -outfmt "7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore length positive qcovs qcovhsp qcovus"
/usr/bin/time blastp  -db ${protein_db_faa} -query ${query} -evalue "${evalue}" -out $TMPDIR$result_name -num_threads 16 -outfmt "6 qseqid sseqid qcovs qcovhsp evalue bitscore"


echo mv file from tmp to final dir
mv ${TMPDIR}${result_name}  ${blast_result_dir}

echo ----end----

#!/bin/bash

module load pyfasta
module load ncbiblastplus
# set -e # exit if command fail

#PARAMETERs FOR PROTEIN EXTRACTION
taxonomy_file="data/taxonomy/taxonomy_virus.txt"
seq_output_dir='data/viral_proteins/'
stat_output_dir='results/stat_viral_protein/'
taxon='ssRNA viruses'
taxon='Alphavirus'
taxon='Viruses'

# taxon="Retro-transcribing viruses"
# taxon='ssRNA viruses'


taxon=${taxon// /_} #replace space by underscore
taxon=${taxon//,/} # replace coma by nothing

protein_db_faa="${seq_output_dir}${taxon}_protein_db.faa"

splitted_fasta_dir="${seq_output_dir}splitted_fasta_files/"
splitted_fasta_dir="${seq_output_dir}"
mkdir -p $splitted_fasta_dir

nb_of_pieces=300

echo Split db in $nb_of_pieces pieces
pyfasta split -n $nb_of_pieces $protein_db_faa


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
sbatch -a 1-${nb_of_pieces} scripts/blast_allvsall_slurm_array.sh

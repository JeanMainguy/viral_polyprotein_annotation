#!/bin/bash

set -e # exit if command fail

#PARAMETERs FOR PROTEIN EXTRACTION
taxonomy_file="data/taxonomy/taxonomy_virus.txt"
seq_output_dir='data/viral_proteins/'
stat_output_dir='results/stat_viral_protein/'
taxon='ssRNA viruses'
taxon='Alphavirus'
taxon='Viruses'
# taxon='1198013'
# taxon="Retro-transcribing viruses"
# taxon='dsDNA viruses, no RNA stage'
# taxon='Caudovirales'
tresholdSP="90" # length threshold in aa for first peptide annotation corresponding to SP
rm -rf /tmp/$USER/
tmpdir=/tmp/$USER/
mkdir -p ${tmpdir}$seq_output_dir
mkdir -p ${tmpdir}$stat_output_dir

echo extraction of viral proteins...

python3 scripts/viral_protein_extraction.py "${taxon}" ${tmpdir}$seq_output_dir $taxonomy_file  ${tresholdSP} ${tmpdir}$stat_output_dir

echo mv file from tmp to final dir
mv ${tmpdir}${seq_output_dir}*  ${seq_output_dir}
mv ${tmpdir}${stat_output_dir}*  $stat_output_dir

rm -rf /tmp/$USER/

echo END

#To extract line of the stat protein with only protein that have peptide annotation
## cat stat_proteins_Viruses.csv | cut -f2-18 | awk '$4 != "0"' |  awk '$8 != "True"' > stat_potential_polyproteins_Viruses.csv

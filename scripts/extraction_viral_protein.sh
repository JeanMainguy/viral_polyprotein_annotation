#!/bin/bash

set -e # exit if command fail

#PARAMETERs FOR PROTEIN EXTRACTION
taxonomy_file="data/taxonomy/taxonomy_virus.txt"
seq_output_dir='data/viral_proteins/'
stat_output_dir='results/stat_viral_protein/'
taxon='ssRNA viruses'
taxon='Alphavirus'
taxon='Viruses'
taxon="Retro-transcribing viruses"
tresholdSP="90" # length threshold in aa for first peptide annotation corresponding to SP
rm -rf /tmp/$USER/
tmpdir=/tmp/$USER/
mkdir -p ${tmpdir}$seq_output_dir
mkdir -p ${tmpdir}$stat_output_dir

echo extraction of viral proteins...

/usr/bin/time python3 scripts/viral_protein_extraction.py "${taxon}" ${tmpdir}$seq_output_dir $taxonomy_file  ${tresholdSP} ${tmpdir}$stat_output_dir

echo mv file from tmp to final dir
mv ${tmpdir}${seq_output_dir}*  ${seq_output_dir}
mv ${tmpdir}${stat_output_dir}*  $stat_output_dir

rm -rf /tmp/$USER/

echo END

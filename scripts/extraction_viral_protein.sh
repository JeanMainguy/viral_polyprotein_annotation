#!/bin/bash

set -e # exit if command fail

#PARAMETERs FOR PROTEIN EXTRACTION
taxonomy_file="data/taxonomy/taxonomy_virus.txt"
# gff_file='data/interpro_results/interproscan-5.30-69.0/domains_viral_sequences.gff3'
seq_output_dir='data/viral_proteins'
stat_output_dir='results/stat_viral_protein'
# taxon='ssRNA viruses'
taxon='Alphavirus'
taxon='Viruses'
# taxon='1198013'
# taxon="Retro-transcribing viruses"
# taxon='dsDNA viruses, no RNA stage'
# taxon='Caudovirales'
tresholdSP="90" # length threshold in aa for first peptide annotation corresponding to SP

real_taxonomy=`realpath ${taxonomy_file}` # /proj/viral_polyprotein_annotation/data/taxonomy/RefSeq_download_date_2018-07-21/taxonomy_virus.txt
real_tax_dir=`dirname $real_taxonomy` #/proj/viral_polyprotein_annotation/data/taxonomy/RefSeq_download_date_2018-07-21/
RefSeq_download_date=`basename $real_tax_dir` # RefSeq_download_date_2018-07-21

rm -rf /tmp/$USER/
tmpdir=/tmp/$USER/
mkdir -p ${tmpdir}$seq_output_dir
mkdir -p ${tmpdir}$stat_output_dir

mkdir -p ${seq_output_dir}/$RefSeq_download_date
mkdir -p ${stat_output_dir}/$RefSeq_download_date

echo extraction of viral proteins...

python3 scripts/viral_protein_extraction.py "${taxon}" ${tmpdir}$seq_output_dir $taxonomy_file  ${tresholdSP} ${tmpdir}$stat_output_dir #$gff_file

echo mv file from tmp to final dir
mv ${tmpdir}${seq_output_dir}/*  ${seq_output_dir}/$RefSeq_download_date
mv ${tmpdir}${stat_output_dir}/*  ${stat_output_dir}/$RefSeq_download_date

echo remove symb link if exist
if [ -L ${seq_output_dir}/$taxon* ];
then
  echo link exist
  find  ${seq_output_dir}/$taxon* -type l -delete
else
  echo link does not exist
fi
echo creation of symblink
real_path_outputdir=`realpath ${seq_output_dir}/$RefSeq_download_date`
ln -s ${real_path_outputdir}/${taxon}*  ${seq_output_dir}/

echo remove symb link of stats files if exist
if [ -L ${stat_output_dir}/stat_protein*$taxon* ];
then
  echo link stat exist
  find  ${stat_output_dir}/stat_*$taxon* -type l -delete
else
  echo link stat does not exist
fi
echo creation of symblink
real_path_outputdir=`realpath ${stat_output_dir}/$RefSeq_download_date`
ln -s ${real_path_outputdir}/stat_*${taxon}* ${stat_output_dir}/

rm -rf /tmp/$USER/

echo END

#To extract line of the stat protein with only protein that have peptide annotation
## cat stat_proteins_Viruses.csv | cut -f2-18 | awk '$4 != "0"' |  awk '$8 != "True"' > stat_potential_polyproteins_Viruses.csv

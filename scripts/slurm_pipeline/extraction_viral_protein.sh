#!/bin/bash

set -e # exit if command fail

################################################################################
## PARAMETERS AND VARIABLES
################################################################################

#PARAMETERs FOR PROTEIN EXTRACTION
taxonomy_file="data/taxonomy/taxonomy_virus.txt"
seq_output_dir='data/viral_proteins'
stat_output_dir='results/stat_viral_protein'

if [ -z ${1+x} ];
then
  taxon='Viruses'
  echo "taxon is not provided." # $taxon is used by default";
  exit 1
else
  taxon=$1;
  echo "taxon provided is $taxon";
fi

if [ -z ${2+x} ];
then
  tresholdSP="90" # length threshold in aa for first peptide annotation corresponding to SP
  echo "tresholdSP treshold Signal Peptide is not provided" #. $tresholdSP is used by default";
  exit 1
else
  tresholdSP=$2;
    echo "tresholdSP provided is $tresholdSP";
fi

if [ -z ${3+x} ] || [ -z ${4+x} ] || [ -z ${5+x} ];
then
  echo "taxonomy_file, seq_output_dir or stat_output_dir are not provided";
  exit 1
else
  taxonomy_file=$3
  seq_output_dir=$4 #'data/viral_proteins'
  stat_output_dir=$5 #'results/stat_viral_protein'

  echo taxonomy_file $taxonomy_file
  echo seq_output_dir $seq_output_dir
  echo stat_output_dir  $stat_output_dir
fi


tmpdir=/tmp/$USER/protein_extraction_step

rm -rf $tmpdir

mkdir -p ${tmpdir}$seq_output_dir
mkdir -p ${tmpdir}$stat_output_dir

mkdir -p ${seq_output_dir}
mkdir -p ${stat_output_dir}


################################################################################
## PYTHON SCRIPT FOR EXTRACTION
################################################################################
echo extraction of viral proteins...

python3 scripts/viral_protein_extraction.py "${taxon}" ${tmpdir}$seq_output_dir $taxonomy_file  ${tresholdSP} ${tmpdir}$stat_output_dir #$gff_file

echo mv file from tmp to final dir
mv ${tmpdir}${seq_output_dir}/*  ${seq_output_dir}/$RefSeq_download_date
mv ${tmpdir}${stat_output_dir}/*  ${stat_output_dir}/$RefSeq_download_date

# taxon=${taxon// /_} #replace space by underscore
# taxon=${taxon//,/} # replace coma by nothing

# echo remove symb link if exist echo ${seq_output_dir}/$taxon*
# if [ -L ${seq_output_dir}/${taxon}_protein_db.faa ];
# then
#   echo link exist
#   find  ${seq_output_dir}/$taxon* -type l -delete
# else
#   echo link does not exist
# fi
# echo creation of symblink
# real_path_outputdir=`realpath ${seq_output_dir}/$RefSeq_download_date`
# ln -s ${real_path_outputdir}/${taxon}*  ${seq_output_dir}/
#
# echo remove symb link of stats files if exist
# if [ -L ${stat_output_dir}/stat_protein*$taxon* ];
# then
#   echo link stat exist
#   find  ${stat_output_dir}/stat_*$taxon* -type l -delete
# else
#   echo link stat does not exist
# fi
# echo creation of symblink
# real_path_outputdir=`realpath ${stat_output_dir}/$RefSeq_download_date`
# ln -s ${real_path_outputdir}/stat_*${taxon}* ${stat_output_dir}/
#
# rm -rf /tmp/$USER/
#
# echo end of viral protein extraction

#To extract line of the stat protein with only protein that have peptide annotation
## cat stat_proteins_Viruses.csv | cut -f2-18 | awk '$4 != "0"' |  awk '$8 != "True"' > stat_potential_polyproteins_Viruses.csv

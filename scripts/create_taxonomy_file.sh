#!/bin/bash

echo creation of the taxonomy index file


if [ -z ${1+x} ];
then
  # ncbi_current="/mirror/ncbi/current/"
  echo "ncbi_current is not provided." # $ncbi_current is used by default";
  exit 1
else
  ncbi_current=$1;
  echo "ncbi database path is $ncbi_current";
fi

if [ -z ${2+x} ];
then
  # ncbi_current="/mirror/ncbi/current/"
  echo "taxonomy output_dir is not provided." # $ncbi_current is used by default";
  exit 1
else
  output_dir=$2;
  echo "output_dir is $output_dir";
fi


taxonomy_index_file="taxonomy_virus.txt"
taxon_id_gencode_file="taxon_id_gencode.txt"
alternative_taxon_id="heterogeneous_taxon_id_taxonomy_virus.txt"

ncbi_refseq_db=${ncbi_current}/"genomes/refseq/viral/"

#Check if we need to recompute the taxonomy file
outputTime=`stat -c %Y ${output_dir}/${taxonomy_index_file}`
ncbi_refseq_dbTime=`stat -c %Y ${ncbi_refseq_db}`
echo $outputTime
echo $ncbi_refseq_dbTime

set -e # exit if command fail.. this command is here because before the stat commands can fail...

TMPDIR=/tmp/$USER/taxonomy_index
mkdir -p $TMPDIR

if [ ! -f ${output_dir}/$taxonomy_index_file ] || [ $ncbi_refseq_dbTime -gt $outputTime ] || [ "$force" == true ];
then

  echo Creation of the genetic code file

  tar -Oxf ${ncbi_current}taxonomy/new_taxdump/new_taxdump.tar.gz nodes.dmp | cut -d$'\t' -f1,13 > ${TMPDIR}/$taxon_id_gencode_file # extract taxon id and the corresponding genetic code

  echo "creation of taxonomy file with for each taxon id the path to the genbank file in refseq db"
  python3 scripts/taxonomy.py ${TMPDIR}/${taxonomy_index_file} \
                              ${ncbi_refseq_db} \
                              ${TMPDIR}/${taxon_id_gencode_file} \
                              ${TMPDIR}/$alternative_taxon_id 2> ${output_dir}taxonomy_file_creation.log

  echo mv taxonomy files in final dir $output_dir
  mkdir -p $output_dir
  mv ${TMPDIR}/${taxonomy_index_file} $output_dir/
  mv ${TMPDIR}/${taxon_id_gencode_file} $output_dir/
  mv ${TMPDIR}/$alternative_taxon_id $output_dir/

else
  echo The RefSeq database is older than $taxonomy_index_file . So no need to compute again the file
fi

echo remove all symbolic link of the taxonomy dir
find data/taxonomy/  -type l -delete

echo create taxonomy symb link files
real_path_outputdir=`realpath ${output_dir}/`
ln -s $real_path_outputdir/*.txt data/taxonomy/

echo END of create taxonomy script

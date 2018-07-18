#!/bin/bash
#
#SBATCH --job-name=viral_taxonomy
#SBATCH --mail-type=ALL
#SBATCH --output=log/%x-%j.out
#SBATCH --error=log/%x-%j.out


if [ ! -f $1 ];
then
  ncbi_current=$1
else
  ncbi_current="/mirror/ncbi/current/"
fi

force=$2

output_file="taxonomy_virus.txt"
taxon_id_gencode_file="taxon_id_gencode.txt"
alternative_taxon_id="heterogeneous_taxon_id_taxonomy_virus.txt"

ncbi_refseq_db=${ncbi_current}"genomes/refseq/viral/"

# Date of the last RefSeq db download
RefSeq_download_date=`stat -c %y ${ncbi_refseq_db} | cut -d' ' -f1`
output_dir="data/taxonomy/RefSeq_download_date_${RefSeq_download_date}/"
mkdir -p $output_dir

#Check if we need to recompute the taxonomy file
outputTime=`stat -c %Y ${output_dir}${output_file}`
ncbi_refseq_dbTime=`stat -c %Y ${ncbi_refseq_db}`
echo $outputTime
echo $ncbi_refseq_dbTime



TMPDIR=/tmp/$USER/


if [ ! -f ${output_dir}$output_file ] || [ $ncbi_refseq_dbTime -gt $outputTime ] || [ "$force" == true ];
then

  echo Creation of the genetic code file
  tar -Oxf ${ncbi_current}taxonomy/new_taxdump/new_taxdump.tar.gz nodes.dmp | cut -d$'\t' -f1,13 > ${TMPDIR}$taxon_id_gencode_file # extract taxon id and the corresponding genetic code

  echo "creation of taxonomy file with for each taxon id the path to the genbank file in refseq db"
  python3 scripts/taxonomy.py ${TMPDIR}${output_file} ${ncbi_refseq_db} ${TMPDIR}${taxon_id_gencode_file} ${TMPDIR}$alternative_taxon_id #2> ${output_dir}taxonomy_file_creation.log
  echo mv taxonomy files in final dir $output_dir
  mv ${TMPDIR}${output_file} $output_dir
  mv ${TMPDIR}${taxon_id_gencode_file} $output_dir
  mv ${TMPDIR}$alternative_taxon_id $output_dir

else
  echo The RefSeq database is older than $output_file . So no need to compute again the file

fi

echo remove all symbolic link of the taxonomy dir
find data/taxonomy/  -type l -delete

echo create taxonomy symb link files
real_path_outputdir=`realpath ${output_dir}`
ln -s $real_path_outputdir/*.txt data/taxonomy/

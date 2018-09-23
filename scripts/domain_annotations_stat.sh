
set -e # exit if command fail
taxonomy_file="data/taxonomy/taxonomy_virus.txt"
sp_treshold='90'
gff_file='data/interpro_results/interproscan-5.30-69.0/domains_viral_sequences.gff3'

stat_protein_file="results/stat_viral_protein/stat_proteins_Viruses.csv"
# stat_protein_file="results/stat_viral_protein/stat_proteins_Picornavirales.csv"

taxon='Viruses'
# taxon='Picornavirales'
stat_output_dir='results/stat_viral_protein'

#Get name of the folder RefSeqdate
real_taxonomy=`realpath ${taxonomy_file}` # /proj/viral_polyprotein_annotation/data/taxonomy/RefSeq_download_date_2018-07-21/taxonomy_virus.txt
real_tax_dir=`dirname $real_taxonomy` #/proj/viral_polyprotein_annotation/data/taxonomy/RefSeq_download_date_2018-07-21/
RefSeq_download_date=`basename $real_tax_dir` # RefSeq_download_date_2018-07-21
mkdir -p ${stat_output_dir}/$RefSeq_download_date


python3 scripts/domains_annotation_stat.py $taxon ${stat_output_dir}/$RefSeq_download_date $taxonomy_file $sp_treshold $gff_file $stat_protein_file


echo remove symb link of stats files if exist
if [ -L ${stat_output_dir}/domain_stat_${taxon}.csv ];
then
  echo link stat exist, they are removed..
  find  ${stat_output_dir}/domain_*$taxon* -type l -delete
else
  echo link stat does not exist
fi
echo creation of new symblink

real_path_outputdir=`realpath ${stat_output_dir}/$RefSeq_download_date`
ln -s ${real_path_outputdir}/domain_*${taxon}* ${stat_output_dir}/

rm -rf /tmp/$USER/

echo END

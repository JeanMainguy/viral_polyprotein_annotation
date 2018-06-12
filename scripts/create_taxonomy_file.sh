

ncbi_current="/mirror/ncbi/current/"

output_file="data/taxonomy/taxonomy_virus.txt"
taxon_id_gencode_file="data/taxonomy/taxon_id_gencode.txt"
alternative_taxon_id="data/taxonomy/heterogeneous_taxon_id_taxonomy_virus.txt"

ncbi_refseq_db=${ncbi_current}"genomes/refseq/viral/"

#Check if we need to recompute the taxonomy file
outputTime=`stat -c %Y ${output_file}`
ncbi_refseq_dbTime=`stat -c %Y ${ncbi_refseq_db}`
echo $outputTime
echo $ncbi_refseq_db

if [ ! -f $output_file ] || [ $ncbi_refseq_dbTime -gt $outputTime ];
then

  echo Creation of the genetic code file
  tar -Oxf ${ncbi_current}taxonomy/new_taxdump/new_taxdump.tar.gz nodes.dmp | cut -d$'\t' -f1,13 > $taxon_id_gencode_file # extract taxon id and the corresponding genetic code

  echo "creation of taxonomy file with for each taxon id the path to the genbank file in refseq db"
  python3 scripts/taxonomy.py ${output_file} ${ncbi_refseq_db} ${taxon_id_gencode_file} $alternative_taxon_id

else
  echo The RefSeq database is older than $output_file . So no need to compute again the file

fi

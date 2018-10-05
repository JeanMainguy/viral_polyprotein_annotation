
tax_index_file="data/taxonomy/RefSeq_download_date_2018-08-13/taxonomy_virus.txt"
taxon=Flaviviridae

grep $taxon $tax_index_file | while read -r line; do
   echo LINE  $line
   
   echo '===================='
   genom_file_path=`echo $line | rev | cut -f1 -d' ' | rev`

   genom_folder_path=`echo $line | rev | cut -f1 -d' ' | cut -d '/' -f 2-7 | rev`
   
   
   echo $genom_folder_path
   echo GENOME 
   mkdir -p genome_db_test/$genom_folder_path
   cp -L $genom_file_path genome_db_test/$genom_folder_path
done







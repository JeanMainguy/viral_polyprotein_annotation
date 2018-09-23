
tax_index_file=data/taxonomy/taxonomy_virus.txt
taxon=Alphavirus

grep $taxon $tax_index_file | while read -r line; do
   echo 'LINE  $line'
   echo '===================='
   path=`echo $line | rev | cut -f1 -d' ' | rev | cut -d '/' -f 1-9`
   echo $path
   
   cp -R $path genome_db_test/ 
done







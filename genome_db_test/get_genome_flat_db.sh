tax_index_file="results/genomes_index/taxonomy_virus.txt"

taxon=Alphavirus
flat_db=flat_db_$taxon/
mkdir -p ${flat_db}


grep $taxon $tax_index_file | while read -r line; do
   echo LINE  $line

   echo '===================='
   genom_file_path=`echo $line | rev | cut -f1 -d' ' | rev`

   # genom_folder_path=`echo $line | rev | cut -f1 -d' ' | cut -d '/' -f 2-7 | rev`
   #
   #
   # echo $genom_folder_path
   # echo GENOME
   # mkdir -p genome_db_test/$genom_folder_path
   cp -L $genom_file_path $flat_db/
done

# NOT COMPRESS
for f in $flat_db/*.gz;
do
    zcat $f > ${f/.gz/}
    rm $f
done

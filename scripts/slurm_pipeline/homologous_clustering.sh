#!/bin/bash
#
#SBATCH --job-name=mcl_clustering
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=ALL
#SBATCH --mem=300M
#SBATCH --output=log/%x-%j.out
#SBATCH --error=log/%x-%j.out
#SBATCH --nice=2000

set -e # exit if command fail

module load mcl

if [ -z ${filtered_blast_result_file+x} ] || [ -z ${inflations+x} ] || [ -z ${clustering_dir+x} ]  || [ -z ${taxon+x} ];
then
  echo "3 arguments need to be supply to the homologous clustering:"
  echo "blast_result, taxon, coverages, evalues_filtering, and inflations"
  echo 'EXIT'
  sleep 10
  exit 1
else
  echo "blast_result provided is $filtered_blast_result_file";
  echo "inflations provided is $inflations";
  echo "clustering_dir provided is $clustering_dir";
  echo "taxon provided is $taxon"
fi


mkdir -p ${TMPDIR}$clustering_dir/
mkdir -p $clustering_dir/

name=${taxon}_$(basename $filtered_blast_result_file)
name=${name%.*}

abc_file="${TMPDIR}${clustering_dir}/${name}.abc"

cut -f1,2,5 $filtered_blast_result_file > $abc_file

mci_file="${TMPDIR}${clustering_dir}/${name}.mci"

seq_tab="${TMPDIR}${clustering_dir}/${name}.tab"

echo mcxload processing
mcxload -abc $abc_file --stream-mirror --stream-neg-log10 -stream-tf 'ceil(200)' -o $mci_file -write-tab $seq_tab #> /dev/null 2>&1

echo mcxload processing end
for inflation in $inflations; #2 #1.2 1.4 1.6 1.8 2 3 5 8; #$(seq 2 2 8);
do
  echo mcl clustering $inflation with file $name ....
  final_result_mcl="${TMPDIR}${clustering_dir}/${name}_I${inflation//./_}.out"

  /usr/bin/time mcl $mci_file -I $inflation  -use-tab $seq_tab -te 4 -o $final_result_mcl #give output >

  echo  ... mcl clustering $inflation with file $name DONE
  mv  $final_result_mcl $clustering_dir/ #mv file in final dir

done

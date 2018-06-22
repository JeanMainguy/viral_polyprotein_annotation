#!/bin/bash
#
#SBATCH --job-name=virus_homologous_clustering
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=ALL
#SBATCH --mem=300M
#SBATCH --output=log/%x-%j.out
#SBATCH --error=log/%x-%j.out

set -e # exit if command fail

module load mcl

force=true # if true even if file exist we recompte them
taxon='Viruses'
# taxon='Alphavirus'
# taxon="Retro-transcribing viruses"

taxon=${taxon// /_} #replace space by underscore
taxon=${taxon//,/} # replace coma by nothing

evalue="1e-5"

tmpdir=/tmp/${USER}_homologousclustering/
rm -rf $tmpdir
blast_result="data/blast_result/Retro-transcribing_viruses_blast_evalue1e-5.out"
blast_result="data/blast_result/Alphavirus_blast_evalue1e-5.out"
blast_result="data/blast_result/Viruses_blast_evalue1e-5.out"

filter_output_dir=data/blast_result/${taxon}_coverage_filtering/

rm -rf $filter_output_dir
mkdir -p $filter_output_dir

mkdir -p ${tmpdir}$filter_output_dir

coverage_min=0

coverage_max=90

coverage_intervalle=10

echo python filter
python3 scripts/filter_blast_result.py $blast_result  ${tmpdir}$filter_output_dir $coverage_min $coverage_max $coverage_intervalle

clustering_dir="data/clustering_result/${taxon}/inflation_test/"
mkdir -p ${tmpdir}$clustering_dir
mkdir -p $clustering_dir

for file in ${tmpdir}$filter_output_dir*;
do
  name=${taxon}_${evalue}_$(basename $file)
  name=${name%.*}
  abc_file="${tmpdir}${clustering_dir}${name}.abc"

  cut -f1,2,5 $file > $abc_file

  mci_file="${tmpdir}${clustering_dir}${name}.mci"

  seq_tab="${tmpdir}${clustering_dir}${name}.tab"

  echo mcxload processing
  mcxload -abc $abc_file --stream-mirror --stream-neg-log10 -stream-tf 'ceil(200)' -o $mci_file -write-tab $seq_tab > /dev/null 2>&1
  echo mcxload processing end
  for inflation in 1.2 1.4 1.6 1.8 2 3 4 6 8; #$(seq 2 2 8);
  do
    if [ ! -f ${clustering_dir}${name}_I${inflation//./_}.out ] || [ "$force" == true ]; then # if the file exist we don't recompute the clustering

      echo mcl clustering $inflation with file $name ....
      final_result_mcl="${tmpdir}${clustering_dir}${name}_I${inflation//./_}.out"
      /usr/bin/time mcl $mci_file -I $inflation  -use-tab $seq_tab -te 4 -o $final_result_mcl #give output >
      echo  ... mcl clustering $inflation with file $name DONE
      mv ${tmpdir}${clustering_dir}${name}_I${inflation//./_}.out $clustering_dir #mv file in final dir
    else
      echo the file ${clustering_dir}${name}_I${inflation//./_}.out exist already. We dont recompute the clustering
    fi
  done

done

echo mv the files from tmp to the final dir
mv ${tmpdir}$filter_output_dir* $filter_output_dir


rm -rf ${tmpdir}

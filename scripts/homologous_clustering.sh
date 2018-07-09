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



blast_result="data/blast_result/Retro-transcribing_viruses_blast_evalue1e-5.out"
blast_result="data/blast_result/Alphavirus_blast_evalue1e-5.out"
blast_result="data/blast_result/Viruses_blast_evalue1e-5.out"
blast_result="data/blast_result/Viruses_1e-5/concat_slurm_array_Viruses_blast_evalue1e-5_30.out"

filter_output_dir=data/${taxon}/evalue_coverage_filtering/

rm -rf $filter_output_dir
# mkdir -p $filter_output_dir

# TMPDIR=/tmp/${USER}_homologousclustering/
mkdir -p ${TMPDIR}$filter_output_dir

coverage_min=30 #30

coverage_max=80 #80

coverage_intervalle=10

evalues='1e-120 1e-140 1e-160 1e-5 1e-10 1e-20 1e-30 1e-40 1e-50 1e-60 1e-70 1e-80' #
echo python filter
python3 scripts/filter_blast_result.py $blast_result ${TMPDIR}$filter_output_dir $coverage_min $coverage_max $coverage_intervalle $evalues

clustering_dir="data/clustering_result/${taxon}/clustering_parameter_variation/"
mkdir -p ${TMPDIR}$clustering_dir
mkdir -p $clustering_dir


for file in ${TMPDIR}$filter_output_dir*;
do
  name=${taxon}_$(basename $file)
  name=${name%.*}
  abc_file="${TMPDIR}${clustering_dir}${name}.abc"

  cut -f1,2,5 $file > $abc_file

  mci_file="${TMPDIR}${clustering_dir}${name}.mci"

  seq_tab="${TMPDIR}${clustering_dir}${name}.tab"

  echo mcxload processing
  mcxload -abc $abc_file --stream-mirror --stream-neg-log10 -stream-tf 'ceil(200)' -o $mci_file -write-tab $seq_tab > /dev/null 2>&1
  # mcxload -abc $abc_file --stream-mirror --stream-neg-log10 -o $mci_file -write-tab $seq_tab > /dev/null 2>&1
  echo mcxload processing end
  for inflation in 1.2 1.4 1.6 1.8 2 3 5 8; #$(seq 2 2 8);
  do
    if [ ! -f ${clustering_dir}${name}_I${inflation//./_}.out ] || [ "$force" == true ]; then # if the file exist we don't recompute the clustering

      echo mcl clustering $inflation with file $name ....
      final_result_mcl="${TMPDIR}${clustering_dir}${name}_I${inflation//./_}.out"
      /usr/bin/time mcl $mci_file -I $inflation  -use-tab $seq_tab -te 4 -o $final_result_mcl #give output >
      echo  ... mcl clustering $inflation with file $name DONE
      mv ${TMPDIR}${clustering_dir}${name}_I${inflation//./_}.out $clustering_dir #mv file in final dir
    else
      echo the file ${clustering_dir}${name}_I${inflation//./_}.out exist already. We dont recompute the clustering
    fi
  done

done

# echo mv the files from tmp to the final dir
# mv ${TMPDIR}$filter_output_dir* $filter_output_dir

echo check that the filtering make sense
wc -l ${TMPDIR}$filter_output_dir*

# rm -rf ${tmpdir}

#!/bin/bash
#
#SBATCH --job-name=viral_proteins_blast_all_vs_all
#SBATCH --cpus-per-task=16
#SBATCH --mail-type=ALL
#SBATCH --mem=300M
#SBATCH --output=log/%x-%j.out
#SBATCH --error=log/%x-%j.out

module load ncbiblastplus
set -e # exit if command fail

#PARAMETERs FOR PROTEIN EXTRACTION
taxonomy_file="data/taxonomy/taxonomy_virus.txt"
seq_output_dir='data/viral_proteins/'
stat_output_dir='results/stat_viral_protein/'
taxon='ssRNA viruses'
taxon='Alphavirus'
taxon='Viruses'
# taxon="Retro-transcribing viruses"
# taxon='ssRNA viruses'
tresholdSP="90" # length threshold in aa for first peptide annotation corresponding to SP

rm -rf /tmp/$USER/
tmpdir=/tmp/$USER/
mkdir -p ${tmpdir}$seq_output_dir
mkdir -p ${tmpdir}$stat_output_dir

echo extraction of viral proteins...

/usr/bin/time python3 scripts/viral_protein_extraction.py "${taxon}" ${tmpdir}$seq_output_dir $taxonomy_file  ${tresholdSP} ${tmpdir}$stat_output_dir


taxon=${taxon// /_} #replace space by underscore
taxon=${taxon//,/} # replace coma by nothing

# blast parameters
evalue="1e-5"
blast_result_dir="data/blast_result/"

result_name="${taxon}_blast_evalue${evalue}.out"
mkdir -p ${tmpdir}$blast_result_dir
result_all_vs_all=${tmpdir}${blast_result_dir}$result_name

protein_db_faa="${tmpdir}${seq_output_dir}${taxon}_protein_db.faa"

echo construct blast db

makeblastdb -in ${protein_db_faa} -dbtype prot

echo Blast all against all

# /usr/bin/time blastp -db ${protein_db_faa} -query ${protein_db_faa} -evalue "${evalue}" -out $result_all_vs_all -num_threads 16 -outfmt "7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore length positive qcovs qcovhsp qcovus"
/usr/bin/time blastp  -db ${protein_db_faa} -query ${protein_db_faa} -evalue "${evalue}" -out $result_all_vs_all -num_threads 16 -outfmt "6 qseqid sseqid qcovs qcovhsp evalue bitscore"

#Filtering of the blast output
coverage_threshold=80
filtered_blast="${tmpdir}${blast_result_dir}${taxon}_blast_evalue${evalue}_covergae${coverage_threshold}.out"
python3 scripts/filter_blast_result.py $result_all_vs_all  $filtered_blast $coverage_threshold
# #mcl parameter
# inflation="2"
# # =${taxon// /_}
# clustering_dir="data/mcl_clustering_result/"
# mkdir -p ${tmpdir}$blast_result_dir
#
# abc_file="${clustering_dir}${taxon}_${evalue}.abc"
# mci_file="${clustering_dir}${taxon}_${evalue}.mci"
# # row_result_mcl="${clustering_dir}${taxon}_${evalue}_I${inflation}.row"
# final_result_mcl="${clustering_dir}${taxon}_${evalue}_I${inflation//./_}.out"
# seq_tab="${clustering_dir}${taxon}_${evalue}.tab"
#
#
# echo cut
# cut -f1,2,11 $results_all_vs_all > $abc_file # We could pipe the blast to the cut directly...
#
# echo mcxload
# mcxload -abc $abc_file --stream-mirror --stream-neg-log10 -stream-tf 'ceil(200)' -o $mci_file -write-tab $seq_tab
#
#
# echo mcl clustering
# /usr/bin/time mcl $mci_file -I $inflation  -use-tab $seq_tab -te 16 -o $final_result_mcl #give output
#
#


echo mv file from tmp to final dir
mv ${tmpdir}${seq_output_dir}*  ${seq_output_dir}
mv ${tmpdir}${stat_output_dir}*  $stat_output_dir
mv ${tmpdir}${blast_result_dir}* ${blast_result_dir}
# rm -rf /tmp/$USER/

echo ----end----

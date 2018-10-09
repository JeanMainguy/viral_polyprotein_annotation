

################################################################################
## FUNCTIONS
################################################################################
exit_if_fail () {
  if [[ $? -ne 0 ]];
  then
    echo FAIL
    exit 1
  fi
}

################################################################################
## INPUT GENOMES
################################################################################
ncbi_db_path="genome_db_test/"
# interpro_path="/mirror/interpro"
TMPDIR=/tmp/${USER}_polyprotein_annotation/


mkdir -p log/
mkdir -p ${TMPDIR}

################################################################################
## VARIABLES
################################################################################

## extraction and basic stat
tresholdSP="90"
taxon='Flaviviridae'
taxon='Alphavirus'
# taxon='ssRNA viruses'
# taxon='Retro-transcribing viruses'
# taxon='Picornavirales'

## blast evalue start
blast_evalue='1e-5'

## filtering and mcl clustering
coverages='60'
evalues_filtering='1e-60' # 1e-50 1e-20' #'1e-140 1e-160'
inflations='1.8'
split_cluster="false"

echo $taxon
echo Parameters C $coverages E $evalues_filtering I $inflations
sleep 1

taxon_name_for_path=${taxon// /_} #replace space by underscore
taxon_name_for_path=${taxon_name_for_path//,/} # replace coma by nothing

################################################################################
echo '## TAXONOMY INDEX FILE CREATION'
################################################################################

#output
taxonomy_index_dir="results/genomes_index/"
mkdir -p $taxonomy_index_dir
taxonomy_file="$taxonomy_index_dir/taxonomy_virus.txt"

bash scripts/create_taxonomy_file.sh $ncbi_db_path $taxonomy_index_dir
exit_if_fail

################################################################################
echo "\n## EXTRACTION VIRAL PROTEINS AND CREATION OF BASIC STAT_FILE\n"
################################################################################

#Output...
sequence_dir="intermediate_data/viral_proteins/${RefSeq_download_date}"
faa_db="$sequence_dir/${taxon_name_for_path}_protein_db.faa"

stat_output_dir="intermediate_data/stat_viral_protein/${RefSeq_download_date}/"
stat_protein_file="${stat_output_dir}/stat_proteins_${taxon_name_for_path}.csv"

if [ ! -f $faa_db ] || [ ! -f $stat_protein_file ] ; then
  echo $taxonomy_file
  echo $seq_output_dir
  echo $sequence_dir
  bash scripts/extraction_viral_protein.sh "$taxon" $tresholdSP $taxonomy_file $sequence_dir $stat_output_dir
  exit_if_fail
else
  echo protein sequence database exists already : $faa_db
fi

################################################################################
echo '\n## BLAST ALL VS ALL\n'
################################################################################

#output
blast_result_dir="intermediate_data/blast_result/${RefSeq_download_date}"

blast_result="${blast_result_dir}/${taxon_name_for_path}_blast_evalue${blast_evalue}.out"

mkdir -p $blast_result_dir

echo construct blast db
makeblastdb -in ${faa_db} -dbtype prot

if [ ! -f ${blast_result} ]; then
  time blastp -db ${faa_db} -query ${faa_db} \
                                  -evalue "${blast_evalue}" -out $blast_result \
                                  -num_threads 4 \
                                  -outfmt "6 qseqid sseqid qcovs qcovhsp evalue bitscore"
else
  echo blast result exist already $blast_result
fi

ls $blast_result
exit_if_fail


################################################################################
echo '\n## FILTER BLAST OUTPUT\n'
################################################################################

blast_filter_dir="$blast_result_dir/blast_result_filter"

mkdir -p $blast_filter_dir
rm -f $blast_filter_dir/*


python3 scripts/filter_blast_result.py $blast_result $blast_filter_dir "$coverages" "$evalues_filtering"
exit_if_fail
echo $blast_filter_dir


################################################################################
echo '\n## CLUSTERING \n'
################################################################################

# clustering_dir="intermediate_data/clustering_result/${taxon_name_for_path}/${RefSeq_download_date}"

cluster_files=()

for file in $blast_filter_dir/*;
  do


  name=${taxon_name_for_path}_$(basename $file)
  name=${name%.*}
  echo Clustering : $name...

  # Intermediate files for the clustering
  abc_file="${TMPDIR}/${name}.abc"
  mci_file="${TMPDIR}/${name}.mci"
  seq_tab="${TMPDIR}/${name}.tab"

  # echo mcxload processing
  cut -f1,2,5 $file > $abc_file
  mcxload -abc $abc_file --stream-mirror --stream-neg-log10 \
                         -stream-tf 'ceil(200)' -o $mci_file \
                         -write-tab $seq_tab #> /dev/null 2>&1

  for inflation in $inflations; #2 #1.2 1.4 1.6 1.8 2 3 5 8; #$(seq 2 2 8);
  do
    # Creation of the clustering folder
    clustering_dir="results/${RefSeq_download_date}/${name}_I${inflation//./_}/clustering/"
    mkdir -p $clustering_dir/

    clustering_result=${clustering_dir}/all_clusters.out
    cluster_files+=("${clustering_result}")

    # if the file exist we don't recompute the clustering
    if [ ! -f $clustering_result ] || [ "$force" == true ]; then
      echo $clustering_result
      echo mcl clustering $inflation with file $name ...

      mcl $mci_file -I $inflation  -use-tab $seq_tab -te 4 -o $clustering_result

    else
      echo the file $clustering_result exist already. We dont recompute the clustering
    fi
  done
done

################################################################################
echo "## IDENTIFY CLUSTER WITH POLYPROTEINS"
################################################################################
cluster_to_aln=()

for cluster_file in "${cluster_files[@]}";
do
  echo $cluster_file
  #Python Output file
  clustering_path=$(dirname $cluster_file)

  clusters_with_polyprotein=${clustering_path}/clusters_with_identified_polyprotein.out


  if [ ! -f ${clusters_with_polyprotein} ] || [ "$force" == true ]; then # we run the python script only if its output file does not exist...
    echo creation of a file containing only cluster with polyproteins
    python3 scripts/cluster_of_interest_identification.py $cluster_file $clusters_with_polyprotein $stat_protein_file
    exit_if_fail
  else
    echo the clusters_with_polyprotein file ${clusters_with_polyprotein} already exist
  fi

  # Split cluster that have framshiffted proteins
  if [ $split_cluster == "true" ]; then
      # Python Output file
      results_files_specific_cluster=$(dirname $clustering_path)
      splitted_dir=${results_files_specific_cluster}_splitted/clustering/
      mkdir -p $splitted_dir
      clusters_with_polyprotein_splitted=${splitted_dir}/clusters_with_identified_polyprotein_splitted.out

     if [ ! -f ${clusters_with_polyprotein_splitted} ] || [ "$force" == true ]; then # we run the python script only if its output file does not exist...
        echo  Split cluster that have framshiffted protein
        python3 scripts/split_cluster_with_framshiffted_protein.py $clusters_with_polyprotein $clusters_with_polyprotein_splitted $taxonomy_file
        exit_if_fail
    else
        echo the clusters_with_polyprotein file ${clusters_with_polyprotein_splitted} already exist
    fi

    cluster_to_aln+=($clusters_with_polyprotein_splitted)
  fi

  cluster_to_aln+=($clusters_with_polyprotein)

done


################################################################################
echo "## ALIGNMENTS"
################################################################################

aln_dirs=()

for cluster_file in ${cluster_to_aln[@]};
do
  echo $cluster_file
  var=0
  alignment_dir="$(dirname $(dirname $cluster_file))/alignment/"
  mkdir -p $alignment_dir

  aln_dirs+=($alignment_dir)

  while read l;
  do
    echo alignment of cluster $var
    cluster_faa_base=seq_cluster${var}

    if [ ! -f ${alignment_dir}${cluster_faa_base}.aln ] || [ "$force" == true ]; then

      echo $l | sed -e 'y/ /\n/' > ${TMPDIR}/ids.txt #replace tab by newline

      # command from http://bioinformatics.cvr.ac.uk/blog/short-command-lines-for-manipulation-fastq-and-fasta-sequence-files/
      #extract faa seq in a new file
      perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' ${TMPDIR}/ids.txt $faa_db > ${TMPDIR}/$cluster_faa_base.faa

      output=${TMPDIR}/${cluster_faa_base}.aln
      echo ALIGNEMENT of $output
      clustalo -i ${TMPDIR}/$cluster_faa_base.faa -o $output --outfmt=clu --threads=4

    else
      echo the file ${alignment_dir}${cluster_faa_base}.aln  exist already. We dont recompute the alignment
    fi
    mv $output ${alignment_dir}
    ((var++)) #increment var to know which cluster we are

  done < $cluster_file


done

################################################################################
echo "## ALIGNMENT ANALYSIS AND PROPAGATION OF CLEAVAGE SITES"
################################################################################

gff_file='data/interpro_results/interproscan-5.30-69.0/domains_viral_sequences.gff3'
# alignment_dir="data/alignment/Viruses_1e-5_coverage90_I2/"

windows="30"

for aln_dir in ${aln_dirs[@]};
do
    stat_output_dir="$(dirname $aln_dir)/alignment_analysis"
    mkdir -p ${stat_output_dir}

    reannotated_genome_dir="$(dirname $aln_dir)/reannotated_genome"
    mkdir -p ${reannotated_genome_dir}

    stat_group_file=${stat_output_dir}/stat_cleavage_site_groups.csv
    alignement_stat_file=${stat_output_dir}/stat_alignments.csv # one line per cluster

    # if [ ! -f $alignement_stat_file ] ; then
      python3 scripts/multiple_alignment_analysis.py $aln_dir "$windows" $stat_group_file $alignement_stat_file $taxonomy_file $reannotated_genome_dir
    # else
    #   echo file exist already $alignement_stat_file
    # fi
done

# DONE
rm -r $TMPDIR
echo DONE

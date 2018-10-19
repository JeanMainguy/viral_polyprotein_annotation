#!/bin/bash

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
genetic_code_path="genome_db_test/taxonomy/new_taxdump/"

genbank_files_db_path="genome_db_test/genomes/refseq/viral/"
RefSeq_structure="True"
# genbank_files_db_path='flat_db_Alphavirus/'
# RefSeq_structure="False"
#Structure of the genbank files database : True or False
# False is a list of genbank files in a folder
# True same structure as in RefSeq:
# ── genomes
#     └── refseq
#         └── viral
#             ├── Aedes_flavivirus
#             │   └── latest_assembly_versions
#             │       └── GCF_000885715.1_ViralProj39601
#             │           └── GCF_000885715.1_ViralProj39601_genomic.gbff.gz


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
# taxon='Alphavirus'

# taxon='ssRNA viruses'
# taxon='Retro-transcribing viruses'
# taxon='Picornavirales'

## blast evalue start
blast_evalue='1e-5'

## filtering and mcl clustering
coverages='60 70'
evalues_filtering='1e-60' # 1e-50 1e-20' #'1e-140 1e-160'
inflations='1.8'
split_cluster="false"

# Get results Folder.. results_[name of the gb databse]_[databse_date]
databse_date="`stat -c %y ${genbank_files_db_path} | cut -d' ' -f1`"
results_folder="results_db_$(basename $genbank_files_db_path)_$databse_date"


echo $taxon
echo Parameters C $coverages E $evalues_filtering I $inflations
echo Results are written in $results_folder

taxon_name_for_path=${taxon// /_} #replace space by underscore
taxon_name_for_path=${taxon_name_for_path//,/} # replace coma by nothing

################################################################################
echo '## GENOME INDEX FILE CREATION'
################################################################################

#output
taxonomy_index_dir="${results_folder}/genomes_index/"
mkdir -p $taxonomy_index_dir

alternative_taxon_id="$taxonomy_index_dir/heterogeneous_taxon_id_taxonomy_virus.txt"
taxon_id_gencode_file="$taxonomy_index_dir/taxon_id_gencode.txt"

if [ ! -f ${taxon_id_gencode_file} ] || [ "$force" == true ];
then
    echo Creation of the genetic code file
     # extract taxon id and the corresponding genetic code
    tar -Oxf ${genetic_code_path}/new_taxdump.tar.gz nodes.dmp | cut -d$'\t' -f1,13 > ${TMPDIR}/taxon_id_gencode_file.tmp
    exit_if_fail

    mv ${TMPDIR}/taxon_id_gencode_file.tmp $taxon_id_gencode_file
fi

taxonomy_file="$taxonomy_index_dir/taxonomy_virus.txt"


#Check if we need to recompute the taxonomy file
outputTime=`stat -c %Y ${taxonomy_file}`
ncbi_refseq_dbTime=`stat -c %Y ${genbank_files_db_path}`
# force=true

if [ ! -f ${taxonomy_file} ] || [ $ncbi_refseq_dbTime -gt $outputTime ] || [ "$force" == true ];
then

  echo "creation of genome index file: taxon id, path to the genbank file, genome taxonomy.. "
  python3 scripts/taxonomy.py ${TMPDIR}/taxonomy_index_file.tmp \
                              ${genbank_files_db_path} \
                              ${taxon_id_gencode_file} \
                              ${TMPDIR}/alternative_taxon_id.tmp \
                              $RefSeq_structure 2> ${taxonomy_index_dir}/taxonomy_file_creation.log
  exit_if_fail

  mv ${TMPDIR}/taxonomy_index_file.tmp $taxonomy_file
  mv ${TMPDIR}/alternative_taxon_id.tmp $alternative_taxon_id

else
  echo The RefSeq database is older than $taxonomy_file . So no need to compute again the file
fi


################################################################################
echo "\n## EXTRACTION VIRAL PROTEINS AND CREATION OF BASIC STAT_FILE\n"
################################################################################

#Output...
sequence_dir="${results_folder}/intermediate_files/viral_proteins/${RefSeq_download_date}"
faa_db="$sequence_dir/${taxon_name_for_path}_protein_db.faa"

stat_output_dir="${results_folder}/intermediate_files/stat_viral_protein/${RefSeq_download_date}/"
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
blast_result_dir="${results_folder}/intermediate_files/blast_result/${RefSeq_download_date}"

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

# clustering_dir="${results_folder}/intermediate_files/clustering_result/${taxon_name_for_path}/${RefSeq_download_date}"

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
    clustering_dir="${results_folder}/${RefSeq_download_date}/${name}_I${inflation//./_}/clustering/"
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
      mv $output ${alignment_dir}
    else
      echo the file ${alignment_dir}${cluster_faa_base}.aln  exist already. We dont recompute the alignment
    fi

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

    reannotated_genome_dir="$(dirname $aln_dir)/reannotated_genomes"
    mkdir -p ${reannotated_genome_dir}
    mkdir -p $TMPDIR/${reannotated_genome_dir}

    stat_group_file=${stat_output_dir}/stat_cleavage_site_groups.csv
    alignement_stat_file=${stat_output_dir}/stat_alignments.csv # one line per cluster

    # if [ ! -f $alignement_stat_file ] ; then
      python3 scripts/multiple_alignment_analysis.py $aln_dir "$windows" $stat_group_file $alignement_stat_file $taxonomy_file $TMPDIR/${reannotated_genome_dir}
      exit_if_fail
    # else
    #   echo file exist already $alignement_stat_file
    # fi
    mv $TMPDIR/${reannotated_genome_dir}/* ${reannotated_genome_dir}/

    if [ -x "$(command -v aha)" ]; then
        echo 'Conversion of the alignment visualisation file in html'
        for f in $stat_output_dir/*.aln;
        do
            cat $f | aha >  ${f//.*/.html}
        done
    fi
done

# DONE
rm -r $TMPDIR
echo DONE

module load clustalomega
module load pyfasta



taxon='Alphavirus'
faa_db="data/viral_proteins/Alphavirus_protein_db.faa"
cluster_file="data/clustering_result/Alphavirus/inflation_test/Alphavirus_1e-5_coverage60_I1_6.out"
# taxon="Retro-transcribing viruses"
# cluster_file='data/clustering_result/Retro-transcribing_viruses/inflation_test/Retro-transcribing_viruses_1e-5_coverage90_I3.out'
# faa_db="data/viral_proteins/Retro-transcribing_viruses_protein_db.faa"

# taxon='Viruses'
# cluster_file='data/clustering_result/Viruses/inflation_test/Viruses_1e-5_coverage90_I2.out'
# faa_db="data/viral_proteins/Viruses_protein_db.faa"


if [ -z "$TMPDIR" ];
then
  TMPDIR=/tmp/$USER/
fi

taxon=${taxon// /_} # replace space by underscore
taxon=${taxon//,/} # replace coma by nothing

evalue="1e-5"

seq_faa_dir=${TMPDIR}cluster_faa_sequences/
# seq_faa_dir=data/alignment/cluster_faa_sequences/
mkdir -p $seq_faa_dir

echo cluster protein extraction and multiplealignment

name_dir=$(basename $cluster_file)
name_dir=${name_dir%.*}
alignement_dir=data/alignment/$name_dir/

mkdir -p $alignement_dir
mkdir -p ${TMPDIR}$alignement_dir

interpro_result_dir=data/interpro_results/$name_dir/
mkdir -p $interpro_result_dir
mkdir -p ${TMPDIR}$interpro_result_dir


echo TMPDIR ${TMPDIR}

var=0
while read l;
do

  if [[ $l = *"Peptide"* ]]; then

    echo $l | sed -e 'y/ /\n/' > ${seq_faa_dir}ids.txt #replace tab by newline
    cluster_faa_base=seq_cluster${var}
    # command from http://bioinformatics.cvr.ac.uk/blog/short-command-lines-for-manipulation-fastq-and-fasta-sequence-files/
    #extract faa seq in a new file
    perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' ${seq_faa_dir}ids.txt $faa_db > ${seq_faa_dir}$cluster_faa_base.faa

    output=${TMPDIR}${alignement_dir}${cluster_faa_base}.aln
    interpro_output=${TMPDIR}${interpro_result_dir}${cluster_faa_base}
    echo ALIGNEMENT of $output

    /usr/bin/time clustalo -i ${seq_faa_dir}$cluster_faa_base.faa -o $output --outfmt=clu  --wrap=150
    /usr/bin/time /localmirror/monthly/interpro/interproscan-5.29-68.0/interproscan.sh --appl PFAM,CDD,ProDom,SMART,ProSiteProfiles -i ${seq_faa_dir}$cluster_faa_base.faa -b ${interpro_output} -f GFF3

  fi
  ((var++)) #increment var to know which cluster we are
done <$cluster_file

mv ${TMPDIR}${alignement_dir}* ${alignement_dir}
mv ${TMPDIR}${interpro_result_dir}* ${interpro_result_dir}
rm -rf ${TMPDIR}

# for aln in $alignement_dir*aln;
# do
#   name=${aln%.*}.faa
#   /apps/emboss/6.6.0/bin/seqret -sequence $aln -sformat1 clustal -osformat2 fasta  -outseq $name
# done

#!/bin/bash
#
#SBATCH --job-name=interproscan
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=ALL
#SBATCH --mem=5GB
#SBATCH --output=log/%x-%j.out
#SBATCH --error=log/%x-%j.out

set -e # exit if command fail

# VARIABLE EXPORTED...

if [ -z "$interpro_dir" ];
then
  interpro_dir="data/interpro_results/test"
fi
echo interpro_dir : $interpro_dir

if [ -z "$faa_db" ];
then
  faa_db="data/viral_proteins/Viruses_protein_db.faa"
fi
echo faa_db : $faa_db

if [ -z "$TMPDIR" ];
then
  TMPDIR=/tmp/$USER
  mkdir -p ${TMPDIR}
fi
echo TMPDIR : $TMPDIR

if [ -z "$SLURM_JOBID" ];
then
  SLURM_JOBID=TEST
fi
echo SLURM_JOBID : $SLURM_JOBID

list_seq_id=${interpro_dir}/complement_new_id_list.txt
final_interpro_result=${interpro_dir}/domains_viral_sequences.gff3


echo Interproscan search with $(wc -l $list_seq_id) sequences
# command from http://bioinformatics.cvr.ac.uk/blog/short-command-lines-for-manipulation-fastq-and-fasta-sequence-files/
# extract faa seq in a new file
perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' $list_seq_id $faa_db > ${TMPDIR}/sequences.faa
interpro_output=${TMPDIR}/domains_viral_sequences_${SLURM_JOBID}
/usr/bin/time /localmirror/monthly/interpro/interproscan*/interproscan.sh -cpu 8 --appl PFAM,CDD,ProDom,SMART,ProSiteProfiles -i ${TMPDIR}/sequences.faa -b ${interpro_output} -f GFF3

# GFF3 interpro output store fasta sequence at the end of the file... we don't want this part as the result will be merge with other
# get line number where FASTA sequence start
line_nb=$(cat ${interpro_output}.gff3 | grep '##FASTA' -n)

echo $line_nb
line_nb=${line_nb%:*}
((line_nb--)) # remove 1 to the line to exclude ##FASTA line
echo $line_nb
head -$line_nb ${interpro_output}.gff3 > ${interpro_output}_no_fasta.gff3

## MERGE NEW FILE TO FINAL RESULT FILE
echo MERGE NEW FILE TO FINAL RESULT FILE
cat ${interpro_output}_no_fasta.gff3 >> $final_interpro_result

## MV ORIGINAL OUTPUT IN FINAL INTERPRO DIR
echo MV ORIGINAL OUTPUT IN FINAL INTERPRO DIR
interpro_dir=$(dirname "${final_interpro_result}")
mv $interpro_output.gff3  $interpro_dir/
cat $list_seq_id >> ${interpro_dir}/seq_header_already_processed.txt

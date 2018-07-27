#!/bin/bash
#
#SBATCH --job-name=multiple_alignement
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=ALL
#SBATCH --mem=4GB
#SBATCH --output=log/%x-%j.out
#SBATCH --error=log/%x-%j.out

# set -e # exit if command fail

module load clustalomega

#Receive cluster_file
#Alignment DIR
# and faa_db
echo multiplealignment dir $alignement_dir
echo cluster_file $cluster_file
echo faa db $faa_db
force=true
if [ -z "$alignement_dir" ];
then
  alignement_dir=data/alignment/Viruses_evalue_1e-20coverage50_I1_8/
  mkdir -p $alignement_dir

fi

if [ -z "$cluster_file" ];
then
  cluster_file='data/clustering_result/Viruses/clustering_parameter_variation/Viruses_evalue_1e-20coverage50_I1_8.out'
fi
if [ -z "$faa_db" ];
then
  faa_db="data/viral_proteins/Viruses_protein_db.faa"
fi

if [ -z "$TMPDIR" ];
then
  TMPDIR=/tmp/$USER
  mkdir -p ${TMPDIR}
fi

echo multiplealignment dir $alignement_dir
echo cluster_file $cluster_file
echo faa db $faa_db


mkdir -p $alignement_dir

echo TMPDIR ${TMPDIR}

var=0
while read l;
do
  cluster_faa_base=seq_cluster${var}

  if [ ! -f ${alignement_dir}${cluster_faa_base}.aln ] || [ "$force" == true ]; then

    echo $l | sed -e 'y/ /\n/' > ${TMPDIR}/ids.txt #replace tab by newline
    # command from http://bioinformatics.cvr.ac.uk/blog/short-command-lines-for-manipulation-fastq-and-fasta-sequence-files/
    #extract faa seq in a new file
    perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' ${TMPDIR}/ids.txt $faa_db > ${TMPDIR}/$cluster_faa_base.faa

    output=${TMPDIR}/${cluster_faa_base}.aln
    echo ALIGNEMENT of $output
    /usr/bin/time clustalo -i ${TMPDIR}/$cluster_faa_base.faa -o $output --outfmt=clu --threads=8

  else
    echo the file ${alignement_dir}${cluster_faa_base}.aln  exist already. We dont recompute the clustering
  fi
  echo $var
  ((var++)) #increment var to know which cluster we are
  echo $var

done <$cluster_file

mv ${TMPDIR}/*.aln ${alignement_dir}
rm -rf ${TMPDIR}

aln_dir_name=`basename ${alignement_dir}`
##Symblink of the aln dir
echo remove symb link of the aln dir global files if exist
if [ -L data/alignment/${aln_dir_name} ];
then
  echo link clustering exist
  find  data/alignment/${aln_dir_name} -type l -delete
else
  echo link clustering does not exist
fi
echo creation of symblink
real_path_outputdir=`realpath ${alignement_dir}`
ln -s ${alignement_dir} data/alignment/

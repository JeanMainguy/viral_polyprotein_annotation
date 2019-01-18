#!/usr/bin/env nextflow


// genetic_code_file = Channel.fromPath("genome_db_test/taxonomy/new_taxdump/new_taxdump.tar.gz")
params.genetic_code_file = "genome_db_test/taxonomy/new_taxdump/new_taxdump.tar.gz"
params.refseq_structure = "True"
params.genbank_files_db_path = "$PWD/genome_db_test2/genomes/refseq/viral/"
genetic_code_file = file(params.genetic_code_file)
// genetic_code_file = file(params.genetic_code_file)

process create_genome_index {
    input:
        file genetic_code_file

    output:
        file 'taxonomy_index_file.txt' into genome_index_file
        stdout db_date

        """
        tar -Oxf ${params.genetic_code_file} nodes.dmp | \
            cut -d\$'\t' -f1,13 > taxon_id_gencode_file

        python3 $PWD/scripts/taxonomy.py taxonomy_index_file.txt \
                                    ${params.genbank_files_db_path} \
                                    taxon_id_gencode_file \
                                    heterogeneous_taxon_id_taxonomy_virus.txt \
                                    ${params.refseq_structure} > python.std

        db=${params.genbank_files_db_path}
        database_date="`stat -c %y \$db | cut -d' ' -f1`"
        db_download_date="\$(basename ${params.genbank_files_db_path})_\$database_date"
        echo \$db_download_date

        """
    }

db_date.subscribe { print "I say..  $it" }

db_date = db_date.map {it -> it.replace('\n', '')}

params.taxon = "Alphavirus"
taxon_name_for_path = params.taxon
                            .replaceAll(/,/, "")
                            .replaceAll(/ /, "_")


genome_index_file.into { genome_index_extraction; genome_index_analysis}

params.thresholdSP ="90" // ## Signal Peptide length threshold in aa

process extract_viral_prot {
     // publishDir '.', saveAs: { it == "foo.*" ? "foos/$it" : "bars/$it" }

    input:
        file genome_index_file from genome_index_extraction
    output:
        file "stat_proteins_${taxon_name_for_path}.csv" into stat_protein_file
        file "list_irrelevant_annotation_${taxon_name_for_path}.txt" into irrelevant_cds_file
        file 'list_annotated_polyprotein.txt' into annotated_polyprotein_list
        file "${taxon_name_for_path}_protein_db.faa" into faa_db
    """
    python3 $PWD/scripts/viral_protein_extraction.py "${params.taxon}" \
                                                ${genome_index_file} \
                                                ${taxon_name_for_path}_protein_db.faa \
                                                --sp_treshold ${params.thresholdSP} \
                                                --stat_output_dir "." \
                                                --irrelevant_annotation_list list_irrelevant_annotation_${taxon_name_for_path}.txt \
                                                --polyprotein_list_file list_annotated_polyprotein.txt
    """

}

params.blast_evalue='1e-5'
params.blast_num_threads = "4"
params.chunkSize = 50

// faa_proteins = Channel.value(faa_db)
// faa_proteins.println()
faa_db.into { fasta_to_split; blast_faa_db; faa_db_aln } // duplicate channel

// faa_proteins = faa_db_test.value()
// faa_proteins.println()

fasta_to_split.splitFasta(by: params.chunkSize).set { fasta }

process blast_all_vs_all {
    input:
        file 'query.fa' from fasta
        file blast_faa_db
    output:
        file "blast_result.out" into blast_result

"""
makeblastdb -in ${blast_faa_db} -dbtype prot
time blastp -db ${blast_faa_db} -query ${blast_faa_db} \
                                  -evalue "${params.blast_evalue}" \
                                  -out "blast_result.out" \
                                  -num_threads ${params.blast_num_threads} \
                                  -outfmt "6 qseqid sseqid qcovs qcovhsp evalue bitscore"
"""
}
//                                  -out "blast_result.out" \

blast_result_merge = blast_result.collectFile(name: "${taxon_name_for_path}_blast_evalue${params.blast_evalue}.out", newLine: true)

params.coverages = '60'
params.evalues_filtering = '1e-60'

coverages = params.coverages.split(' ')
evalues_filtering = params.evalues_filtering.split(' ')

process filter_blast_result {
    input:
        file blast_result_merge
        each coverage from coverages
        each evalue_cutoff from evalues_filtering
    output:
        set val("evalue${evalue_cutoff}_coverage${coverage}"), file("evalue_${evalue_cutoff}coverage${coverage}.out") into filtered_results


    """
    python3 $PWD/scripts/filter_blast_result.py $blast_result_merge \
                                           "." \
                                           "$coverage" "$evalue_cutoff"

    """
}

process mcl_file_preparation {
    input:
        set val(parameter), file('blast_result.out') from filtered_results

    output:
        set val(parameter), file('seq.mci'), file("seq.tab") into formated_mcl_files

    """
      cut -f1,2,5 blast_result.out > seq.abc
      mcxload -abc seq.abc --stream-mirror --stream-neg-log10 \
                             -stream-tf 'ceil(200)' -o seq.mci \
                             -write-tab seq.tab
    """

}

params.inflations = '1.8'
inflations = params.inflations.split(' ')
params.mcl_num_threads = "4"

process mcl_clustering {
    input:
        set val(parameter), file('seq.mci'), file("seq.tab")  from formated_mcl_files
        each inflation from inflations
    output:
        set val("${parameter}_I${inflation}"), file('all_clusters.out') into clustering_results


    """
        mcl seq.mci -I $inflation -use-tab seq.tab -te ${params.mcl_num_threads} -o all_clusters.out
    """
}

process identify_polyprotein_clusters {
    input:
        set val(parameter), file('cluster_file') from clustering_results
        file stat_protein_file
    output:
        set val(parameter), file('clusters_with_identified_polyprotein.out') into polyprotein_clustering
    """
    python3 $PWD/scripts/cluster_of_interest_identification.py cluster_file \
                                                        clusters_with_identified_polyprotein.out \
                                                        $stat_protein_file
    """
}


// ADD PROCESS SPLIT CLUSTER WITH FRAMESHOFTED as an option

// split line of clustering result
// one element of the channel corresponds to one cluster
// filter to remove singleton
polyprotein_clustering.splitText().filter { it -> it[1] =~ /.*\t.*/}.set {clusters_of_polyprotein}

params.aln_num_threads = 2

process mltp_alignement {

    input:
        file faa_db from faa_db_aln
        set val(parameter), val(cluster_line) from clusters_of_polyprotein
    output:
        set val(parameter), file("cluster.aln") into alignments

    """
    echo "$cluster_line" | sed 's/\\t\\t*/\\n/g' > cluster_ids.txt # replace tab by newline

    # command from http://bioinformatics.cvr.ac.uk/blog/short-command-lines-for-manipulation-fastq-and-fasta-sequence-files/
    # extract faa seq in a new file
    # escape back slash and dollar sign with back slash
    perl -ne 'if(/^>(\\S+)/){\$c=\$i{\$1}}\$c?print:chomp;\$i{\$_}=1 if @ARGV' cluster_ids.txt \
                                                                        $faa_db > cluster.faa

    clustalo -i cluster.faa -o cluster.aln --outfmt=clu --threads=${params.aln_num_threads}


    """

}

all_aln_by_param = alignments.groupTuple() // groups aln file by the parameter

params.window = 30
params.confidence_score_treshold = 4

process alignment_analysis {
    // "nextflow_result_${db_date}/${taxon_name_for_path}/${parameter}/reannotated_genomes/"
    db_date = db_date.replaceAll(/\n/, "")
    publishDir "nextflow_result_${db_date}/${taxon_name_for_path}/${parameter}/reannotated_genomes/", pattern: '*.gff'
    publishDir "nextflow_result_${db_date}/${taxon_name_for_path}/${parameter}/reannotated_genomes/", pattern: '*.gbff'

    publishDir "nextflow_result_${db_date}/${taxon_name_for_path}/${parameter}/alignment_analysis/", pattern: 'visu*.aln'
    publishDir "nextflow_result_${db_date}/${taxon_name_for_path}/${parameter}/alignment_analysis/", pattern: 'stat*.csv'


    input:
        set val(parameter), file('cluster*.aln') from all_aln_by_param
        file genome_index_file from genome_index_analysis
        val db_date
    output:
        // file '*.csv' into stat_files
        file '*.gff' into gff_files
        file '*.gbff' into gb_files
        file "stat*.csv" into stat_files
        set val(parameter),  file('visu*.aln') into visuaisation_files

    """
    python3 $PWD/scripts/multiple_alignment_analysis.py  "." $genome_index_file \
                                                    -w ${params.window} \
                                                    -t ${params.confidence_score_treshold} \
                                                    --sp_treshold ${params.thresholdSP}
                                                    ## --blacklist black_list_file
                                                    ## --interpro_domains gff_domain_file.csv

    """
}


process convert_aln_to_html {
    publishDir "nextflow_result_${db_date}/${taxon_name_for_path}/${parameter}/alignment_analysis/", pattern: '*.html'

    input:
        set val(parameter),  file(aln) from visuaisation_files
        val db_date
    output:
        file "*.html" into html_visuaisation_files


    """
    python3 $PWD/scripts/html_conversion.py $aln
    """
}

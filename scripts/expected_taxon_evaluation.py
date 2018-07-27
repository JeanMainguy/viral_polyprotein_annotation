#!/usr/bin/env python3
import taxonomy as tax
import logging, sys
import taxonomic_homogeneity_checking as tax_hom
import csv, os, re

def write_unexpected_annotated_genome(unexpected_annotated_genomes, unexpected_annotated_genome_file, taxonomy_file):

    header = [
    "taxon_id",
    "protein_ids",
    "taxonomy",
    "gb_file"
    ]
    with open(unexpected_annotated_genome_file, "w") as out:
        out.write('\t'.join(header)+"\n")
        for dic_info  in tax.getAllRefseqFromTaxonIdList(list(unexpected_annotated_genomes), taxonomy_file):
            taxon_id = dic_info['taxon_id']
            line = [taxon_id, '|'.join(unexpected_annotated_genomes[taxon_id]),dic_info['taxonomy'], dic_info['gb_file']]

            out.write('\t'.join(line)+"\n")


def get_genomes_clustered_with_poly(cluster_file, annotated_genomes, annotated_genomes_filtered):
    clustered_genome_with_poly = set()
    clustered_genome_with_expected_poly = set()
    with open(cluster_file, 'r') as fl:
        #665102|YP_003622544.1|381	1046572|YP_004958227.1|360	536084|YP_001960955.1|356
        for i, cluster_line in enumerate(fl):
            cluster_has_polyprotein = False

            cluster_elements = cluster_line.rstrip().split("\t")

            genomes_id_in_cluster = {element.split('|')[0]:element.split('|')[1] for element in cluster_elements}
            genomes_id_in_cluster_set = set(genomes_id_in_cluster)
            for genome_id, protein_id in genomes_id_in_cluster.items():

                # print(genome_id, protein_id)
                if genome_id in annotated_genomes_filtered and protein_id in annotated_genomes_filtered[genome_id]:
                    clustered_genome_with_poly |= genomes_id_in_cluster_set
                    clustered_genome_with_expected_poly |= genomes_id_in_cluster_set
                    break
                elif genome_id in annotated_genomes and protein_id in annotated_genomes[genome_id]:
                    clustered_genome_with_poly |= genomes_id_in_cluster_set
                    break
        return clustered_genome_with_poly, clustered_genome_with_expected_poly, i+1

def compute_expected_evaluation_on_cluster(cluster_file, annotated_genomes, annotated_genomes_filtered, set_genomes_expected):
    # we retrieve the coverage value and the inflation from the file name
    re_result = re.search("evalue_(1e[-\d]+)coverage(\d+)_I([\d_]+).out", cluster_file)
    #if the regex doesn't work an expception will be raised
    if re_result:
        evalue=  re_result.group(1)
        coverage = re_result.group(2)
        inflation = re_result.group(3)
    else:
        raise NameError('parsing cluster file name failed %s' % cluster_file)

    clustered_genome_with_poly, clustered_genome_with_expected_poly, nb_cluster = get_genomes_clustered_with_poly(cluster_file, annotated_genomes, annotated_genomes_filtered)

    cluster_iter = tax_hom.clusterFileParser(cluster_file, annotated_genomes)

    expected_and_clustered_intersection = set_genomes_expected &  clustered_genome_with_poly
    expected_and_clustered_union = set_genomes_expected | clustered_genome_with_poly


    # The traditional F-measure or balanced F-score (F1 score) is the harmonic mean of precision and recall
    # https://en.wikipedia.org/wiki/F1_score
    # print(precision)
    # print(sensitivity)
    precision =len(expected_and_clustered_intersection)/len(clustered_genome_with_poly)
    sensitivity = len(expected_and_clustered_intersection)/len(set_genomes_expected)
    F1_score = 2*(precision*sensitivity/(precision + sensitivity))

    expected_and_clustered_intersection_filtered = set_genomes_expected &  clustered_genome_with_expected_poly
    expected_and_clustered_union_filtered = set_genomes_expected | clustered_genome_with_expected_poly


    # The traditional F-measure or balanced F-score (F1 score) is the harmonic mean of precision and recall
    # https://en.wikipedia.org/wiki/F1_score
    # print(precision)
    # print(sensitivity)
    precision_filtered =len(expected_and_clustered_intersection_filtered)/len(clustered_genome_with_expected_poly)
    sensitivity_filtered = len(expected_and_clustered_intersection_filtered)/len(set_genomes_expected)
    F1_score_filtered = 2*(precision_filtered*sensitivity_filtered/(precision_filtered + sensitivity_filtered))

    dict_to_write = {"nb_cluster":nb_cluster,
                    "Evalue":evalue,
                    "coverage":coverage,
                    "inflation":inflation,
                    "nb_expected_taxid":len(set_genomes_expected),

                    "genome_clustered_with_annotated":len(clustered_genome_with_poly),
                    "expected_and_clustered_intersection": len(expected_and_clustered_intersection),
                    "expected_and_clustered_union":len(expected_and_clustered_union),
                    "precision": precision,
                    'sensitivity': sensitivity,
                    "F1_score":F1_score,

                    "genome_clustered_with_expecteted_annotated":len(clustered_genome_with_expected_poly),
                    "expected_and_clustered_intersection_filtered": len(expected_and_clustered_intersection_filtered),
                    "expected_and_clustered_union_filtered":len(expected_and_clustered_union_filtered),
                    "precision_filtered": precision_filtered,
                    'sensitivity_filtered': sensitivity_filtered,
                    "F1_score_filtered":F1_score_filtered
                    }
    return dict_to_write


def getGenomeTaxonmy(taxonomy_file):
    leaves_taxonomy = {}
    with open(taxonomy_file, 'r') as taxon_handle:
        for l in taxon_handle:
            genome_id, genome_name, taxonomy, genetic_code, gbff_file = l.split("\t")
            taxonomy = taxonomy.split(';')
            leaves_taxonomy[genome_id] = tuple(taxonomy)
    return leaves_taxonomy

def getExpectedGenomeList(taxonomy_file, expected_taxons):
    '''
    Get genome id list in which we expect to find polyprotein
    '''
    expected_tax_ids = []
    with open(taxonomy_file, 'r') as taxon_handle:
        for l in taxon_handle:
            genome_id, genome_name, taxonomy, genetic_code, gbff_file = l.split("\t")
            taxonomy = taxonomy.split(';')
            for taxon in taxonomy:
                if taxon in expected_taxons:
                    expected_tax_ids.append(genome_id)
                    continue
            # leaves_taxonomy[genome_id] = tuple(taxonomy)
    return expected_tax_ids

def getExpectedTaxon(expected_taxon_file):
    expected_taxons=[]
    with open(expected_taxon_file, 'r') as fl:
        for l in fl:
            taxon, nb_peptides = l.rstrip().split('\t')
            if nb_peptides != '0':
                expected_taxons.append(taxon)
        return expected_taxons




def setComparison(E, A):
    """
    Compare two set of taxon id E and A:
    Output:
        tuple (difference_E, difference_A, intersection )
    overlaping and non-overlapping number of elements
    """
    return (len(E - A), len(A - E), len(A & E))
    # return (E - A, A - E, A & E)

if __name__ == '__main__':
    """
    How well genome represented in cluster with
    at least one identified polyprotein match
    with taxon where we expect to find polyprotein?
    """
    unclassified_term_file = "data/taxonomy/unclassified_terms.txt"
    taxonomy_file = "data/taxonomy/taxonomy_virus.txt"
    expected_taxon_file='etc/viruses_w_polyproteins.txt'
    stat_protein_file='results/stat_viral_protein/stat_proteins_Viruses.csv'
    unexpected_annotated_genome_file = "results/stat_viral_protein/unexpected_annotated_genomes.csv"


    expected_taxons = getExpectedTaxon(expected_taxon_file)
    expected_tax_ids = getExpectedGenomeList(taxonomy_file, expected_taxons)

    cluster_file="data/clustering_result/Viruses/inflation_test/Viruses_evalue_1e-20coverage50_I2.out"
    cluster_dir='data/clustering_result/Viruses/clustering_parameter_variation'
    cluster_files = (os.path.join(cluster_dir, f) for f in os.listdir(cluster_dir) if f.endswith('.out'))

    nb_genome_expected = len(expected_tax_ids)
    set_genomes_expected = set(expected_tax_ids)

    print("NUMBER OF GENOMES WHERE WE EXPECT POLYPROTEIN",len(expected_tax_ids))
    annotated_genomes = tax_hom.getAnnotatedProteinTaxId(stat_protein_file)
    print("NUMBER OF GENOMES ANNOTATED", len(annotated_genomes))
    t = setComparison(set(expected_tax_ids), set(annotated_genomes))
    print(t)

    # annotated_genomes_filtered = {k:v for k,v in annotated_genomes.items() if k not in unexpected_annotated_genome}
    # unexpected_annotated_genomes = {k:v for k,v in annotated_genomes.items() if k in unexpected_annotated_genome}

    expected_annotated_genome =  set(annotated_genomes) - set(annotated_genomes) | set_genomes_expected # union
    annotated_genomes_filtered = {k:v for k,v in annotated_genomes.items() if k in expected_annotated_genome}
    # will be used then to select cluster of interest
    unexpected_annotated_genomes = {k:v for k,v in annotated_genomes.items() if k not in expected_annotated_genome}

    print('annotated_genomes_filtered', len(annotated_genomes_filtered))
    print('unexpected_annotated_genomes', len(unexpected_annotated_genomes))
    write_unexpected_annotated_genome(unexpected_annotated_genomes, unexpected_annotated_genome_file, taxonomy_file)
    expected_annotated = set(expected_tax_ids) &  set(annotated_genomes)
    # annotated_genomes_filtered = expected_annotated
    print('NUMBER OF GENOMES EXPECTED AND ANNOTATED', len(expected_annotated))
    print('NUMBER OF GENOMES UNEXPECTED AND ANNOTATED', len(unexpected_annotated_genomes))
    # expected_and_annotated_genomes = {}
    # for genome in annotated_genomes:
    #     if genome in expected_annotated:
    #         expected_and_annotated_genomes[genome] = annotated_genomes[genome]

    ##initiate output file
    output_file = 'results/clustering_evaluation/expected_taxon_cluster_evaluation.csv'

    fl_out = open(output_file, 'w')
    header = ["nb_cluster",
                        "Evalue",
                        "coverage",
                        "inflation",
                        "nb_expected_taxid",

                        "genome_clustered_with_annotated",
                        "expected_and_clustered_intersection",
                        "expected_and_clustered_union",
                        "precision",
                        'sensitivity',
                        "F1_score",

                        "genome_clustered_with_expecteted_annotated",
                        "expected_and_clustered_intersection_filtered",
                        "expected_and_clustered_union_filtered",
                        "precision_filtered",
                        'sensitivity_filtered',
                        "F1_score_filtered"]

    csv_output = csv.DictWriter(fl_out, fieldnames=header, delimiter='\t')
    csv_output.writeheader()


    for nb_file, cluster_file in enumerate(cluster_files):

        print(nb_file, cluster_file)
        dict_to_write = compute_expected_evaluation_on_cluster(cluster_file, annotated_genomes, annotated_genomes_filtered, set_genomes_expected)
        csv_output.writerow(dict_to_write)
        # print(t)
        # print('nb of cluster with polyprotein annotated and expected', ie)
        # print('number of genome clustered with annotated and expected genome', len(clustered_genome_with_poly_exp))
        # t = setComparison(set(expected_tax_ids), clustered_genome_with_poly_exp)
        # print(t)

    fl_out.close()

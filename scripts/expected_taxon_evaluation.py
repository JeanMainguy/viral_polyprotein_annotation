import logging, sys
import taxonomic_homogeneity_checking as tax_hom
import csv, os, re


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

    expected_annotated = set(expected_tax_ids) &  set(annotated_genomes)
    print('NUMBER OF GENOMES EXPECTED AND ANNOTATED', len(expected_annotated))

    # expected_and_annotated_genomes = {}
    # for genome in annotated_genomes:
    #     if genome in expected_annotated:
    #         expected_and_annotated_genomes[genome] = annotated_genomes[genome]

    ##initiate output file
    output_file = 'results/clustering_evaluation/expected_taxon_cluster_evaluation.csv'

    fl_out = open(output_file, 'w')
    header = [  "nb_cluster",
                    "nb_cluster_with_poly",
                    "nb_clustered_taxid_with_poly",
                    "nb_expected_taxid",
                    "complement_clustered_taxid_with_poly",
                    "complement_expected_taxid",
                    "expected_and_clustered_intersection",
                    "Evalue",
                    "coverage",
                    "inflation",
                    "precision",
                    'sensitivity',
                    "F1_score"]

    csv_output = csv.DictWriter(fl_out, fieldnames=header, delimiter='\t')
    csv_output.writeheader()


    for nb_file, cluster_file in enumerate(cluster_files):
        print(nb_file)
        # we retrieve the coverage value and the inflation from the file name
        re_result = re.search("evalue_(1e[-\d]+)coverage(\d+)_I([\d_]+).out", cluster_file)
        #if the regex doesn't work an expception will be raised
        if re_result:
            evalue=  re_result.group(1)
            coverage = re_result.group(2)
            inflation = re_result.group(3)
        else:
            raise NameError('parsing cluster file name failed %s' % cluster_file)


        cluster_iter = tax_hom.clusterFileParser(cluster_file, annotated_genomes)
        # cluster_iter_exp = tax_hom.clusterFileParser(cluster_file, expected_and_annotated_genomes)
        clustered_genome_with_poly = set()
        i = 0
        nb_cluster = 0
        for cluster in cluster_iter:
            nb_cluster += 1
            if cluster['nb_polyprotein'] > 0:
                i += 1
                clustered_genome_with_poly |= set(cluster['genome_ids']) # add the two set
        # ie = 0
        # clustered_genome_with_poly_exp = set()
        # for cluster in cluster_iter_exp:
        #     if cluster['nb_polyprotein'] > 0:
        #         ie += 1
        #         clustered_genome_with_poly_exp |= set(cluster['genome_ids']) # add the two set
        # print('nb of cluster with polyprotein annotated', i)
        # print('number of genome clustered with annotated', len(clustered_genome_with_poly))
        t = setComparison(set(expected_tax_ids), clustered_genome_with_poly)
        expected_and_clustered_intersection = set_genomes_expected &  clustered_genome_with_poly
        expected_and_clustered_union = set_genomes_expected | clustered_genome_with_poly

        precision = len(expected_and_clustered_intersection)/len(expected_and_clustered_union)
        sensitivity = len(expected_and_clustered_intersection)/len(set_genomes_expected)

        # The traditional F-measure or balanced F-score (F1 score) is the harmonic mean of precision and recall
        # https://en.wikipedia.org/wiki/F1_score
        # print(precision)
        # print(sensitivity)
        F1_score = 2*(precision*sensitivity/(precision + sensitivity))

        dict_to_write = {"nb_cluster":nb_cluster,
                        "nb_cluster_with_poly":i,
                        "nb_clustered_taxid_with_poly":len(clustered_genome_with_poly),
                        "nb_expected_taxid":nb_genome_expected,
                        "complement_clustered_taxid_with_poly":len(clustered_genome_with_poly - set_genomes_expected),
                        "complement_expected_taxid": len(set_genomes_expected - clustered_genome_with_poly ),
                        "expected_and_clustered_intersection": len(expected_and_clustered_intersection),
                        "Evalue":evalue,
                        "coverage":coverage,
                        "inflation":inflation,
                        "precision": len(expected_and_clustered_intersection)/len(expected_and_clustered_union),
                        'sensitivity': len(expected_and_clustered_intersection)/len(set_genomes_expected),
                        "F1_score":F1_score
                        }

        csv_output.writerow(dict_to_write)
        # print(t)
        # print('nb of cluster with polyprotein annotated and expected', ie)
        # print('number of genome clustered with annotated and expected genome', len(clustered_genome_with_poly_exp))
        # t = setComparison(set(expected_tax_ids), clustered_genome_with_poly_exp)
        # print(t)

    fl_out.close()

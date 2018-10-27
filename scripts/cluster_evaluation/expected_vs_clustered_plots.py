# library
from pandas import DataFrame, read_csv
import matplotlib.pyplot as plt
import matplotlib
from matplotlib_venn import venn2, venn3, _venn3
from numpy import log10
import expected_taxon_evaluation as exp_eval
import logging, sys
import taxonomic_homogeneity_checking as tax_hom
import csv, os, re

def cluster_expected_vs_clustered(evalue, coverage, inflation):
    clustered_with_all_annotated, sensitivity, precision, F1_score = get_clustered_genome_set(cluster_dir, taxon, evalue, coverage, inflation, annotated_genomes_filtered, expected_tax_ids)
    if not clustered_with_all_annotated:
        return
    sets = [set_genomes_expected, clustered_with_all_annotated]
    labels=('Expected Genomes', 'Clustered Genomes')
    title = 'Expected Genomes vs Genomes Clustered with an expected and annotated genome'
    caption = f"\nClustering parameters: Evalue:{evalue} Coverage:{coverage} Inflation:{inflation}\nSensitivity:{round(sensitivity, 3)}  precision:{round(precision, 3)}  F-measure:{round(F1_score, 3)}"
    name = f'venn_expected_vs_clustered_ann_expected_E{evalue}_C{coverage}_I{inflation}.png'
    path_to_img = os.path.join(output_dir,name)

    getVenn2plot(sets, labels, title, caption,  path_to_img)

def getVenn2plot(sets, labels, title, caption,   path_to_img):
    fig = plt.figure(figsize=(7, 7))
    plt.xlim(0, 3)
    plt.ylim(0, 3)
    fig.text(0.15,.05, caption,fontsize=12)
    v = venn2(subsets = sets, set_labels = labels)


    v.get_patch_by_id('11').set_color('purple')
    v.get_patch_by_id('01').set_color('blue')

    plt.title(title)

    plt.savefig(path_to_img, bbox_inches='tight', transparent=False)
    plt.show()

def getVenn3plot(sets, labels, title,  path_to_img):

    v = venn3(subsets = sets, set_labels = labels)
    print(_venn3.compute_venn3_subsets(sets[0], sets[1], sets[2]))
    # v.get_patch_by_id('11').set_color('purple')
    # v.get_patch_by_id('01').set_color('blue')

    plt.title(title)

    plt.savefig(path_to_img, bbox_inches='tight')
    plt.show()

def get_expected_vs_annotated(expected_tax_ids, annotated_genomes ):
    print("NUMBER OF GENOMES WHERE WE EXPECT POLYPROTEIN",len(expected_tax_ids))

    print("NUMBER OF GENOMES ANNOTATED", len(annotated_genomes))
    t = exp_eval.setComparison(set(expected_tax_ids), set(annotated_genomes))
    print(t)
    expected_annotated = set(expected_tax_ids) &  set(annotated_genomes)
    print('NUMBER OF GENOMES EXPECTED AND ANNOTATED', len(expected_annotated))
    venn2(subsets = (t), set_labels = ('Expected Genomes', 'Annotated Genomes'))

    name='venn_expected_vs_annotated.png'
    path_to_img = os.path.join(output_dir,name)
    plt.savefig(path_to_img, bbox_inches='tight')
    plt.show()

def get_clustered_genome_set(cluster_dir, taxon, evalue, coverage, inflation, annotated_genomes, expected_tax_ids):


    file_name = f'{taxon}_evalue_{evalue}coverage{coverage}_I{inflation}.out'
    cluster_file = os.path.join(cluster_dir, file_name)
    print(cluster_file)
    print("Clustering file", file_name)
    if not os.path.exists(cluster_file):
        print('file', cluster_file)
        print("The file hasn't been found in", cluster_dir)
        return set(), None, None, None


    cluster_iter = tax_hom.clusterFileParser(cluster_file, annotated_genomes)

    clustered_genome_with_poly = set()
    i = 0
    nb_cluster = 0
    for cluster in cluster_iter:
        nb_cluster += 1
        if cluster['nb_polyprotein'] > 0:
            i += 1
            clustered_genome_with_poly |= set(cluster['genome_ids']) # add the two set

    t = exp_eval.setComparison(set(expected_tax_ids), clustered_genome_with_poly)
    expected_and_clustered_intersection = set_genomes_expected &  clustered_genome_with_poly
    expected_and_clustered_union = set_genomes_expected | clustered_genome_with_poly

    precision =len(expected_and_clustered_intersection)/len(clustered_genome_with_poly)
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
                    "precision": precision,
                    'sensitivity': sensitivity,
                    "F1_score":F1_score
                    }

    for k,v in dict_to_write.items():
        print(k, v)

    return clustered_genome_with_poly, sensitivity, precision, F1_score

if __name__ == '__main__':
    # csv_file = "results/clustering_evaluation/expected_taxon_cluster_evaluation.csv"
    # First way to call the 2 group Venn diagram:
    unclassified_term_file = "data/taxonomy/unclassified_terms.txt"
    taxonomy_file = "data/taxonomy/taxonomy_virus.txt"
    expected_taxon_file='etc/viruses_w_polyproteins.txt'
    stat_protein_file='results/stat_viral_protein/stat_proteins_Viruses.csv'

    cluster_dir='data/clustering_result/Viruses/clustering_parameter_variation'
    # cluster_files = (os.path.join(cluster_dir, f) for f in os.listdir(cluster_dir) if f.endswith('.out'))

    output_dir="results/figures/venn_diagram"

    taxon = "Viruses"

    try:
        inflation = sys.argv[1]
    except:
        inflation = "1_4"
    try:
        evalue = sys.argv[2]
    except:
        evalue = 1e-60
    try:
        coverage = sys.argv[3]
    except:
        coverage = 20

    expected_taxons = exp_eval.getExpectedTaxon(expected_taxon_file)
    expected_tax_ids = exp_eval.getExpectedGenomeList(taxonomy_file, expected_taxons)

    annotated_genomes = tax_hom.getAnnotatedProteinTaxId(stat_protein_file)

    nb_genome_expected = len(expected_tax_ids)
    set_genomes_expected = set(expected_tax_ids)

    expected_annotated_genome =  set(annotated_genomes) - set(annotated_genomes) | set_genomes_expected # union
    annotated_genomes_filtered = {k:v for k,v in annotated_genomes.items() if k in expected_annotated_genome}
    # will be used then to select cluster of interest
    unexpected_annotated_genome = {k:v for k,v in annotated_genomes.items() if k not in expected_annotated_genome}
    print("LEN unexpected_annotated genome", len(unexpected_annotated_genome))

    cluster_expected_vs_clustered(1e-50, 20, 2)
    cluster_expected_vs_clustered(1e-50, 60, '1_4')
    cluster_expected_vs_clustered(1e-60, 50, 3)
    cluster_expected_vs_clustered(1e-160, 50, 2)
    while 1:
        print("Current parameters:")
        print(f'Inflation:{inflation}\nevalue:{evalue}\ncoverage:{coverage}')
        print("\n 1.choose parameter\n",
        '2.Expected vs Annotated\n',
        "3.Venn 3 with expected, clustered from all annotated genome, clustered from expected annotated genome\n",
        "4.Expected vs Clustering\n",
        "5.from new refSeq\n",
        '0.Exit\n')
        rep = input()
        if rep == '1':
            inflation_c = input('inflation:')
            evalue_c = input('evalue:')
            coverage_c = input('coverage:')

            evalue = evalue_c if evalue_c != '' else evalue
            coverage = coverage_c if coverage_c != '' else coverage
            inflation = inflation_c if inflation_c != '' else inflation
            continue

        elif rep == '2':
            get_expected_vs_annotated(expected_tax_ids, annotated_genomes)

        elif rep == "0":
            break

        elif rep == '3':
            print("CLUSTERED FROM EXPECTED ANNOTATED")
            clustered_from_expected_annotated, sensitivity, precision, F1_score = get_clustered_genome_set(cluster_dir, taxon, evalue, coverage, inflation, annotated_genomes_filtered, expected_tax_ids)
            print("CLUSTERED FROM ALL ANNOTATED")
            clustered_from_all_annotated, sensitivity, precision, F1_score = get_clustered_genome_set(cluster_dir, taxon, evalue, coverage, inflation, annotated_genomes, expected_tax_ids)
            if not clustered_from_all_annotated or not clustered_from_expected_annotated:
                continue

            sets = [set_genomes_expected, clustered_from_all_annotated, clustered_from_expected_annotated]
            labels=('Expected Genomes', 'Clustered Genomes all', "Clustered Genome from expected_annotated")
            title = f"\nClustering parameters: Evalue:{evalue} Coverage:{coverage} Inflation:{inflation}"
            name = f'venn_expected_vs_clust_ann_vs_clust_expected_ann_E{evalue}_C{coverage}_I{inflation}.png'
            path_to_img = os.path.join(output_dir,name)
            # getVenn3plot(sets, labels, title,  path_to_img)
        elif rep == '4':
            cluster_expected_vs_clustered(evalue, coverage, inflation)

        elif rep == '5':
            new_cluster_dir = "data/clustering_result/Viruses/RefSeq_download_date_2018-07-21"

            clustered_from_expected_annotated, sensitivity, precision, F1_score = get_clustered_genome_set(new_cluster_dir,
                                                                        taxon,
                                                                        evalue,
                                                                        coverage,
                                                                        inflation,
                                                                        annotated_genomes_filtered,
                                                                        expected_tax_ids)
            sets = [set_genomes_expected, clustered_from_expected_annotated]
            if not clustered_from_expected_annotated:
                continue
            labels=('Expected Genomes', 'Clustered Genomes')
            title = 'Expected Genomes vs Genomes clustered with an expected and annotated genome'
            caption = f"Clustering parameters: Evalue:{evalue} Coverage:{coverage} Inflation:{inflation.replace('_', '.')}"
            name = f'RefSeq_download_date_2018-07-21_venn_expected_vs_clustered_E{evalue}_C{coverage}_I{inflation}.png'
            path_to_img = os.path.join(output_dir,name)
            getVenn2plot(sets, labels, title, caption,  path_to_img)

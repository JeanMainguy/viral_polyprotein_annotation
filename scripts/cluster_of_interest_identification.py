#!/usr/bin/env python3
import sys, csv

def getAnnotatedProteinTaxId(stat_protein_file, colname_filter = 'polyprotein_outline', filter_value='True'):
    annotated_genomes = {}
    with open(stat_protein_file) as csvfile:
        reader = csv.DictReader(csvfile, delimiter='\t')
        for row in reader:
            # print(row)
            if row[colname_filter] == filter_value:
                annotated_genomes.setdefault(row['taxon_id'], []).append(row['protein_id'])
                # print(row['taxon_id'], row['protein_id'])
    return annotated_genomes


def identify_annotated_cluster(annotated_genomes, cluster_file, cluster_of_interest_out_file):
    cluster_out = open(cluster_of_interest_out_file, 'w')

    with open(cluster_file, 'r') as fl:
        #665102|YP_003622544.1|381	1046572|YP_004958227.1|360	536084|YP_001960955.1|356
        poly_cluster_nb = 0
        for i, cluster_line in enumerate(fl):
            cluster_has_polyprotein = False

            cluster_elements = cluster_line.rstrip().split("\t")

            cluster_with_polyprotein = []
            for element in cluster_elements:
                genome_id, protein_id, length = element.split('|')

                if genome_id in annotated_genomes and protein_id in annotated_genomes[genome_id]:
                    cluster_has_polyprotein = True
                #     cluster_with_polyprotein.append(genome_id+'|'+protein_id+'|Polyprotein')
                # else:
                #     cluster_with_polyprotein.append(genome_id+'|'+protein_id+'|NONE')
            if cluster_has_polyprotein:
                # cluster_out.write('\t'.join(cluster_with_polyprotein)+"\n")
                cluster_out.write(cluster_line)
                poly_cluster_nb += 1
    cluster_out.close()

if __name__ == '__main__':
    cluster_file=  sys.argv[1]
    cluster_of_interest_out_file =  sys.argv[2]
    stat_protein_file =  sys.argv[3]

    annotated_genomes = getAnnotatedProteinTaxId(stat_protein_file)
    identify_annotated_cluster(annotated_genomes, cluster_file, cluster_of_interest_out_file)

#!/usr/bin/env python3

import taxonomy as tax
import viral_genome_classes as obj
import viruses_statistics as stat
import multiple_alignment_analysis as analysis
import visualisation_genome as fct_visu

import sys,collections

PROTEIN_SIZE_MAX=1200 # max polyprotein len is 20000 but only few protien are concerned..
SCREEN_SIZE = 250

def visualisation_protein(cds, segment, nb_line):
    display_len = int(nb_line)*SCREEN_SIZE
    conversion = (len(cds)/3)/display_len

    compatible_dico = collections.OrderedDict()

    compatible_dico['match'] = fct_visu.buildCompatibleGroup(set(cds.matchs))[::-1]
    compatible_dico['cds'] = [[cds]]
    compatible_dico['pep'] = fct_visu.buildCompatibleGroup(cds.peptides)
    compatible_dico['unannotated_region'] = fct_visu.buildCompatibleGroup(set(cds.unannotated_region))

    strings = fct_visu.get_final_strings(cds, 1, 1, compatible_dico)
    print(strings)


if __name__ == '__main__':

    try:
        tax_id_protein_id = sys.argv[1] # 11768|NP_047256.1
    except IndexError:
        tax_id_protein_id = "11768|NP_047256.1"

    try:
        fct_visu.SCREEN_SIZE = int(sys.argv[2])
        print('new SCREEN SIZE', fct_visu.SCREEN_SIZE)
    except IndexError:
        pass

    try:
        nb_line = sys.argv[3]
    except IndexError:
        nb_line = 1

    # we don't get result from a tuple because the length may included
    # or other info after the protein id sep by '|'
    tax_id_protein_id_splited = tax_id_protein_id.split('|')
    print(tax_id_protein_id_splited)
    taxon_id = tax_id_protein_id_splited[0]
    protein_id = tax_id_protein_id_splited[1]

    sp_treshold=90
    taxonomy_file="data/taxonomy/taxonomy_virus.txt"
    expected_file =  "data/taxonomy/polyprotein_expectation_by_taxon.csv"

    gff_file = 'data/interpro_results/interproscan-5.30-69.0/domains_viral_sequences.gff3'

    cds, segment = analysis.getCdsObject(taxon_id, protein_id, taxonomy_file,gff_file, sp_treshold)


    print('*-'*100)
    print(f'VISUALISATION OF THE PROTEIN {taxon_id}  FROM THE TAXON {taxon_id}')
    print('*-'*100)

    visualisation_protein(cds, segment, 1)

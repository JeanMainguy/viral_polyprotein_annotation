import taxonomy as tax

import os
import logging
import sys
import csv
import viral_protein_extraction as ext


def getProteinStat(cds, taxon_id, taxonomy):
    is_sub_protein = True if cds.parental_prot else False

    # Check if the cds has a signal peptide scheme
    # Meaning that it has from 1 to 2 pep annotations and one of it annotation start at the 5' of the cds

    potential_signal_P_first_pep_len = "None"

    cds_start = cds.start if cds.bp_obj.strand == 1 else cds.end
    potential_signal_P_first_pep_type = "None"
    for pep in cds.peptides:
        pep_start = pep.start if cds.bp_obj.strand == 1 else pep.end

        if cds_start == pep_start and 0 < len(cds.peptides) <= 2:
            potential_signal_P_first_pep_len = len(pep)
            potential_signal_P_first_pep_type = pep.bp_obj.type

    # Write dict info
    dict_info = {
        "First_node": taxonomy[1],
        "taxon_id": taxon_id,
        "strand": cds.bp_obj.strand,
        'protein_id': cds.protein_id,
        "polyprotein_outline": cds.polyprotein,
        "has_peptide": True if len(cds.peptides) else False,
        "nb_peptides": str(len(cds.peptides)),
        "nb_unannotated_part": str(len(cds.unannotated_region)),
        "nb_cleavage_sites": str(len(cds.cleavage_sites)),
        "is_sub_protein": is_sub_protein,
        "non_polyprotein_explanation": cds.non_polyprotein_explanation,
        "potential_signal_P_first_pep_len": potential_signal_P_first_pep_len,
        "taxonomy": str(taxonomy),
        "potential_signal_P_first_pep_type": potential_signal_P_first_pep_type,
        "len_protein": len(cds)
    }
    return dict_info


def writeGenomeStat(taxon_id, nb_cds, nb_peptide, handle_stat_genome, taxonomy, nb_polyprotein):
    #  ["taxon_id", "nb_protein", "has_peptide", "nb_peptides"]
    has_peptide = "TRUE" if nb_peptide else 'FALSE'
    has_polyprotein = True if nb_polyprotein else False
    line = [taxon_id, str(nb_cds), str(has_polyprotein), str(nb_polyprotein),
            has_peptide, str(nb_peptide), ";".join(taxonomy)]
    handle_stat_genome.write("\t".join(line)+"\n")


def initiateBasicStatFile(taxon, output_dir):
    taxon = taxon.replace(',', '').replace(' ', '_')

    stat_file_prot = "stat_proteins_{}.csv".format(taxon)
    stat_file_genome = "stat_genomes_{}.csv".format(taxon)

    handle_stat_prot = open(os.path.join(output_dir, stat_file_prot), "w")
    # segment.taxon_id, cds.protein_id, has_peptide, len(cds.peptides), len(cds.cleavage_sites), is_sub_protein
    protein_header = ["taxon_id",
                      "protein_id",
                      "polyprotein_outline",
                      "has_peptide",
                      "nb_peptides",
                      "nb_unannotated_part",
                      "nb_cleavage_sites",
                      "is_sub_protein",
                      "non_polyprotein_explanation",
                      "potential_signal_P_first_pep_len",
                      "potential_signal_P_first_pep_type",
                      "len_protein",
                      "strand",
                      "taxonomy",
                      "First_node"]

    csv_writer_protein = csv.DictWriter(handle_stat_prot, fieldnames=protein_header, delimiter='\t')
    csv_writer_protein.writeheader()

    # #PEPTIDE
    # handle_stat_pep = open(os.path.join(output_dir, stat_file_peptides), "w")
    # #segment.taxon_id, cds.protein_id, has_peptide, len(cds.peptides), len(cds.cleavage_sites), is_sub_protein
    # pep_header = [  "taxon_id",
    #                 "peptide_id",
    #                 "feature_type",
    #                 "strand",
    #                 "len",
    #                 "nb_proteins_annotated_by_pep",
    #                 "protein_ids",
    #                 "polyprotein_outline",
    #                 "nb_protein_annotated_by_pep",
    #                 "cover_by_domain",
    #                 "overlapped_by_domain", # are domains overlapping the peptide annotation
    #                 "length_covered_by_fully_included_domains", # length covered by domains that are fully included
    #                 "length_covered_by_all_domains", # length covered by domains that are fully included or not
    #                 "First_node",
    #                 "taxonomy" ]

    # csv_writer_peptides = csv.DictWriter(handle_stat_pep, fieldnames=pep_header, delimiter='\t')
    # csv_writer_peptides.writeheader()

    handle_stat_genome = open(os.path.join(output_dir, stat_file_genome), "w")
    # taxon
    header = ["taxon_id",
              "nb_protein",
              "has_annotated_poly",
              "nb_polyprotein",
              "has_peptide",
              "nb_peptides",
              "taxonomy"]
    handle_stat_genome.write('\t'.join(header)+'\n')

    files_to_close = [handle_stat_genome, handle_stat_prot]

    writer_stat_dict = {"protein": csv_writer_protein,
                        "genome": handle_stat_genome}

    return writer_stat_dict, files_to_close


if __name__ == '__main__':

    # logging.basicConfig(filename='log/viral_protein_statistic.log',level=logging.INFO)
    logging.basicConfig(level=logging.INFO)

    try:
        taxon = sys.argv[1]  # if not given then we don't compute any stat
    except IndexError:
        taxon = "Viruses"
    try:
        stat_output_dir = sys.argv[2]
    except IndexError:
        stat_output_dir = 'test'

    taxonomy_file = "data/taxonomy/taxonomy_virus.txt"
    sp_treshold = 90

    gbff_iter = tax.getAllRefseqFromTaxon(taxon, taxonomy_file)

    taxon = taxon.replace(',', '').replace(' ', '_')

    writer_stat_dict, files_to_close = initiateBasicStatFile(taxon, stat_output_dir)

    handle_prot = False  # no protein extraction

    i = None
    for i, gb_dico in enumerate(gbff_iter):
        if (i+1) % 1000 == 0:
            print(i+1, 'genomes treated')
        gb_file = gb_dico['gb_file']
        genetic_code = gb_dico['genetic_code']
        taxon_id = gb_dico['taxon_id']
        print(taxon_id)
        ext.extractAndStat(gb_file, handle_prot, writer_stat_dict,
                           genetic_code, taxon_id, sp_treshold)

    if i is None:
        print("No genome have been found with taxon:", taxon)
    else:
        print('Analysis completed for ', i+1, "genomes")

    for f in files_to_close:
        f.close()

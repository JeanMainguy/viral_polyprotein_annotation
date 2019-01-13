#!/usr/bin/env python3
import taxonomy as tax
import viral_genome_classes as obj
import viruses_statistics as stat
import visualisation_protein as visu

import os
import logging
import sys
from Bio import SeqIO
# from Bio.SeqRecord import SeqRecord
# from Bio.Seq import Seq
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter


def extract_cleavage_site_sequences(cds_list, extract_window):
    half_window = int(extract_window/2)
    print("kalf wind", half_window)
    for cds in cds_list:
        # organism = cds.segment.organis
        genetic_code = cds.genetic_code

        protein_sequence = cds.getSequenceAA(cds.segment.record, genetic_code)

        for cs in cds.cleavage_sites:
            seq = protein_sequence[cs.start_aa(cds)-half_window:cs.end_aa(cds)+half_window-1]

            header = f'{cds.segment.taxon_id}|{cds.protein_id}|Cleavage_site_{cs.number}'
            # description = f"from {cs.start_aa(cds)} to {cs.end_aa(cds)}"

            print(f'>{header}')
            print(seq)


def extractAndStat(gb_file, handle_prot, writer_stat_dict, genetic_code, taxon_id, sp_treshold):
    nb_peptide = 0
    nb_cds = 0

    genome = obj.gb_file_parser(gb_file, taxon_id, sp_treshold)

    for i, segment in enumerate(genome.segments):
        taxonomy = segment.record.annotations['taxonomy']
        nb_cds += len(segment.cds)
        nb_polyprotein = 0

        for p, cds in enumerate(segment.cds):
            key = str(int(len(cds)/3))
            header = '|'.join([taxon_id, cds.protein_id, key])
            if cds.polyprotein:
                nb_polyprotein += 1
                polypro_list_fl.write(header+'\n')

            # WRITE PROTEIN SEQUENCE
            if handle_prot:  # if we want to write the sequences it the writer otherwise False
                seq_to_write = cds.getSequenceRecord(segment.organism, header, genetic_code)
                SeqIO.write(seq_to_write, handle_prot, "fasta")

            # identification of poorly annotated cds
            if len(cds.unannotated_region) > 1 and irrelevant_annotatation_file:
                string_info = f'\nThe cds {cds.protein_id} from taxon:{cds.segment.taxon_id} seems poorly annotated\n'
                string_info += f'It has {len(cds.unannotated_region)} unnanotated part\n'
                string_info += visu.visualisation_protein(cds, 1, len(cds))[:-4]
                string_info += "\n" + ":/"*50 + "\n"
                irrelevant_annotatation_file.write(string_info)

                # BLACK LIST
                black_list_writer.write(f"{taxon_id}|{cds.protein_id}\n")

            # WRITE PROTEIN STAT
            if writer_stat_dict:  # may be False if we don't want to write stat
                info_dict = stat.getProteinStat(cds, taxon_id, taxonomy)
                writer_stat_dict['protein'].writerow(info_dict)

                nb_peptide += len(cds.peptides)

    if writer_stat_dict:
        # WRITE GENOME STAT
        stat.writeGenomeStat(taxon_id, nb_cds, nb_peptide,
                             writer_stat_dict['genome'], taxonomy, nb_polyprotein)


def parse_arguments():

    parser = ArgumentParser(description="Extract and stat viral proteins",
                            formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument("taxon", type=str,
                        help="")
    parser.add_argument("genome_index_file", type=str,
                        help="Index file of the genomes")

    parser.add_argument("protein_db_file", type=str,
                        help="file where fasta file of protein sequences is stored")

    parser.add_argument("--sp_treshold", type=int, default=90,
                        help="Signal peptide length threshold.")

    parser.add_argument('--stat_output_dir', type=str, default=False,
                        help="stat output dir if empty no stat are recorded")

    parser.add_argument("--polyprotein_list_file", type=str, default="annotated_polyprotein_list.txt",
                        help="list of polyproteins")

    parser.add_argument("--irrelevant_annotation_list", type=str, default="irrelevant_annotation_list.txt",
                        help="list of irrelevant polyprotein annotations")

    args = parser.parse_args()
    return args


if __name__ == '__main__':

    try:
        logging.basicConfig(filename='log/viral_protein_extraction_and_stat.log',
                            level=logging.INFO)
    except FileNotFoundError:
        logging.basicConfig(filename='viral_protein_extraction_and_stat.log',
                            level=logging.INFO)

    args = parse_arguments()

    taxon = args.taxon
    protein_db = args.protein_db_file

    taxonomy_file = args.genome_index_file

    sp_treshold = args.sp_treshold

    stat_output_dir = args.stat_output_dir

    irrelevant_annotation_list_fl = args.irrelevant_annotation_list
    polyprotein_list_file = args.polyprotein_list_file
    # try:
    #     stat_output_dir = sys.argv[5]  # if not given then we don't compute any stat
    #     print(f'Write genomes, proteins and peptides in {stat_output_dir} statistics')
    # except IndexError:
    #     print("No statistics...")
    #     stat_output_dir = False
    # try:
    #     irrelevant_annotation_list_fl = sys.argv[6]
    # except IndexError:
    #     irrelevant_annotation_list_fl = "test/irrelevant_annotation_list.txt"
    #
    # try:
    #     polyprotein_list_file = sys.argv[7]
    # except IndexError:
    #     polyprotein_list_file = "test/annotated_polyprotein_list.txt"

    gbff_iter = tax.getAllRefseqFromTaxon(taxon, taxonomy_file)

    taxon = taxon.replace(',', '').replace(' ', '_')

    if stat_output_dir:
        writer_stat_dict, files_to_close = stat.initiateBasicStatFile(taxon, stat_output_dir)
        irrelevant_annotatation_file = open(os.path.join(
            stat_output_dir, f'irrelevant_annotation_identification_{taxon}.txt'), "w")
        files_to_close.append(irrelevant_annotatation_file)
    else:
        files_to_close = []
        writer_stat_dict = False
        irrelevant_annotatation_file = False

    # protein_db = "{}_protein_db.faa".format(taxon)

    handle_prot = open(protein_db, "w")
    files_to_close.append(handle_prot)

    print("Creation of a list of irrelevant mat_peptide annotations")
    black_list_writer = open(irrelevant_annotation_list_fl, "w")
    files_to_close.append(black_list_writer)

    print("Creation of a list of annotated polyproteins")
    polypro_list_fl = open(polyprotein_list_file, "w")
    files_to_close.append(polypro_list_fl)

    i = None
    for i, gb_dico in enumerate(gbff_iter):
        if (i+1) % 1000 == 0:
            print(i+1, 'genomes analysed')
        # print(gb_dico)
        gb_file = gb_dico['gb_file']
        genetic_code = gb_dico['genetic_code']
        taxon_id = gb_dico['taxon_id']

        extractAndStat(gb_file, handle_prot, writer_stat_dict, genetic_code, taxon_id, sp_treshold)

    if i is None:
        raise NameError("No genome have been found with taxon:", taxon)
    else:
        print(f'Analysis completed for {i+1} genomes')

    for f in files_to_close:
        f.close()

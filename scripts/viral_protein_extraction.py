#!/usr/bin/env python3
import taxonomy as tax
import viral_genome_classes as obj
import viruses_statistics as stat
import parser_interpro_results as do
import visualisation_protein as visu

import os
import gzip
import logging
import sys
from Bio import SeqIO
# from Bio.SeqRecord import SeqRecord
# from Bio.Seq import Seq


def extract_cleavage_site_sequences(cds_list, extract_window):
    half_window = int(extract_window/2)
    print("kalf wind", half_window)
    for cds in cds_list:
        organism = cds.segment.organism
        genetic_code = cds.genetic_code

        protein_sequence = cds.getSequenceAA(cds.segment.record, genetic_code)

        for cs in cds.cleavage_sites:
            seq = protein_sequence[cs.start_aa(cds)-half_window:cs.end_aa(cds)+half_window-1]

            header = f'{cds.segment.taxon_id}|{cds.protein_id}|Cleavage_site_{cs.number}'
            description = f"from {cs.start_aa(cds)} to {cs.end_aa(cds)}"

            print(f'>{header}')
            print(seq)

    input()

    # seq_to_write = SeqRecord(Seq(seq, generic_protein) ,id=header, description = description)
    #
    # SeqIO.write(seq_to_write, handle_prot,"fasta")


def extractAndStat(gb_file, handle_prot, writer_stat_dict, genetic_code, taxon_id, sp_treshold):
    nb_peptide = 0
    nb_cds = 0

    genome = obj.gb_file_parser(gb_file, taxon_id, sp_treshold)

    for i, segment in enumerate(genome.segments):
        taxonomy = segment.record.annotations['taxonomy']
        nb_cds += len(segment.cds)
        nb_polyprotein = 0

        for p, cds in enumerate(segment.cds):
            if cds.polyprotein:
                nb_polyprotein += 1

            # WRITE PROTEIN SEQUENCE
            if handle_prot:  # if we want to write the sequences it the writer otherwise False
                key = str(int(len(cds)/3))
                header = '|'.join([taxon_id, cds.protein_id, key])
                seq_to_write = cds.getSequenceRecord(segment.organism, header, genetic_code)
                SeqIO.write(seq_to_write, handle_prot, "fasta")

            # identification of poorly annotated cds
            if len(cds.unannotated_region) > 1:
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

        #
        # # WRITE PEPTIDE STAT
        # if writer_stat_dict:
        #     for pep in segment.peptides:
        #         pep_info_dict = stat.getPepStat(pep, taxon_id, taxonomy)
        #         writer_stat_dict['peptide'].writerow(pep_info_dict)

    if writer_stat_dict:
        # WRITE GENOME STAT
        stat.writeGenomeStat(taxon_id, nb_cds, nb_peptide,
                             writer_stat_dict['genome'], taxonomy, nb_polyprotein)

        # # WRITE DOMAIN STAT
        # domain_header = writer_stat_dict['domain'].fieldnames
        # for domain in genome.matchs:
        #     info_dict = stat.getDomainStat(taxon_id, domain, domain_header)
        #     writer_stat_dict['domain'].writerow(info_dict)
        #
        # # WRITE CLEAVAGE SITE STAT
        # cs_header = writer_stat_dict['cleavage_site'].fieldnames
        # for cs in segment.cleavage_sites:
        #     # print(cs)
        #     info_dict = stat.getCleavageSiteStat(taxon_id, cs, cs_header)
        #     writer_stat_dict['cleavage_site'].writerow(info_dict)


if __name__ == '__main__':
    from time import clock
    START_TIME = clock()

    logging.basicConfig(filename='log/viral_protein_extraction_and_stat.log',
                        level=logging.INFO)

    taxon = sys.argv[1]
    output_dir = sys.argv[2]

    taxonomy_file = sys.argv[3]
    sp_treshold = int(sys.argv[4])

    try:
        stat_output_dir = sys.argv[5]  # if not given then we don't compute any stat
        print(f'Write genomes, proteins and peptides in {stat_output_dir} statistics')
    except:
        print("No statistics...")
        stat_output_dir = False
    try:
        irrelevant_annotation_list_fl = sys.argv[6]
    except:
        irrelevant_annotation_list_fl = "test/irrelevant_annotation_list.txt"

    gbff_iter = tax.getAllRefseqFromTaxon(taxon, taxonomy_file)

    taxon = taxon.replace(',', '').replace(' ', '_')

    if stat_output_dir:
        writer_stat_dict, files_to_close = stat.initiateBasicStatFile(taxon, stat_output_dir)
    else:
        files_to_close = []
        writer_stat_dict = False

    protein_db = "{}_protein_db.faa".format(taxon)

    handle_prot = open(os.path.join(output_dir, protein_db), "w")
    files_to_close.append(handle_prot)

    irrelevant_annotatation_file = open(os.path.join(
        stat_output_dir, f'irrelevant_annotation_identification_{taxon}.txt'), "w")
    files_to_close.append(handle_prot)

    print("Creation of a list of irrelevant mat_peptide annotations")
    black_list_writer = open(irrelevant_annotation_list_fl, "w")
    files_to_close.append(black_list_writer)

    i = None
    for i, gb_dico in enumerate(gbff_iter):
        if (i+1) % 1000 == 0:
            print(i+1, 'genomes analysed')
        # print(gb_dico)
        gb_file = gb_dico['gb_file']
        genetic_code = gb_dico['genetic_code']
        taxon_id = gb_dico['taxon_id']

        extractAndStat(gb_file, handle_prot, writer_stat_dict, genetic_code, taxon_id, sp_treshold)

    if i == None:
        raise NameError("No genome have been found with taxon:", taxon)
    else:
        print(f'Analysis completed for {i+1} genomes')

    for f in files_to_close:
        f.close()

    print('\nPROGRAM RUN TIME:%6.2f' % (clock()-START_TIME), 'seconds.')
    print('-'*20)

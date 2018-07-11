#!/usr/bin/env python3

import taxonomy as tax
import viral_genome_classes as obj
import viruses_statistics as stat


import os, gzip, logging
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import sys, csv, re
from operator import attrgetter



def extractProteins(gb_file, handle_prot, writer_dict, genetic_code, sp_treshold):
    nb_peptide = 0
    nb_cds = 0
    genome = obj.Genome( gb_file)
    # print(gb_file)
    with gzip.open(gb_file, "rt") as handle:
        for i, record in enumerate(SeqIO.parse(handle, "genbank")):


            segment = obj.Segment(record, gb_file)
            genome.segments.append(segment)

            segment.getMatpeptidesAndPolyproteins()
            # print(segment.taxon_id)

            if not segment.peptides: # if no peptide annotation we don't need to do the next step of the loop
                continue


            segment.checkPeptideRedundancy() #remove the redundant peptide
            segment.checkSubPeptides()
            segment.associatePepWithProt(sp_treshold)

            # segment.checkForSlippage()
            segment.identifySubProtein()

            segment.getCleavageSites()


    for i, segment in enumerate(genome.segments):
        taxon_id = segment.taxon_id
        taxonomy = record.annotations['taxonomy']
        nb_cds += len(segment.cds)

        for p, cds in enumerate(segment.cds):
            # annotation = "Peptide" if cds.polyprotein else None
            # key = '{}|{}'.format(int(len(cds)/3), annotation)

            key=str(int(len(cds)/3))
            header = '|'.join([taxon_id, cds.protein_id, key])

            seq_to_write = cds.getSequenceRecord(segment.organism, header, segment.record, genetic_code)

            SeqIO.write(seq_to_write, handle_prot,"fasta")

            if  writer_dict: # may be False if we don't want to write stat
                info_dict = stat.getProteinStat(cds, segment.taxon_id, taxonomy)
                writer_dict['protein'].writerow(info_dict)

                nb_peptide += len(cds.peptides)

        if writer_dict:
            for pep in segment.peptides:
                pep_info_dict = stat.getPepStat(pep, segment.taxon_id, taxonomy)
                writer_dict['peptide'].writerow(pep_info_dict)


    if writer_dict:
        stat.writeGenomeStat(taxon_id, nb_cds, nb_peptide, writer_dict['genome'], taxonomy)


if __name__ == '__main__':
    from time import clock; START_TIME = clock()

    logging.basicConfig(filename='log/viral_protein_extraction.log',level=logging.INFO)

    taxon = sys.argv[1]
    output_dir = sys.argv[2]

    taxonomy_file = sys.argv[3]
    sp_treshold = int(sys.argv[4])
    try:
        stat_output_dir = sys.argv[5] # if not given then we don't compute any stat
    except:
        stat_output_dir = False

    gbff_iter = tax.getAllRefseqFromTaxon(taxon, taxonomy_file)

    taxon = taxon.replace(',', '').replace(' ', '_')

    if stat_output_dir:
        csv_writer_peptides, csv_writer_protein, handle_stat_genome, files_to_close = stat.initiateStatFile(taxon, stat_output_dir)
        writer_dict = {"peptide":csv_writer_peptides, "protein":csv_writer_protein, "genome":handle_stat_genome}
    else:
        files_to_close = []
        writer_dict = False


    protein_db = "{}_protein_db.faa".format(taxon)

    handle_prot = open(os.path.join(output_dir,protein_db), "w")



    i=None
    for i, gb_dico in enumerate(gbff_iter):
        # print(gb_dico)
        gb_file = gb_dico['gb_file']
        genetic_code = gb_dico['genetic_code']

        extractProteins(gb_file, handle_prot, writer_dict, genetic_code, sp_treshold)
        # if i%100 == 0:
        #     print(i)
        if (i+1)%1000 == 0:
            print(i+1, 'genomes analysed')
    if i == None:
        # print("No genome have been found with taxon:",taxon)
        raise NameError("No genome have been found with taxon:",taxon)
    else:
        print('Analysis completed for ', i+1, "genomes")


    handle_prot.close()
    for f in files_to_close:
        f.close()
        f.close()


    print('\nPROGRAM RUN TIME:%6.2f'%(clock()-START_TIME), 'seconds.');
    print('-'*20)

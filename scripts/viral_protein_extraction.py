#!/usr/bin/env python3

import taxonomy as tax
import viral_genome_classes as obj
import viruses_statistics as stat


import os, gzip, logging, sys
from Bio import SeqIO
# from Bio.SeqRecord import SeqRecord
# from Bio.Seq import Seq




def extractAndStat(gb_file, handle_prot, writer_stat_dict, genetic_code, taxon_id, sp_treshold):
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
            segment.associatePepWithProt()

            # segment.checkForSlippage()
            segment.identifySubProtein()

            segment.getCleavageSites()
            segment.identifyAnnotatedPolyproteins(sp_treshold)


    for i, segment in enumerate(genome.segments):
        taxonomy = record.annotations['taxonomy']
        nb_cds += len(segment.cds)
        nb_polyprotein = 0

        for p, cds in enumerate(segment.cds):
            if cds.polyprotein:
                nb_polyprotein += 1

            #WRITE PROTEIN SEQUENCE
            if handle_prot: # if we want to write the sequences it the writer otherwise False
                key=str(int(len(cds)/3))
                header = '|'.join([taxon_id, cds.protein_id, key])
                seq_to_write = cds.getSequenceRecord(segment.organism, header, segment.record, genetic_code)
                SeqIO.write(seq_to_write, handle_prot,"fasta")


            #WRITE PROTEIN STAT
            if  writer_stat_dict: # may be False if we don't want to write stat
                info_dict = stat.getProteinStat(cds, taxon_id, taxonomy)
                writer_stat_dict['protein'].writerow(info_dict)

                nb_peptide += len(cds.peptides)
        # WRITE PEPTIDE STAT
        if writer_stat_dict:
            for pep in segment.peptides:
                pep_info_dict = stat.getPepStat(pep, taxon_id, taxonomy)
                writer_stat_dict['peptide'].writerow(pep_info_dict)

    # WRITE GENOME STAT
    if writer_stat_dict:
        stat.writeGenomeStat(taxon_id, nb_cds, nb_peptide, writer_stat_dict['genome'], taxonomy, nb_polyprotein)


if __name__ == '__main__':
    from time import clock; START_TIME = clock()

    logging.basicConfig(filename='log/viral_protein_extraction_and_stat.log',level=logging.INFO)

    taxon = sys.argv[1]
    output_dir = sys.argv[2]

    taxonomy_file = sys.argv[3]
    sp_treshold = int(sys.argv[4])
    try:
        stat_output_dir = sys.argv[5] # if not given then we don't compute any stat
        print(f'Write genomes, proteins and peptides in {stat_output_dir} statistics')
    except:
        print("No statistics...")
        stat_output_dir = False

    gbff_iter = tax.getAllRefseqFromTaxon(taxon, taxonomy_file)

    taxon = taxon.replace(',', '').replace(' ', '_')

    if stat_output_dir:
        csv_writer_peptides, csv_writer_protein, handle_stat_genome, files_to_close = stat.initiateStatFile(taxon, stat_output_dir)
        writer_stat_dict = {"peptide":csv_writer_peptides, "protein":csv_writer_protein, "genome":handle_stat_genome}
    else:
        files_to_close = []
        writer_stat_dict = False


    protein_db = "{}_protein_db.faa".format(taxon)

    handle_prot = open(os.path.join(output_dir,protein_db), "w")
    files_to_close.append(handle_prot)


    i=None
    for i, gb_dico in enumerate(gbff_iter):
        if (i+1)%1000 == 0:
            print(i+1, 'genomes analysed')
        # print(gb_dico)
        gb_file = gb_dico['gb_file']
        genetic_code = gb_dico['genetic_code']
        taxon_id = gb_dico['taxon_id']

        extractAndStat(gb_file, handle_prot, writer_stat_dict, genetic_code, taxon_id, sp_treshold)

    if i == None:
        raise NameError("No genome have been found with taxon:",taxon)
    else:
        print('Analysis completed for ', i+1, "genomes")

    for f in files_to_close:
        f.close()


    print('\nPROGRAM RUN TIME:%6.2f'%(clock()-START_TIME), 'seconds.');
    print('-'*20)

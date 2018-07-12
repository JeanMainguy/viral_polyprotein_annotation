import taxonomy as tax
import viral_genome_classes as obj

import os, gzip, logging
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import sys, csv, re
from Bio.Alphabet import generic_protein
from Bio.SeqFeature import SeqFeature, FeatureLocation
from operator import attrgetter
import viral_protein_extraction as ext


def getPepStat(pep, taxon_id, taxonomy):
    #Check if the cds has a signal peptide scheme
    #Meaning that it has from 1 to 2 pep annotations and one of it annotation start at the 5' of the cds

    # start/end position inverse if strand -1
    signal_peptide_outline = False
    start_as_prot = False
    for cds in pep.polyproteins:
        cds_start = cds.start if cds.bp_obj.strand == 1 else  cds.end
        pep_start = pep.start if cds.bp_obj.strand == 1 else pep.end

        if cds_start == pep_start:
            start_as_prot = True
            if 0 < len(cds.peptides) <= 2:
                signal_peptide_outline = True

    peptide_id = "Unknown" if 'protein_id' not in pep.bp_obj.qualifiers else pep.bp_obj.qualifiers['protein_id'][0]
    #Write dict info
    dict_info = {
        "First_node":taxonomy[1],
        "taxon_id":taxon_id,
        "peptide_id":peptide_id,
        "feature_type":pep.bp_obj.type ,
        "strand":pep.bp_obj.strand,
        "start_as_prot":start_as_prot,
        "len":len(pep),
        "signal_peptide_outline":signal_peptide_outline,
        "inclued_in_nb_prot":len(pep.polyproteins),
        "taxonomy":taxonomy
        }

    return dict_info

def getProteinStat(cds, taxon_id, taxonomy):
    has_peptide = True if cds.peptides else False
    is_sub_protein = True if cds.parental_prot else False

    #Check if the cds has a signal peptide scheme
    #Meaning that it has from 1 to 2 pep annotations and one of it annotation start at the 5' of the cds

    # start/end position inverse if strand -1
    non_polyprotein_explanation = cds.non_polyprotein_explanation
    len_first_peptide = "None"

    cds_start = cds.start if cds.bp_obj.strand == 1 else cds.end
    first_pep_type = "None"
    for pep in cds.peptides:
        pep_start = pep.start if cds.bp_obj.strand == 1 else pep.end

        if cds_start == pep_start and  0 < len(cds.peptides) <= 2:
            len_first_peptide = len(pep)
            first_pep_type = pep.bp_obj.type

    #Write dict info
    dict_info = {
        "First_node":taxonomy[1],
        "taxon_id":taxon_id,
        "strand":cds.bp_obj.strand,
        'protein_id':cds.protein_id,
        "polyprotein_outline":cds.polyprotein,
        "has_peptide":True if len(cds.peptides) else False,
        "nb_peptides":str(len(cds.peptides)),
        "nb_unannotated_part":str(len(cds.unannotated_region)),
        "nb_cleavage_sites":str(len(cds.cleavage_sites)),
        "is_sub_protein":is_sub_protein,
        "non_polyprotein_explanation":non_polyprotein_explanation,
        "len_first_peptide":len_first_peptide,
        "taxonomy":str(taxonomy),
        "first_pep_type":first_pep_type,
	    "len_protein":len(cds)
            }
    return dict_info



def writeGenomeStat(taxon_id, nb_cds, nb_peptide, handle_stat_genome, taxonomy, nb_polyprotein):
    #  ["taxon_id", "nb_protein", "has_peptide", "nb_peptides"]
    has_peptide = "TRUE" if nb_peptide else 'FALSE'
    has_polyprotein = True if nb_polyprotein else False
    line = [taxon_id,str(nb_cds), str(has_polyprotein), str(nb_polyprotein),  has_peptide, str(nb_peptide), ";".join(taxonomy)]
    handle_stat_genome.write("\t".join(line)+"\n")



def initiateStatFile(taxon, output_dir ):
    taxon = taxon.replace(',', '').replace(' ', '_')


    stat_file_prot = "stat_proteins_{}.csv".format(taxon)
    stat_file_genome = "stat_genomes_{}.csv".format(taxon)
    stat_file_peptides = 'stat_peptides_{}.csv'.format(taxon)

    handle_stat_prot = open(os.path.join(output_dir, stat_file_prot), "w")
    #segment.taxon_id, cds.protein_id, has_peptide, len(cds.peptides), len(cds.cleavage_sites), is_sub_protein
    protein_header = ["taxon_id",
                    "protein_id",
                    "polyprotein_outline",
                    "has_peptide",
                    "nb_peptides",
                    "nb_unannotated_part",
                    "nb_cleavage_sites",
                    "is_sub_protein",
                    "non_polyprotein_explanation",
                    "len_first_peptide",
                    "first_pep_type",
		            "len_protein",
                    "strand",
                    "taxonomy",
                    "First_node"]

    csv_writer_protein = csv.DictWriter(handle_stat_prot, fieldnames=protein_header, delimiter='\t')
    csv_writer_protein.writeheader()

    #PEPTIDE
    handle_stat_pep = open(os.path.join(output_dir, stat_file_peptides), "w")
    #segment.taxon_id, cds.protein_id, has_peptide, len(cds.peptides), len(cds.cleavage_sites), is_sub_protein
    pep_header = [  "taxon_id",
                    "peptide_id",
                    "feature_type",
                    "strand",
                    "start_as_prot",
                    "len",
                    "nb_peptides",
                    "signal_peptide_outline",
                    "inclued_in_nb_prot",
                    "First_node",
                    "taxonomy" ]

    csv_writer_peptides = csv.DictWriter(handle_stat_pep, fieldnames=pep_header, delimiter='\t')
    csv_writer_peptides.writeheader()


    handle_stat_genome = open(os.path.join(output_dir, stat_file_genome), "w")
    #taxon
    header = ["taxon_id","nb_protein", "has_annotated_poly","nb_polyprotein",  "has_peptide", "nb_peptides", "taxonomy"]
    handle_stat_genome.write('\t'.join(header)+'\n')

    files_to_close = [handle_stat_genome, handle_stat_prot, handle_stat_pep]

    return csv_writer_peptides, csv_writer_protein, handle_stat_genome, files_to_close

if __name__ == '__main__':

    # logging.basicConfig(filename='log/viral_protein_statistic.log',level=logging.INFO)
    logging.basicConfig(level=logging.INFO)

    try:
        taxon = sys.argv[1] # if not given then we don't compute any stat
    except:
        taxon = "Viruses"
    try:
        stat_output_dir = sys.argv[2] # if not given then we don't compute any stat
    except:
        stat_output_dir = 'test'

    taxonomy_file ="data/taxonomy/taxonomy_virus.txt"
    sp_treshold = 90

    gbff_iter = tax.getAllRefseqFromTaxon(taxon, taxonomy_file)

    taxon = taxon.replace(',', '').replace(' ', '_')

    csv_writer_peptides, csv_writer_protein, handle_stat_genome, files_to_close = initiateStatFile(taxon, stat_output_dir )
    writer_stat_dict = {"peptide":csv_writer_peptides, "protein":csv_writer_protein, "genome":handle_stat_genome}
    handle_prot = False # no protein extraction

    i=None
    for i, gb_dico in enumerate(gbff_iter):
        if (i+1)%1000 == 0:
            print(i+1, 'genomes treated')
        gb_file = gb_dico['gb_file']
        genetic_code = gb_dico['genetic_code']
        taxon_id = gb_dico['taxon_id']

        ext.extractAndStat(gb_file, handle_prot, writer_stat_dict, genetic_code, taxon_id, sp_treshold)

    if i == None:
        print("No genome have been found with taxon:",taxon)
    else:
        print('Analysis completed for ', i+1, "genomes")

    for f in files_to_close:
        f.close()

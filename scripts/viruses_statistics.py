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

def get_csv_dico(self, header):
    # get dico info to write csv. header need to follow method or attribute of the obj
    dico = {}
    for attribute_name in header:
        if hasattr(self, attribute_name):
            attribute = getattr(self, attribute_name)

            if callable(attribute):
                dico[attribute_name] = attribute()
            else:
                dico[attribute_name] = attribute

    return dico

def getNonOverlappingCoveragePositionsPairs(positions):
    """
    listposition is a list of tuple (start, end) that potentially  overlaped each other
    We want to find pair of positions that represent the coverage...
    """
    positions = sorted(positions, key=lambda x: x[0], reverse=False) # sort positions by start
    try:
        s, e = positions.pop()
    except IndexError:
        return []

    non_overlapping_pairs = []

    for new_s, new_e in positions: # positions is sorted by starts
        if s < new_s < new_e: # the two pair are overlapping
            e = max(new_e, e)
        else: # the two pait are not overlapping then s and e are saved and replace by new s et new e
            non_overlapping_pairs.append((s, e))
            s = new_s
            e = new_e

    non_overlapping_pairs.append((s, e))
    return non_overlapping_pairs


def getPepStat(pep, taxon_id, taxonomy):
    #Check if the cds has a signal peptide scheme
    #Meaning that it has from 1 to 2 pep annotations and one of it annotation start at the 5' of the cds

    # start/end position inverse if strand -1
    polyprotein_outline = any((True for p in pep.polyproteins if p.polyprotein))

    fully_included_positions = getNonOverlappingCoveragePositionsPairs(pep.fully_included_domains.values())
    len_fully_coverage = sum((end - start+1 for start, end in fully_included_positions))

    included_positions = getNonOverlappingCoveragePositionsPairs(list(pep.partially_included_domains.values())+list(pep.fully_included_domains.values()))
    len_all_coverage = sum((end - start+1 for start, end in included_positions))


    peptide_id = "Unknown" if 'protein_id' not in pep.bp_obj.qualifiers else pep.bp_obj.qualifiers['protein_id'][0]

    #Write dict info
    dict_info = {
        "First_node":taxonomy[1],
        "taxon_id":taxon_id,
        "peptide_id":peptide_id,
        "feature_type":pep.bp_obj.type ,
        "strand":pep.bp_obj.strand,
        "len":len(pep),
        "polyprotein_outline":polyprotein_outline,
        "nb_proteins_annotated_by_pep":len(pep.polyproteins),
        "protein_ids":'|'.join((p.protein_id for p in pep.polyproteins)),
        "nb_protein_annotated_by_pep":len(pep.polyproteins),
        "taxonomy":taxonomy,
        "cover_by_domain": True if pep.fully_included_domains else False,
        "overlapped_by_domain": True if pep.partially_included_domains else False,
        "length_covered_by_fully_included_domains":len_fully_coverage, # length covered by domains that are fully included
        "length_covered_by_all_domains":len_all_coverage, # length covered by domains that are fully included or not
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

def getDomainStat(taxon_id, domain, domain_header):
    domain_dico = get_csv_dico(domain, domain_header)
    try:
        protein_id = domain.protein.protein_id
    except AttributeError:
        protein_id = "Unknown"
        logging.warning(f'Domain annotation {domain.name} has no protein iD in {taxon_id} !!!')
    domain_dico.update({"taxon_id":taxon_id,
                        "protein_id":protein_id})
    return domain_dico

def getCleavageSiteStat(taxon_id, cs, cs_header):
    domain_dico = get_csv_dico(cs, cs_header)

    domain_dico.update({"taxon_id":taxon_id,
                    "protein_id": "|".join((p.protein_id for p in cs.proteins))})
    return domain_dico

def initiateStatFile(taxon, output_dir ):
    taxon = taxon.replace(',', '').replace(' ', '_')


    stat_file_prot = "stat_proteins_{}.csv".format(taxon)
    stat_file_genome = "stat_genomes_{}.csv".format(taxon)
    stat_file_peptides = 'stat_peptides_{}.csv'.format(taxon)
    stat_file_domains =  'stat_domains_{}.csv'.format(taxon)
    stat_file_cleavage_sites =  'stat_cleavage_sites_{}.csv'.format(taxon)

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
                    "len",
                    "nb_proteins_annotated_by_pep",
                    "protein_ids",
                    "polyprotein_outline",
                    "nb_protein_annotated_by_pep",
                    "cover_by_domain",
                    "overlapped_by_domain", # are domains overlapping the peptide annotation
                    "length_covered_by_fully_included_domains", # length covered by domains that are fully included
                    "length_covered_by_all_domains", # length covered by domains that are fully included or not
                    "First_node",
                    "taxonomy" ]

    csv_writer_peptides = csv.DictWriter(handle_stat_pep, fieldnames=pep_header, delimiter='\t')
    csv_writer_peptides.writeheader()


    handle_stat_genome = open(os.path.join(output_dir, stat_file_genome), "w")
    #taxon
    header = ["taxon_id","nb_protein", "has_annotated_poly","nb_polyprotein",  "has_peptide", "nb_peptides", "taxonomy"]
    handle_stat_genome.write('\t'.join(header)+'\n')



    # INTERPRO DOMAINS annotation stat
    genome_header = ['taxon_id', 'organism', "taxon_of_expectation",
                "expected_number_of_peptide", "peptide", "number_of_final_peptide",
                "expected_number_reached", "polyprotein_fully_covered", "unannotated_parts", "relevant_annotation"] # genome header

    handle_stat_do = open(os.path.join(output_dir, stat_file_domains), "w")

    domain_header = ["taxon_id",
                    "protein_id",
                    "match_id",
                    "method" ,
                    "score" ,
                    "Dbxref" ,
                    "name",
                    "signature_desc" ,
                    "start_in_prot",
                    "end_in_prot",
                    "start",
                    "end",
                    "duplicated",
                    "overlapping",
                    "OverlappingDistance",
                    "left_overlaps_peptide",
                    "right_overlaps_peptide"]

    csv_writer_domains = csv.DictWriter(handle_stat_do, fieldnames=domain_header, delimiter='\t')
    csv_writer_domains.writeheader()

    handle_stat_cs = open(os.path.join(output_dir, stat_file_cleavage_sites), "w")

    cs_header = ["taxon_id",
                "protein_id",
                "start",
                "end",
                "peptide_composition", # is it made from 2 mat_peptide or 1 mat_peptide and an unannotated region?
                "overlapped", # is it overlapped by domain annotations?
                "OverlappingDistance", # minimal distance of overlapping
                "overlap_distance_left",
                "overlap_distance_right"]

    csv_writer_cleavage_sites = csv.DictWriter(handle_stat_cs, fieldnames=cs_header, delimiter='\t')
    csv_writer_cleavage_sites.writeheader()

    files_to_close = [handle_stat_genome, handle_stat_prot, handle_stat_pep, handle_stat_cs, handle_stat_do]

    writer_stat_dict = {"peptide":csv_writer_peptides,
                        "protein":csv_writer_protein,
                        "genome":handle_stat_genome,
                        "domain":csv_writer_domains,
                        "cleavage_site":csv_writer_cleavage_sites}

    return writer_stat_dict, files_to_close

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

    gff_file='data/interpro_results/interproscan-5.30-69.0/domains_viral_sequences.gff3'


    taxonomy_file ="data/taxonomy/taxonomy_virus.txt"
    sp_treshold = 90

    gbff_iter = tax.getAllRefseqFromTaxon(taxon, taxonomy_file)

    taxon = taxon.replace(',', '').replace(' ', '_')

    writer_stat_dict, files_to_close = initiateStatFile(taxon, stat_output_dir )

    handle_prot = False # no protein extraction

    i=None
    for i, gb_dico in enumerate(gbff_iter):
        if (i+1)%1000 == 0:
            print(i+1, 'genomes treated')
        gb_file = gb_dico['gb_file']
        genetic_code = gb_dico['genetic_code']
        taxon_id = gb_dico['taxon_id']
        print(taxon_id)
        ext.extractAndStat(gb_file, handle_prot, writer_stat_dict, genetic_code, taxon_id, sp_treshold, gff_file)

    if i == None:
        print("No genome have been found with taxon:",taxon)
    else:
        print('Analysis completed for ', i+1, "genomes")

    for f in files_to_close:
        f.close()

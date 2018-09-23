#!/usr/bin/env python3
import taxonomy as tax
import viral_genome_classes as obj
import viruses_statistics as stat
import parser_interpro_results as do
import os, gzip, logging, sys, csv
from Bio import SeqIO
import cluster_of_interest_identification
import subprocess

def problematic_overlapping_events_identification(domain_dict,
                                                    acceptable_overlapping_threshold_percentage,
                                                    acceptable_overlapping_threshold_aa,
                                                    blacklist_threshold_ratio):
    """
    Find domain annotations that seem to always overlap (use of blacklist_threshold).
    In order to ignored them in the rest of the analysis (blacklist)
    All the other domain annotations that are found overlapping more than the acceptable_overlapping_threshold are consider
    as problematic domain annotation and are written in a file

    INPUT:
        domain_dict: dictionnary with as key: domain name and as value a dictionnary
                    containing the folowing content:
                    'overlapping_dist':[],
                    'overlapping_dist_prct':[],
                    'not_overlapping_count':0}

        acceptable_overlap_threshold : threshold of the overlapping length domain
                                           (in aa or in percentage) which is consider
                                           as accepatble
        blacklist_threshold : proportion or number of time
                              a domain annotation has
                              to be found overlapping to be
                              ignored by the program
    """

    blacklisted_domains= {}
    problematic_domains = {}
    safe_domains = {}


    for domain, info in domain_dict.items():
        # SET UP nb of acceptable overlapping event
        nb_acceptable_overlap = 0
        nb_overlap = 0
        for dist, prct in zip(info['overlapping_dist'], info['overlapping_dist_prct']):
            if int(dist) <  acceptable_overlapping_threshold_aa or float(prct) < acceptable_overlapping_threshold_percentage:
                nb_acceptable_overlap += 1
            else:
                nb_overlap += 1

        overlapping_time_ratio = nb_overlap/info['total']
        if overlapping_time_ratio >= blacklist_threshold_ratio:
            print('//'*20)
            print('BLACKLISTED DOMAIN')
            print('acceptable overlapping',nb_acceptable_overlap )
            print('OVERLAPPING',nb_overlap )
            print(domain)
            [print(k,v) for k, v in info.items()]
            blacklisted_domains[domain] = info

        elif  nb_overlap != 0:
            print('*!!*'*20)
            print("overlapping_time_ratio", overlapping_time_ratio)
            print('PROBLEMATIC DOMAIN')
            print('acceptable overlapping',nb_acceptable_overlap )
            print('OVERLAPPING',nb_overlap )
            print(domain)
            [print(k,v) for k, v in info.items()]
            problematic_domains[domain]=info
        else:
            safe_domains[domain]=info


    print('blacklisted',len(blacklisted_domains))
    print('safe_domain', len(safe_domains))
    print('problematic domain', len(problematic_domains))


def get_domain_overlapping_dict(domain_stat_file):
    """
    From the stat file of domain, compute a dictionary with overlapping info for each domain annotation.
    """
    domain_dico = {}
    with open(domain_stat_file) as fl:

        reader = csv.DictReader(fl, delimiter='\t')

        for row in reader:
            if row['name'] not in domain_dico:
                domain_dico[row['name']] = {'overlapping_dist':[],
                                            'overlapping_dist_prct':[],
                                            'not_overlapping_count':0,
                                            "total":0}

            if row['overlaps_cleavage_site'] == 'True':
                domain_dico[row['name']]["overlapping_dist"].append(row['overlapping_distance'])
                domain_dico[row['name']]["overlapping_dist_prct"].append(row['overlapping_size_percentage'])

            else:
                domain_dico[row['name']]["not_overlapping_count"] += 1
            domain_dico[row['name']]['total'] += 1
    return domain_dico


def getPolyproteinDomainsStat(gb_file, writer_stat_dict, genetic_code, taxon_id, sp_treshold, gff_file):
    nb_peptide = 0
    nb_cds = 0

    genome = obj.gb_file_parser(gb_file, taxon_id, sp_treshold)
    if any((True for s in genome.segments if s.peptides)):
        do.getMatchObject(genome, gff_file)
        do.associateDomainsWithPolyprotein(genome)

    for i, segment in enumerate(genome.segments):
        taxonomy = segment.record.annotations['taxonomy']
        nb_cds += len(segment.cds)

        # WRITE PEPTIDE STAT
        if writer_stat_dict:
            for pep in segment.peptides:
                pep_info_dict = getPepStat(pep, taxon_id, taxonomy)
                writer_stat_dict['peptide'].writerow(pep_info_dict)


    if writer_stat_dict:
        # WRITE DOMAIN STAT
        domain_header = writer_stat_dict['domain'].fieldnames
        for domain in genome.matchs:
            if domain.protein.polyprotein: # we are interested only to get overlapping data on CDS that have mat_peptide annotations
                dict_to_write = getDomainStat(taxon_id, domain, domain_header)
                writer_stat_dict['domain'].writerow(dict_to_write)

        # WRITE CLEAVAGE SITE STAT
        cs_header = writer_stat_dict['cleavage_site'].fieldnames
        for cs in segment.cleavage_sites:
            # print(cs)
            info_dict = getCleavageSiteStat(taxon_id, cs, cs_header)
            writer_stat_dict['cleavage_site'].writerow(info_dict)


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



def getPepStat(pep, taxon_id, taxonomy):
    #Check if the cds has a signal peptide scheme
    #Meaning that it has from 1 to 2 pep annotations and one of it annotation start at the 5' of the cds

    # start/end position inverse if strand -1
    polyprotein_outline = any((True for p in pep.polyproteins if p.polyprotein))

    # fully_included_positions = getNonOverlappingCoveragePositionsPairs(pep.fully_included_domains.values())
    # len_fully_coverage = sum((end - start+1 for start, end in fully_included_positions))

    # included_positions = getNonOverlappingCoveragePositionsPairs(list(pep.partially_included_domains.values())+list(pep.fully_included_domains.values()))
    # len_all_coverage = sum((end - start+1 for start, end in included_positions))


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
        "cover_by_domain": True if pep.included_domains['fully'] else False,
        "overlapped_by_domain": True if pep.included_domains['partially'] else False,
        # "length_covered_by_fully_included_domains":len_fully_coverage, # length covered by domains that are fully included
        # "length_covered_by_all_domains":len_all_coverage, # length covered by domains that are fully included or not
        }

    return dict_info


def getDomainStat(taxon_id, domain, domain_header):
    """
    Information about the domain annotation and the protein
    in case the domain is not overlapping, return the line with general info
    otherwise when the annotation overlaps it return a line giving info about the
    sum of the biggest overlaping events on each side (left/right)
    An overlapping event is a cleavage site overlapped by the domain
    the distance of overlapping correspond to the minimal distance of the domain
    on the right or on the left of the cleavage site
    """
    domain_dico = get_csv_dico(domain, domain_header)

    protein_id = domain.protein.protein_id

    domain_dico.update({"taxon_id":taxon_id,
                        "protein_id":protein_id,

                        "overlaps_cleavage_site": True if domain.overlapped_cleavage_sites else False , #overlaps a cleavage site or not?
                        "nb_of_cleavage_sites_overlapped": len(domain.overlapped_cleavage_sites),

                        "protein_fully_covered":True if len(domain.protein.unannotated_region) == 0 else False, # is the CDS fully covered by mat_peptide annotations or not?
                        "nb_peptide_in_protein":len(domain.protein.peptides),
                        'nb_unannotated_part_in_protein':len(domain.protein.unannotated_region),
                        "polyprotein_outline":domain.protein.polyprotein,
                        "non_polyprotein_explanation":domain.protein.non_polyprotein_explanation
                        })

    # just one line when domain does not overlap a cs
    if not domain.overlapped_cleavage_sites:
        return domain_dico

    # write line containing info about the biggest overlapping event
    lines= []
    [print(k,v) for k, v  in domain.overlapped_cleavage_sites.items()]

    distance_max = {'right':0,"left":0 }
    for overlapped_cs in domain.overlapped_cleavage_sites.values():
        # cs_obj: {'overlapping_side': 'left', 'overlapping_distance': 2}
        if overlapped_cs['overlapping_distance'] > distance_max[overlapped_cs['overlapping_side']]:
             distance_max[overlapped_cs['overlapping_side']] = overlapped_cs['overlapping_distance']

    domain_dico['overlapping_distance'] = sum(distance_max.values())
    domain_dico['overlapping_size_percentage'] = domain_dico['overlapping_distance']/len(domain) * 100
    print()
    print('final line')
    [print(k,v) for k, v  in domain_dico.items()]
    # input()
    print('///'*20)



    return domain_dico


def getCleavageSiteStat(taxon_id, cs, cs_header):
    domain_dico = get_csv_dico(cs, cs_header)



    domain_dico.update({"taxon_id":taxon_id,
                    "protein_id": "|".join((p.protein_id for p in cs.proteins)),
                    "overlapped": True if cs.overlapping_domains else False})
    return domain_dico


def getNonOverlappingCoveragePositionsPairs(positions):
    """
    listposition is a list of tuple (start, end) that potentially  overlaped each other
    We want to find pair of positions that represent the coverage...
    """
    positions = sorted(positions, key=lambda x: x[0], reverse=False) # sort positions by start

    try:
        s, e = positions[0]
        positions.remove((s, e))
    except IndexError:
        return []

    non_overlapping_pairs = []
    for new_s, new_e in positions: # positions is sorted by starts
        if s <= new_s <= e: # the two pair are overlapping
            e = max(new_e, e)
        else: # the two pait are not overlapping then s and e are saved and replace by new s et new e
            non_overlapping_pairs.append((s, e))
            s = new_s
            e = new_e

    non_overlapping_pairs.append((s, e))
    return non_overlapping_pairs




def initiateDomainStatFile(taxon, output_dir ):
    taxon = taxon.replace(',', '').replace(' ', '_')

    stat_file_peptides = 'domain_peptides_stat_{}.csv'.format(taxon)
    stat_file_domains =  'domain_stat_{}.csv'.format(taxon)
    stat_file_cleavage_sites =  'domain_cleavage_sites_stat_{}.csv'.format(taxon)

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

    #
    # genome_header = ['taxon_id', 'organism', "taxon_of_expectation",
    #             "expected_number_of_peptide", "peptide", "number_of_final_peptide",
    #             "expected_number_reached", "polyprotein_fully_covered", "unannotated_parts", "relevant_annotation"] # genome header

    # INTERPRO DOMAINS annotation stat
    handle_stat_do = open(os.path.join(output_dir, stat_file_domains), "w")

    domain_header = ["taxon_id",
                    "protein_id",
                    "polyprotein_outline",
                    "non_polyprotein_explanation",
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
                    "overlaps_cleavage_site",
                    'nb_of_cleavage_sites_overlapped',
                    "duplicated",
                    "overlaps_cleavage_site",
                    "OverlappingDistance",
                    "left_overlaps_peptide",
                    "right_overlaps_peptide",
                    "protein_fully_covered",
                    "nb_peptide_in_protein",
                    'nb_unannotated_part_in_protein',
                    'overlapping_distance',
                    'overlapping_side',
                    'overlapping_size_percentage']

    csv_writer_domains = csv.DictWriter(handle_stat_do, fieldnames=domain_header, delimiter='\t')
    csv_writer_domains.writeheader()

    handle_stat_cs = open(os.path.join(output_dir, stat_file_cleavage_sites), "w")

    cs_header = ["taxon_id",
                "protein_id",
                "start",
                "end",
                "peptide_composition", # is it made from 2 mat_peptide or 1 mat_peptide and an unannotated region?
                "overlapped" # is it overlapped by domain annotations?
                ]

    csv_writer_cleavage_sites = csv.DictWriter(handle_stat_cs, fieldnames=cs_header, delimiter='\t')
    csv_writer_cleavage_sites.writeheader()

    files_to_close = [handle_stat_pep, handle_stat_cs, handle_stat_do]

    writer_stat_dict = {"peptide":csv_writer_peptides,
                        "domain":csv_writer_domains,
                        "cleavage_site":csv_writer_cleavage_sites}

    return writer_stat_dict, files_to_close



if __name__ == '__main__':

    from time import clock; START_TIME = clock()

    logging.basicConfig(filename='log/domain_stats.log',level=logging.INFO)
    logging.basicConfig(level=logging.INFO)
    # $taxon $stat_output_dir $taxonomy_file $sp_treshold $gff_file $stat_protein_file

    taxon = sys.argv[1]

    stat_output_dir = sys.argv[2] # if not given then we don't compute any stat

    print(f"Write domains stat on {taxon} in {stat_output_dir}")

    taxonomy_file = sys.argv[3]

    sp_treshold = int(sys.argv[4]) #90

    gff_file= sys.argv[5] #'data/interpro_results/interproscan-5.30-69.0/domains_viral_sequences.gff3'
    stat_protein_file =  sys.argv[6] #"results/stat_viral_protein/stat_proteins_Viruses.csv"

    annotated_taxons = cluster_of_interest_identification.getAnnotatedProteinTaxId(stat_protein_file, "polyprotein_outline", "True")

    print('annotated_taxons', len(annotated_taxons))



    gbff_iter = tax.getAllRefseqFromTaxonIdList(annotated_taxons, taxonomy_file)
    # print(len(list(gbff_iter)))
    # input()

    taxon = taxon.replace(',', '').replace(' ', '_')

    if stat_output_dir:
        writer_stat_dict, files_to_close = initiateDomainStatFile(taxon, stat_output_dir)
    else:
        files_to_close = []
        writer_stat_dict = False


    i=None
    for i, gb_dico in enumerate(gbff_iter):
        # print(gb_dico)
        if (i+1)%100 == 0:
            print(i+1, 'genomes analysed')
        # print(gb_dico)
        gb_file = gb_dico['gb_file']
        genetic_code = gb_dico['genetic_code']
        taxon_id = gb_dico['taxon_id']

        getPolyproteinDomainsStat(gb_file, writer_stat_dict, genetic_code, taxon_id, sp_treshold, gff_file)

    if i == None:
        raise NameError("No genome have been found with taxon:",taxon)
    else:
        print(f'Analysis completed for {i+1} genomes')

    for f in files_to_close:
        f.close()

    print('\nPROGRAM RUN TIME:%6.2f'%(clock()-START_TIME), 'seconds.');
    print('-'*20)

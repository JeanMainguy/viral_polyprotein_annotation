#!/usr/bin/env python3
import taxonomy as tax
import viral_genome_classes as obj
import viruses_statistics as stat
import parser_interpro_results as do
import os, gzip, logging, sys, csv
from Bio import SeqIO
import cluster_of_interest_identification


def getPolyproteinDomainsStat(gb_file, writer_stat_dict, genetic_code, taxon_id, sp_treshold, gff_file):
    nb_peptide = 0
    nb_cds = 0

    genome = obj.gb_file_parser(gb_file, taxon_id, sp_treshold)
    if any((True for s in genome.segments if s.peptides)):
        do.getMatchObject(genome, gff_file)
        do.associateMatchWithPolyprotein(genome)

    for i, segment in enumerate(genome.segments):
        taxonomy = segment.record.annotations['taxonomy']
        nb_cds += len(segment.cds)
        nb_polyprotein = 0

        for p, cds in enumerate(segment.cds):
            if cds.polyprotein:
                nb_polyprotein += 1

            # #WRITE PROTEIN STAT
            # if  writer_stat_dict: # may be False if we don't want to write stat
            #     info_dict = stat.getProteinStat(cds, taxon_id, taxonomy)
            #     writer_stat_dict['protein'].writerow(info_dict)
            #
            #     nb_peptide += len(cds.peptides)


        # WRITE PEPTIDE STAT
        if writer_stat_dict:
            for pep in segment.peptides:
                pep_info_dict = stat.getPepStat(pep, taxon_id, taxonomy)
                writer_stat_dict['peptide'].writerow(pep_info_dict)


    if writer_stat_dict:
        # WRITE GENOME STAT
        # stat.writeGenomeStat(taxon_id, nb_cds, nb_peptide, writer_stat_dict['genome'], taxonomy, nb_polyprotein)

        # WRITE DOMAIN STAT
        domain_header = writer_stat_dict['domain'].fieldnames
        for domain in genome.matchs:
            info_dict = stat.getDomainStat(taxon_id, domain, domain_header)
            writer_stat_dict['domain'].writerow(info_dict)

        # WRITE CLEAVAGE SITE STAT
        cs_header = writer_stat_dict['cleavage_site'].fieldnames
        for cs in segment.cleavage_sites:
            # print(cs)
            info_dict = stat.getCleavageSiteStat(taxon_id, cs, cs_header)
            writer_stat_dict['cleavage_site'].writerow(info_dict)

def initiateDomainStatFile(taxon, output_dir ):
    taxon = taxon.replace(',', '').replace(' ', '_')

    stat_file_peptides = 'stat_peptides_domains_{}.csv'.format(taxon)
    stat_file_domains =  'stat_domains_{}.csv'.format(taxon)
    stat_file_cleavage_sites =  'stat_cleavage_sites_domains{}.csv'.format(taxon)

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
                    "fully included",
                    'partially included',
                    "duplicated",
                    "overlapping",
                    "OverlappingDistance",
                    "left_overlaps_peptide",
                    "right_overlaps_peptide",
                    "protein_fully_covered",
                    "nb_peptide_in_protein",
                    'nb_unannotated_part_in_protein']

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

    files_to_close = [handle_stat_pep, handle_stat_cs, handle_stat_do]

    writer_stat_dict = {"peptide":csv_writer_peptides,
                        "domain":csv_writer_domains,
                        "cleavage_site":csv_writer_cleavage_sites}

    return writer_stat_dict, files_to_close



if __name__ == '__main__':
    from time import clock; START_TIME = clock()

    logging.basicConfig(filename='log/domain_stats.log',level=logging.INFO)

    try:
        taxon = sys.argv[1] # if not given then we don't compute any stat
    except:
        taxon = "Viruses"
    try:
        stat_output_dir = sys.argv[2]
    except:
        stat_output_dir = 'test'

    print(f"Write domains stat on {taxon} in {stat_output_dir}")
    taxonomy_file ="data/taxonomy/taxonomy_virus.txt"
    sp_treshold = 90

    gff_file='data/interpro_results/interproscan-5.30-69.0/domains_viral_sequences.gff3'
    stat_protein_file = "results/stat_viral_protein/stat_proteins_Viruses.csv"

    annotated_taxons = cluster_of_interest_identification.getAnnotatedProteinTaxId(stat_protein_file, "has_peptide", "True")

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

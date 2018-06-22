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
from numpy import std

def getCdsObject(taxon_id, protein_id, taxonomy_file):
    #TODO this step should be optimize !
    try:
        gb_file_dict = next(tax.getAllRefseqFromTaxon(taxon_id, taxonomy_file))
    except StopIteration:
        return None
    gb_file = gb_file_dict['gb_file']

    genome = obj.Genome( gb_file)
    with gzip.open(gb_file, "rt") as handle:
        for i, record in enumerate(SeqIO.parse(handle, "genbank")):


            segment = obj.Segment(record, gb_file)
            genome.segments.append(segment)

            segment.getMatpeptidesAndPolyproteins()
            cds_of_interest = None
            for cds in segment.cds:
                if cds.protein_id == protein_id:
                    cds_of_interest = cds
                    break
            if not cds_of_interest:
                continue

            segment.checkPeptideRedundancy() #remove the redundant peptide
            segment.checkSubPeptides()
            segment.associatePepWithProt()

            segment.identifySubProtein()

            segment.getCleavageSites()
            break #we have found the protein so no need to continue the parsing

    if not cds.cleavage_sites:
        logging.warning('protein has no cleavage site ... :-§')
    return cds

def convertAlignmentToList(alignment_seq):
    """
    Take fasta sequence of the alignment so the sequence with gap symbol
    and make a list with indice of the ungap sequence and None when it is a gap
    the following sequence :
    --ABC-D--E
    would become :
    [None, None, 0, 1, 2, 3, None,4, None, None, 5]
    """
    gap="-"
    list=[]
    indice_seq = 0
    for i, e in enumerate(alignment_seq):
        if e != gap:
            indice_seq += 1
            list.append(indice_seq)

        else:
            list.append(None)
    return list

def findCloseSites(site, site_start_in_aln, other_cds, window):

    site.neighboring_sites = set()

    for cds in other_cds:
        seq_part = cds.aligned_sequence[site_start_in_aln-window-1:site_start_in_aln+window+1]
        aln_part_list = cds.aln_list[site_start_in_aln-window-1:site_start_in_aln+window+1]
        # print(aln_part_list)
        # print(seq_part)

        for position in aln_part_list:

            if position in cds.cleavage_site_positions:

                adj_site = cds.cleavage_site_positions[position]
                site.neighboring_sites.add(adj_site)
                # print("correpsond to a CS", adj_site.taxon_id, adj_site.start_aa(cds) )
                # print(seq_part)
                # print(aln_part_list)


def checkCleavageSiteAcrossAlignment(cds_list, window):
    all_sites = set()
    for i, cds in enumerate(cds_list):
        all_sites |= set(cds.cleavage_sites)
        # print()
        # print(cds.number, cds.cleavage_site_positions)
        for site in cds.cleavage_sites:
            site.cds_of_aln = cds
            site_start = site.start_aa(cds) # to be in base 0 we need the -1
            site_start_in_aln = cds.aln_list.index(site_start)
            # print('CLEAVAGE SITE:::: from', site.taxon_id, site.start_aa(cds))
            # print(cds.aligned_sequence[site_start-window-1:site_start+window+1])
            findCloseSites(site, site_start_in_aln,  cds_list[:i]+cds_list[i+1:], window)

    #make group of cleavage sites
    grouped_sites = set()
    group_list = []
    for site in all_sites:
        if site not in grouped_sites:
            group = set()
            groupCleavageSites(group, site)
            grouped_sites |= group #add site that are grouped in a the set to not process them again
            group_list.append(group)
    print('NB_GROUP', len(group_list))
    return group_list

def groupCleavageSites(group, site):
    #recursivity
    group.add(site)
    for adj_site in site.neighboring_sites:
        if adj_site not in group:
            groupCleavageSites(group, adj_site)

def getStatOnGroup(group):
    mean_positions = []
    protein_ids = set()
    for site in group:
        cds = site.cds_of_aln
        site_start = site.start_aa(cds)
        site_end =site.end_aa(cds)
        protein_ids.add(cds.protein_id)

        mean_position_in_aln = sum([cds.aln_list.index(site_start), cds.aln_list.index(site_end)])/2.0
        mean_positions.append(mean_position_in_aln)

    nb_protein_in_group = len(protein_ids)
    nb_of_cleavage_site = len(group)
    nb_cleavage_from_the_same_protein =  nb_of_cleavage_site - nb_protein_in_group
    standard_dev = std(mean_positions)

    return [ nb_of_cleavage_site, nb_protein_in_group, nb_cleavage_from_the_same_protein, standard_dev, round(standard_dev)]

def processAlignmentFile(alignment_file, taxonomy_file, window, output_file, cluster_nb):
    cds_list = []
    cds_annotated = 0
    with open(alignment_file, "rU") as handle:
        for record in SeqIO.parse(handle, "clustal") :

            # 11764|YP_009109691.1|564.0|Peptide
            taxon_id, protein_id, length, hasPeptide =  record.id.split('|')

            cds = getCdsObject(taxon_id, protein_id, taxonomy_file)
            if cds is None:
                continue
            cds.cleavage_site_positions = {s.start_aa(cds):s for s in cds.cleavage_sites}
            cds.aligned_sequence = record.seq
            cds.aln_list = convertAlignmentToList(cds.aligned_sequence)

            cds_list.append(cds)
            if cds.cleavage_sites:
                cds_annotated += 1
            # print(cds.aligned_sequence)
            # print(cds.aln_list)

    groups_of_cs = checkCleavageSiteAcrossAlignment(cds_list, window)
    header = ["cluster_nb",'nb_seq', 'nb_seq_with_annotation', "nb_of_cleavage_site", 'nb_protein_in_group', 'nb_cleavage_from_the_same_protein', "standard_deviation", 'round_std', "group_completness", "completness_category"]
    general_info = [cluster_nb, len(cds_list), cds_annotated ]
    with open(output_file, 'w') as writer:
        writer.write("\t".join(header)+"\n")
        for group in groups_of_cs:
            stat_info = getStatOnGroup(group)
            nb_protein_in_group = stat_info[1]
            completeness = (nb_protein_in_group / cds_annotated )*100
            if completeness < 25:
                category = 'from 0 to 25%'
            elif 25 <= completeness < 50:
                category = 'from 25% to 50%'
            elif 50 <= completeness < 75:
                category = 'from 50% to 75%'
            elif 75 <= completeness:
                category = 'from 75% to 100%'
            stat_info.append(completeness)
            stat_info.append(category)
            list_info_group = general_info + stat_info
            list_info_group = [str(e) for e in list_info_group]
            writer.write("\t".join(list_info_group)+"\n")


if __name__ == '__main__':

    logging.basicConfig(filename='log/alignment_analysis.log', level=logging.WARNING)
    # alignment_file= "data/alignment/Retro-transcribing_viruses_1e-5_coverage90_I3/seq_cluster16.aln"
    alignment_file=  sys.argv[1] #'data/alignment/Viruses_1e-5_coverage90_I2/seq_cluster1037.aln'
    taxonomy_file =  sys.argv[3]#"data/taxonomy/taxonomy_virus.txt"
    output_file =  sys.argv[2]#'data/alignment/Viruses_1e-5_coverage90_I2/seq_cluster1037_stat.csv'
    window =  int(sys.argv[4])#20

    re_result = re.search("cluster(\d+).aln", alignment_file)
    cluster_nb = "NA" if not re_result else re_result.group(1)
    processAlignmentFile(alignment_file, taxonomy_file, window, output_file, cluster_nb)

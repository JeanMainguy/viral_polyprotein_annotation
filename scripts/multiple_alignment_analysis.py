#!/usr/bin/env python3

import taxonomy as tax
import viral_genome_classes as obj
import viruses_statistics as stat
import parser_interpro_results as do
import visualisation_alignment as view_aln

import os, gzip, logging
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import sys, csv, re
from operator import attrgetter
from numpy import std, mean
import operator

def getCdsObject(taxon_id, protein_id, taxonomy_file, gff_file, sp_treshold):
    #TODO this step should be optimize !
    try:
        gb_file_dict = next(tax.getAllRefseqFromTaxon(taxon_id, taxonomy_file))
    except StopIteration:
        return None
    gb_file = gb_file_dict['gb_file']

    genome = obj.gb_file_parser(gb_file, taxon_id, sp_treshold)
    do.getMatchObject(genome, gff_file)
    do.associateMatchWithPolyprotein(genome)
    for segment in genome.segments:
        do.getDomainOverlappingInfo(segment)

    for segment in genome.segments:
        for cds in segment.cds:
            if cds.protein_id == protein_id:
                cds_of_interest = cds
                if not cds.cleavage_sites:
                    logging.warning('protein has no cleavage site ... :-§')
                return cds, segment

def convertAlignmentToList(alignment_seq):
    """
    Take fasta sequence of the alignment so the sequence with gap symbol
    and make a list with indice of the ungap sequence and None when it is a gap
    the following sequence :
    --ABC-D--E
    would become :
    [None, None, 0, 1, 2, None,3, None, None, 4]
    """
    gap="-"
    list=[]
    indice_seq = 0
    for i, e in enumerate(alignment_seq):
        if e != gap:

            list.append(indice_seq)
            indice_seq += 1
        else:
            list.append(None)
    return list

def findCloseSites(site, site_start_in_aln, cds_list, window):

    site.neighboring_sites = set()

    for cds in cds_list:
        seq_part = cds.aligned_sequence[site_start_in_aln-int(window/2)+1:site_start_in_aln+int(window/2)+1]

        aln_part_list = cds.aln_list[site_start_in_aln-int(window/2)+1:site_start_in_aln+int(window/2)+1]
        # print(aln_part_list)
        # print(seq_part)

        for position in aln_part_list:

            if position in cds.cleavage_site_positions:

                adj_site = cds.cleavage_site_positions[position]
                if adj_site != site:
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
            site_start = site.start_aa(cds) -1 # to be in base 0 we need the -1

            site_start_in_aln = cds.aln_list.index(site_start)
            site.start_in_aln = site_start_in_aln
            # print(cds.aln_list[site_start_in_aln-window-1:site_start_in_aln+window+1])
            findCloseSites(site, site_start_in_aln,  cds_list, window)

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
    # sort list of group by the position in aln of the first site of each group
    return sorted(group_list, key=lambda x: list(x)[0].start_in_aln)


def groupCleavageSites(group, site):
    #recursivity
    group.add(site)
    for adj_site in site.neighboring_sites:
        if adj_site not in group:
            groupCleavageSites(group, adj_site)


def getStatOnGroup(group, cds_annotated):
    mean_positions = []
    protein_ids = set()

    for site in group:
        cds = site.cds_of_aln
        site_start = site.start_aa(cds)
        site_end =site.end_aa(cds)
        protein_ids.add(cds.protein_id)

        mean_position_in_aln = sum([cds.aln_list.index(site_start), cds.aln_list.index(site_end)])/2.0
        mean_positions.append(mean_position_in_aln)
    standard_dev = std(mean_positions)
    nb_of_cleavage_site = len(group)
    nb_protein_in_group = len(protein_ids)

    completeness = (nb_protein_in_group / cds_annotated )*100

    if completeness < 25:
        category = 'from 0 to 25%'
    elif 25 <= completeness < 50:
        category = 'from 25% to 50%'
    elif 50 <= completeness < 75:
        category = 'from 50% to 75%'
    elif 75 <= completeness:
        category = 'from 75% to 100%'

    dico_info = {
    'nb_protein_in_group' : nb_protein_in_group,
    'nb_of_cleavage_site' : nb_of_cleavage_site,
    'nb_cleavage_from_the_same_protein' :  nb_of_cleavage_site - nb_protein_in_group,
    'standard_deviation' : standard_dev,
    "round_std" : round(standard_dev),
    "group_completness" : completeness,
    "completness_category": category,
    "average_position_in_aln":mean(mean_positions)
    }
    #[ nb_of_cleavage_site, nb_protein_in_group, nb_cleavage_from_the_same_protein, standard_dev, round(standard_dev)]
    return dico_info


def getPositionInSeq(cds, position, window):
    position_in_seq = []
    print([a for a in str(cds.aligned_sequence[position-3:position+3])])
    print(str(cds.aligned_sequence[position]))
    print(cds.aln_list[position-3:position+3])
    print(cds.aln_list[position])
    if cds.aln_list[position]: # when the position corresponf to a gap then it is None
        position_in_seq.append(cds.aln_list[position])

    for i in range(1, int(window/2)):

        if cds.aln_list[position + i]:
            position_in_seq.append(cds.aln_list[position + i])
            if len(position_in_seq) == 2:
                break

        if  cds.aln_list[position - i]:

            position_in_seq.append(cds.aln_list[position - i])
            if len(position_in_seq) == 2:
                break

    position_in_seq.sort()

    if len(position_in_seq) < 2:
        print("No aa in the window for the cds")
        print('could no propagate the cleavage site')
        print(cds.aligned_sequence[position-int(window/2)+1:position+int(window/2)+1])
    print(position_in_seq)
    if 399 in position_in_seq:
        input()
    return position_in_seq


def propagate_cleavage_sites(group, general_info, cds_list):
    window = general_info['window']
    average_po = int(general_info['average_position_in_aln']) -1
    print(average_po)
    iter_cds_non_annotated = (cds for cds in cds_list if not cds.polyprotein)

    for blank_cds in iter_cds_non_annotated:
        # print(blank_cds)
        # print(blank_cds.aln_list[average_po-window:average_po+window])
        position_cs = getPositionInSeq(blank_cds, average_po, window)
        start, end = position_cs[0]+1, position_cs[1]+1
        new_site = obj.PredictedCleavageSite(start, end, blank_cds, group, 'Unknown')


def processAlignmentFile(alignment_file, taxonomy_file, windows, csv_writer, cluster_nb, sp_treshold, gff_file):
    cds_list = []
    cds_annotated = 0
    with open(alignment_file, "rU") as handle:
        for record in SeqIO.parse(handle, "clustal") :
            print(record)

            # 11764|YP_009109691.1|564
            taxon_id, protein_id, length =  record.id.split('|')
            cds, segment = getCdsObject(taxon_id, protein_id, taxonomy_file, gff_file, sp_treshold )
            if cds is None:
                logging.warning(f"TAXONID:{taxon_id} has not been found in taxonomy file")
                continue
            cds.cleavage_site_positions = {s.start_aa(cds):s for s in cds.cleavage_sites}
            cds.aligned_sequence = record.seq
            cds.aln_list = convertAlignmentToList(cds.aligned_sequence)

            cds_list.append(cds)
            if cds.cleavage_sites:
                cds_annotated += 1
            # print(cds.aligned_sequence)
            # print(cds.aln_list)
    groups_of_cs_wind = {}
    for window in windows:
        groups_of_cs_wind[window] = checkCleavageSiteAcrossAlignment(cds_list, window)

    general_info = {"cluster_nb":cluster_nb,
                    'nb_seq':len(cds_list),
                    'nb_seq_with_annotation':cds_annotated}

    for window, groups_of_cs in groups_of_cs_wind.items():
        general_info["window"] = window
        for i, group in enumerate(groups_of_cs):
            #print(group)
            group_info = getStatOnGroup(group, cds_annotated)

            group_info.update(general_info)

            propagate_cleavage_sites(group, group_info, cds_list)

            print(group_info)

            print('='*20)
            print('Group',i)
            {print(k, v) for k, v in group_info.items()}

            if csv_writer:
                csv_writer.writerow(group_info)
        #    writer.write("\t".join(list_info_group)+"\n")

    view_aln.visualisation(cds_list, alignment_file, 150)

def initiate_ouput(output_file):

    header = ["cluster_nb",
    'nb_seq',
    'nb_seq_with_annotation',
    "nb_of_cleavage_site",
    'nb_protein_in_group',
    'nb_cleavage_from_the_same_protein',
    "standard_deviation",
    'round_std',
    "average_position_in_aln",
    "group_completness",
    "completness_category",
    "window"]

    handle_cs_out = open(output_file, "w")
    csv_writer = csv.DictWriter(handle_cs_out, fieldnames=header, delimiter='\t')
    csv_writer.writeheader()
    return csv_writer, handle_cs_out



if __name__ == '__main__':

    logging.basicConfig(filename='log/alignment_analysis.log', level=logging.INFO)
    # alignment_file= "data/alignment/Retro-transcribing_viruses_1e-5_coverage90_I3/seq_cluster16.aln"

    try:
        alignment_file=  sys.argv[1] #'data/alignment/Viruses_1e-5_coverage90_I2/seq_cluster1037.aln'
    except:
        alignment_file=  'data/alignment/RefSeq_download_date_2018-07-21/Viruses_evalue_1e-60coverage20_I1_4/seq_cluster3.aln'

    try:
        taxonomy_file =  sys.argv[3]#"data/taxonomy/taxonomy_virus.txt"
        output_file =  sys.argv[2]#'data/alignment/Viruses_1e-5_coverage90_I2/seq_cluster1037_stat.csv'
        windows_input = sys.argv[4]

    except:
        taxonomy_file =  "data/taxonomy/taxonomy_virus.txt"
        output_file =  'test/aln_analysis.csv'
        windows_input = "40"

    gff_file = 'data/interpro_results/interproscan-5.30-69.0/domains_viral_sequences.gff3'

    # window includ the cleavage in the middle
    # cleavage sites have a length of 2
    # then a window will always be even and >= 2
    # to take of that we add 1 to uneven window
    windows =  {int((int(w)+int(w)%2) ) for w in windows_input.split(' ')} #10,20,30
    assert min(windows) >= 0

    print('window used to analyse cleavage sites',  windows)

    sp_treshold = 90

    re_result = re.search("cluster(\d+).aln", alignment_file)
    cluster_nb = "NA" if not re_result else re_result.group(1)

    print("PROCESS of CLUSTER ", cluster_nb)
    csv_writer, handle_cs_out =  initiate_ouput(output_file)
    processAlignmentFile(alignment_file, taxonomy_file, windows, csv_writer, cluster_nb, sp_treshold, gff_file)
    handle_cs_out.close()

#!/usr/bin/env python3

import taxonomy as tax
import viral_genome_classes as obj
import viruses_statistics as stat
import parser_interpro_results as do
import visualisation_alignment as view_aln
import visualisation_protein as view_prot
import viral_protein_extraction as extract

from os import path, listdir
import gzip
import logging
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import sys
import csv
import re
from operator import attrgetter
from numpy import std, mean
import operator


def getCdsObject(taxon_prot_ids, taxonomy_file, interpro_domains_file, sp_treshold):

    gb_dict_iter = tax.getAllRefseqFromTaxonIdList(taxon_prot_ids, taxonomy_file)

    cds_list = []
    for gb_file_dict in gb_dict_iter:
        cds_found = False
        taxon_id = gb_file_dict['taxon_id']
        gb_file = gb_file_dict['gb_file']

        genome = obj.gb_file_parser(gb_file, taxon_id, sp_treshold)
        if interpro_domains_file:  # interpro_domains_file can be set as None and the domain annotation part is not computed
            do.getMatchObject(genome, interpro_domains_file)
            do.associateDomainsWithPolyprotein(genome)
            # for segment in genome.segments:
            #     do.getDomainOverlappingInfo(segment)

        for segment in genome.segments:

            for cds in segment.cds:
                if cds.protein_id in taxon_prot_ids[taxon_id]:
                    cds_found = True
                    cds.genetic_code = gb_file_dict['genetic_code']
                    cds_list.append(cds)

    if len(cds_list) != sum([len(protein_ids) for protein_ids in taxon_prot_ids.values()]):
        # print(len(cds_list), '!=' ,sum([len(protein_ids) for protein_ids in taxon_prot_ids.values()]))
        set_requested = set()
        for taxon_id, protein_ids in taxon_prot_ids.items():
            set_requested |= {f'{taxon_id}|{protein_id}' for protein_id in protein_ids}

        set_retrieved = set()

        for cds in cds_list:
            set_retrieved.add(f'{cds.segment.taxon_id}|{cds.protein_id}')

        print(set_retrieved)
        print(set_requested - set_retrieved)

    return cds_list


def convertAlignmentToList(alignment_seq):
    """
    Take fasta sequence of the alignment so the sequence with gap symbol
    and make a list with indice of the ungap sequence and None when it is a gap
    the following sequence :
    --ABC-D--E
    would become :
    [None, None, 0, 1, 2, None, 3, None, None, 4]
    """
    gap = "-"
    list = []
    indice_seq = 0
    for i, e in enumerate(alignment_seq):
        if e != gap:

            list.append(indice_seq)
            indice_seq += 1
        else:
            list.append(None)
    return list


def findCloseSites(site, site_start_in_aln, site_end_in_aln,  cds_list, window):

    size_site_in_aln = site_end_in_aln - site_start_in_aln + 1
    # site.neighboring_sites = set()

    for cds in cds_list:
        aln_part_list = cds.aln_list[site_start_in_aln -
                                     int(window/2)+1:site_start_in_aln+(size_site_in_aln-1)+int(window/2)]

        for position in aln_part_list:

            if position in cds.cleavage_site_positions:

                adj_site = cds.cleavage_site_positions[position]
                if adj_site != site:
                    site.neighboring_sites.add(adj_site)
                    adj_site.neighboring_sites.add(site)
                # print("correpsond to a CS", adj_site.taxon_id, adj_site.start_aa(cds) )
                # print(seq_part)
                # print(aln_part_list)


def groupCleavageSitesAcrossAlignment(cds_list, window):
    all_sites = set()
    annotated_cds_list = [cds for cds in cds_list if cds.polyprotein]
    # Reinitialisation of cleavage_site obj
    for cds in cds_list:
        cds.predicted_cleavage_sites = []
        for site in cds.cleavage_sites:
            site.cds_of_aln = []
            site.start_in_aln = None
            site.end_in_aln = None
            site.neighboring_sites = set()

    for i, cds in enumerate(annotated_cds_list):
        all_sites |= set(cds.cleavage_sites)
        # print()
        # print(cds.number, cds.cleavage_site_positions)

        for site in cds.cleavage_sites:
            site.cds_of_aln.append(cds)
            site_start = site.start_aa(cds) - 1  # to be in base 0 we need the -1
            site_end = site.end_aa(cds) - 1
            site_start_in_aln = cds.aln_list.index(site_start)
            site_end_in_aln = cds.aln_list.index(site_end)
            site.start_in_aln = site_start_in_aln
            site.end_in_aln = site_end_in_aln

            findCloseSites(site, site_start_in_aln, site_end_in_aln,  annotated_cds_list, window)

    # make group of cleavage sites
    grouped_sites = set()
    group_list = []
    for site in all_sites:
        if site not in grouped_sites:
            new_group = set()
            groupCleavageSites(new_group, site)

            grouped_sites |= new_group  # add site that are grouped in a the set to not process them again
            group_list.append(new_group)

    # print('NB_GROUP', len(group_list))

    nb_cs_in_group = sum([len(group) for group in group_list])

    # print("len(all_sites) ",len(all_sites) )
    # print('nb_cs_in_group',nb_cs_in_group )
    assert nb_cs_in_group == len(all_sites), 'Cs grouped twice in two different groupe'
    # sort list of group by the position in aln of the first site of each group
    return sorted(group_list, key=lambda x: list(x)[0].start_in_aln)


def groupCleavageSites(group, site):
    # recursivity
    group.add(site)
    for adj_site in site.neighboring_sites - group:
        groupCleavageSites(group, adj_site)


def getStatOnGroup(group, nb_cds_annotated, absent_annotated_cds_penalty):
    mean_positions = []
    positions = []
    protein_in_group_list = []

    for i, site in enumerate(group):
        for cds in site.cds_of_aln:
            """
            The same cleavage site may attached to more than one cds
            It happends when two cds share an identical part due to a ribosomal frameshift/readthrough and
            the cleavage sites in this shared part are identical in the 2 proteins.
            The cleavage sites are not duplicated during the parsing but have more than one cds in
            their attribute 'proteins'.
            But here to compute stat on the group we consider the cleavaged site the number of time it appears in
            the group to make it less confusing
            """
            if not cds.polyprotein:
                print(cds)
                print(site)
                input()
            site_start = site.start_aa(cds)
            site_end = site.end_aa(cds)
            protein_in_group_list += site.cds_of_aln
            positions += [cds.aln_list.index(site_start), cds.aln_list.index(site_end)]
            mean_position_in_aln = sum([cds.aln_list.index(site_start),
                                        cds.aln_list.index(site_end)])/2.0
            mean_positions.append(mean_position_in_aln)

    standard_dev = std(mean_positions)
    average_position_in_aln = mean(mean_positions)
    deviation_length = [abs(average_position_in_aln - po) for po in mean_positions]
    sum_of_deviation = sum(deviation_length)

    nb_protein_in_group = len(set(protein_in_group_list))
    # nb of time a protein appears more than 1 time in the group
    nb_cleavage_from_the_same_protein = len(protein_in_group_list) - nb_protein_in_group

    nb_of_absent_cds = nb_cds_annotated - nb_protein_in_group
    confidence_score = (1 + sum_of_deviation + nb_of_absent_cds *
                        absent_annotated_cds_penalty)/(nb_protein_in_group**2)  # nb_cds_annotated)#**2
    nb_of_cleavage_site = len(group)

    completeness = (nb_protein_in_group / nb_cds_annotated)*100

    # print(f'average_position_in_aln {average_position_in_aln}')
    # print(f'mean position : {mean_positions}')
    # print(f"deviation_length : {deviation_length}")
    # print(f'confidence_score =    (1 + sum_of_deviation + nb_of_absent_cds * absent_annotated_cds_penalty )/(nb_protein_in_group**2 )')
    # print(f'confidence_score =    (1 + {sum_of_deviation} + {nb_of_absent_cds} * {absent_annotated_cds_penalty} )/({nb_protein_in_group**2} ) = {confidence_score}')

    if completeness <= 25:
        category = 'from 0 to 25%'
    elif 25 < completeness <= 50:
        category = 'from 25% to 50%'
    elif 50 < completeness <= 75:
        category = 'from 50% to 75%'
    elif 75 < completeness:
        category = 'from 75% to 100%'

    dico_info = {
        'nb_protein_in_group': nb_protein_in_group,
        'nb_of_cleavage_site': nb_of_cleavage_site,
        'nb_cleavage_from_the_same_protein': nb_cleavage_from_the_same_protein,
        'standard_deviation': standard_dev,
        "round_std": round(standard_dev),
        "group_completeness": completeness,
        "completeness_category": category,
        "average_position_in_aln": average_position_in_aln,
        "position_max": max(positions)-1,
        'position_min': min(positions)-1,
        'confidence_score': confidence_score
    }
    #[ nb_of_cleavage_site, nb_protein_in_group, nb_cleavage_from_the_same_protein, standard_dev, round(standard_dev)]
    return dico_info


def getPositionInSeq(cds, position, window):
    position_in_seq = []
    # print([a for a in str(cds.aligned_sequence[position-3:position+3])])
    # print(str(cds.aligned_sequence[position]))
    # print(cds.aln_list[position-3:position+3])
    # print(cds.aln_list[position])
    if cds.aln_list[position]:  # when the position corresponf to a gap then it is None
        position_in_seq.append(cds.aln_list[position])

    for i in range(1, int(window/2)):
        if len(cds.aln_list) > position+i and cds.aln_list[position + i]:
            position_in_seq.append(cds.aln_list[position + i])
            if len(position_in_seq) == 2:
                break

        if position-i >= 0 and cds.aln_list[position - i]:
            position_in_seq.append(cds.aln_list[position - i])
            if len(position_in_seq) == 2:
                break

    position_in_seq.sort()

    if len(position_in_seq) == 0:
        # print("No aa in the window for the cds")
        # print('could no propagate the cleavage site')
        # print(cds.aligned_sequence[position-int(window/2)+1:position+int(window/2)+1])
        return None
    elif len(position_in_seq) == 1:
        # print("Only one aa in the window for the cds")
        # print('could no propagate the cleavage site completely')
        # print(cds.aligned_sequence[position-int(window/2)+1:position+int(window/2)+1])
        return position_in_seq.append(position_in_seq[0])
    return position_in_seq


def propagate_cleavage_sites(group, general_info, cds_list, window):

    average_po = int(general_info['average_position_in_aln'])
    confidence_score = general_info['confidence_score']
    iter_cds_non_annotated = (cds for cds in cds_list if not cds.polyprotein)

    for blank_cds in iter_cds_non_annotated:
        # print(blank_cds)
        # print(blank_cds.aln_list[average_po-window:average_po+window])
        position_cs = getPositionInSeq(blank_cds, average_po, window)
        if position_cs:
            start, end = position_cs[0], position_cs[1]
            # blank cds is notify during the creation of the obj
            new_site = obj.PredictedCleavageSite(start, end, blank_cds, group, confidence_score)


def filter_cds_cleavage_sites(cds_list):
    # stamp cleavage sites that
    # are made from only one mature peptide border
    for cds in cds_list:
        quality_sum = 0
        for cs in cds.cleavage_sites:

            if sum([1 for p in cs.peptides if p.__class__.__name__ == "Peptide"]) == 2:
                cs.quality = True
                quality_sum += 1
            else:
                cs.quality = False
        if cds.cleavage_sites:
            cds.annotation_quality = quality_sum/len(cds.cleavage_sites)*100


def parse_alignment_file(alignment_file):

    seq_aln_dict = {}
    taxon_prot_ids = {}
    with open(alignment_file, "rU") as handle:
        for record in SeqIO.parse(handle, "clustal"):
            # 11764|YP_009109691.1|564
            taxon_id, protein_id, length = record.id.split('|')
            seq_aln_dict[protein_id] = record.seq
            taxon_prot_ids.setdefault(taxon_id, set()).add(protein_id)

    return taxon_prot_ids, seq_aln_dict


def visualisation_of_processed_aln(cds_list, alignment_file, group_info_list, display_line_size=180):

    cds_annotated = [cds for cds in cds_list if cds.polyprotein]
    # VISUALISATION OF THE ALIGNMENT
    len_cds_max = max((len(cds) for cds in cds_list))
    # for cds in cds_annotated:
    # for cds in cds_list:
    #print(view_prot.visualisation_protein(cds, 1, len_cds_max))

    for cds in cds_list:
        cds.matchs = []
    view_aln.visualisation(cds_list, alignment_file, display_line_size, group_info_list)


def analyse_cleavage_site_groups(cds_list, window):

    groups_of_cs = groupCleavageSitesAcrossAlignment(cds_list, window)

    group_info_list = []

    nb_cds_annotated = sum([1 for cds in cds_list if cds.polyprotein])

    absent_annotated_cds_penalty = window

    for i, group in enumerate(groups_of_cs):
        group_info = getStatOnGroup(group,  nb_cds_annotated, absent_annotated_cds_penalty)
        group_info['group_index'] = i
        # group_info["cleavage_sites"] = group
        group_info_list.append(group_info)

    return group_info_list, groups_of_cs


def compute_alignment_stat(group_info_list, parameters, cds_list, aln_csv_writer, group_site_csv_writer):

    cds_annotated = [cds for cds in cds_list if cds.polyprotein]
    cleavages_sites_per_seq = [len(cds.cleavage_sites) for cds in cds_annotated]

    general_info = {'nb_sequences': len(cds_list),
                    'nb_seq_with_annotation': len(cds_annotated),
                    'cleavages_sites_per_seq': cleavages_sites_per_seq}
    general_info.update(parameters)

    # print('number of group', len(group_info_list))
    nb_group_with_valid_score = sum(
        [1 for grp in group_info_list if grp['confidence_score'] < parameters["confidence_score_threshold"]])
    # print( f'nb group that have a completeness > of {completeness_threshold} is {nb_group_with_valid_score}')
    # print(f'cleavage site per seq: {cleavages_sites_per_seq}' )

    if general_info['nb_seq_with_annotation'] == 1:
        # print('SINGLE....')
        aln_validity = 'single annotated polyprotein'

    elif sum([1 for grp in group_info_list if grp['group_completeness'] == 100]) == len(group_info_list):
        # print('PERFECT!!!')
        # print("All group have a good confidence score")
        # print('We can propagate')
        aln_validity = 'All cleavage sites groups have a good score'

    elif min(cleavages_sites_per_seq) <= nb_group_with_valid_score:
        # print("GOOD ENOUGH...")
        # print('We can propagate')
        aln_validity = 'Number of valid cleavage_site groups >= minimum number of cleavage site/cds annotated'

    elif nb_group_with_valid_score > 0:
        # print("OK.. ")
        # print('At least one group is valid')
        aln_validity = 'At least one cleavage site group is valid'
    elif nb_group_with_valid_score == 0:
        # print('NOPE...')
        # print('No cleavage site groups are valid')
        aln_validity = 'None of the cleavage site groups are valid'

    aln_dict_stat = {'nb_unannotated_seq': general_info['nb_sequences'] - general_info['nb_seq_with_annotation'],
                     "nb_group_of_cleavage_sites": len(group_info_list),
                     "max_cleavage_site_per_sequence": max(cleavages_sites_per_seq),
                     "min_cleavage_site_per_sequence": min(cleavages_sites_per_seq),
                     "aln_validity": aln_validity,
                     "nb_group_with_valid_score": nb_group_with_valid_score
                     }
    aln_dict_stat.update(general_info)

    if aln_csv_writer:
        aln_csv_writer.writerow(aln_dict_stat)

    if group_site_csv_writer:
        for group_info in group_info_list:
            group_info.update(general_info)
            group_info['aln_validity'] = aln_validity

            group_info['cleavages_sites_per_seq'] = None
            print(f'group {group_info["group_index"]} score: {group_info["confidence_score"]}')
            group_site_csv_writer.writerow(group_info)


def initiate_ouput(cs_group_stat_file, output_stat_aln):

    # File stat on cleavage site
    cs_header = ["cluster_nb",
                 "nb_sequences",
                 "nb_seq_with_annotation",
                 "window",
                 "nb_protein_in_group",
                 "nb_of_cleavage_site",
                 "nb_cleavage_from_the_same_protein",
                 "standard_deviation",
                 "round_std",
                 "group_completeness",
                 "completeness_category",
                 "average_position_in_aln",
                 "position_max",
                 "position_min",
                 "group_index",
                 'aln_validity',
                 "confidence_score",
                 "confidence_score_threshold",
                 'cluster_with_isoforms_splitted',
                 'evalue',
                 'coverage',
                 'inflation',
                 "cleavages_sites_per_seq"]

    handle_cs_out = open(cs_group_stat_file, "w")
    group_site_csv_writer = csv.DictWriter(handle_cs_out, fieldnames=cs_header, delimiter='\t')
    group_site_csv_writer.writeheader()

    # File stat on the Alignment
    aln_header = ['cluster_nb',
                  'nb_sequences',
                  'nb_seq_with_annotation',
                  'nb_unannotated_seq',
                  "nb_group_of_cleavage_sites",
                  "window",
                  'nb_group_with_valid_score',
                  "max_cleavage_site_per_sequence",
                  "min_cleavage_site_per_sequence",
                  "aln_validity",

                  "confidence_score_threshold",
                  'cluster_with_isoforms_splitted',
                  'evalue',
                  'coverage',
                  'inflation',
                  "cleavages_sites_per_seq"]

    handle_aln_out = open(output_stat_aln, "w")
    aln_csv_writer = csv.DictWriter(handle_aln_out, fieldnames=aln_header, delimiter='\t')
    aln_csv_writer.writeheader()
    file_handles = [handle_aln_out, handle_cs_out]

    return group_site_csv_writer, aln_csv_writer, file_handles


def parse_aln_directory_name(alignment_dir):
    dir_info = {}
    if alignment_dir.endswith('splitted'):
        dir_info['cluster_with_isoforms_splitted'] = True
    else:
        dir_info['cluster_with_isoforms_splitted'] = False

    re_result = re.search("_evalue_([\de-]+)coverage(\d+)_I(\d_{0,1}\d*)", alignment_dir)
    dir_info['evalue'] = "NA" if not re_result else re_result.group(1)
    dir_info['coverage'] = "NA" if not re_result else re_result.group(2)
    dir_info['inflation'] = "NA" if not re_result else re_result.group(3)

    return dir_info


def main():
    re_result = re.search("cluster(\d+).aln", alignment_file)
    cluster_nb = "NA" if not re_result else re_result.group(1)
    print("PROCESS of CLUSTER ", cluster_nb)
    taxon_prot_ids, seq_aln_dict = parse_alignment_file(alignment_file)

    for window in windows:
        cds_list = getCdsObject(taxon_prot_ids, taxonomy_file, interpro_domains_file, sp_treshold)

        if sum((1 for cds in cds_list if cds.polyprotein)) == 0:
            logging.warning(f'{alignment_file} has no identified polyprotein with annotation')
            return

        if len(cds_list) != len(seq_aln_dict):
            # print('len(cds_list) !=  len(seq_aln_dict) ')
            # print(len(cds_list),len(seq_aln_dict) )
            logging.critical(
                f'Not all cds from the alignment have been retrieved in {alignment_file}')
            return

        for cds in cds_list:
            cds.aligned_sequence = seq_aln_dict[cds.protein_id]
            cds.cleavage_site_positions = {s.start_aa(cds): s for s in cds.cleavage_sites}
            cds.aln_list = convertAlignmentToList(cds.aligned_sequence)

        group_info_list, cs_group_list = analyse_cleavage_site_groups(cds_list, window)

        for group_info in group_info_list:
            group_of_cs = cs_group_list[group_info["group_index"]]
            if group_info['confidence_score'] < confidence_score_threshold:
                propagate_cleavage_sites(group_of_cs, group_info, cds_list, window)
            #     group_info['valid'] = True
            # else:
            #     group_info['valid'] = False

        if not aln_csv_writer:
            print(f'VISUALISATION WITH WINDOW {window}')

            for cds in cds_list:
                print(cds.segment.record.annotations['taxonomy'])

            for group in group_info_list:
                print(group)

            visualisation_of_processed_aln(cds_list, alignment_file,
                                           group_info_list, display_line_size)

        if aln_csv_writer:
            parameters["cluster_nb"] = cluster_nb
            parameters['window'] = window
            compute_alignment_stat(group_info_list, parameters, cds_list,
                                   aln_csv_writer, group_site_csv_writer)


if __name__ == '__main__':

    # logging.basicConfig(filename='log/alignment_analysis.log', level=logging.INFO)
    # alignment_file= "data/alignment/Retro-transcribing_viruses_1e-5_coverage90_I3/seq_cluster16.aln"

    try:
        # 'data/alignment/Viruses_1e-5_coverage90_I2/seq_cluster1037.aln'
        alignment_file_or_dir = sys.argv[1]
    except:
        alignment_file_or_dir = 'data/alignment/Viruses/Viruses_evalue_1e-60coverage60_I1_8/seq_cluster0.aln'
        # alignment_file_or_dir = 'data/alignment/Viruses/RefSeq_download_date_2018-07-21/Viruses_evalue_1e-40coverage40_I2/seq_cluster30.aln'
    try:
        # 'data/alignment/Viruses_1e-5_coverage90_I2/seq_cluster1037_stat.csv'
        windows_input = sys.argv[2]
    except:
        windows_input = "30"

    try:
        cs_group_stat_file = sys.argv[3]
        output_stat_aln = sys.argv[4]
    except:
        cs_group_stat_file = None
        output_stat_aln = None
        cs_group_stat_file = "test/test_cs_result.csv"
        output_stat_aln = "test/test_aln_result.csv"
    try:
        taxonomy_file = sys.argv[5]  # "data/taxonomy/taxonomy_virus.txt"
    except:
        taxonomy_file = "data/taxonomy/taxonomy_virus.txt"

    try:
        interpro_domains_file = sys.argv[6]
    except:
        interpro_domains_file = None

    sp_treshold = 90
    display_line_size = 75  # 140/2
    confidence_score_threshold = 4

    parameters = {"confidence_score_threshold": confidence_score_threshold}

    if path.isdir(alignment_file_or_dir):
        alignment_files = (path.join(alignment_file_or_dir, f)
                           for f in listdir(alignment_file_or_dir) if f.endswith('.aln'))
        info_dir = parse_aln_directory_name(alignment_file_or_dir)
        parameters.update(info_dir)

    elif path.isfile(alignment_file_or_dir) and alignment_file_or_dir.endswith('.aln'):
        alignment_files = [alignment_file_or_dir]
    else:
        raise ValueError('file or directory provided is not correct')

    print("OUTPUT:")
    print(cs_group_stat_file)
    print(output_stat_aln)

    # window includ the cleavage in the middle
    # cleavage sites have a length of 2
    # then a window will always be even and >= 2
    # consequently we add 1 to odd window

    windows = {int((int(w)+int(w) % 2)) for w in windows_input.split(' ')}  # 10,20,30
    assert min(windows) >= 0
    windows = list(windows)
    windows.sort()
    print('window used to analyse cleavage sites',  windows)

    if cs_group_stat_file:
        print("INITIATE OUTPUT FILES")
        group_site_csv_writer, aln_csv_writer, file_handles = initiate_ouput(
            cs_group_stat_file, output_stat_aln)
    else:
        group_site_csv_writer = False
        aln_csv_writer = False
        file_handles = []

    for alignment_file in alignment_files:
        print(alignment_file)
        main()  # variable needed in main() are global

    for fl in file_handles:
        fl.close()

#!/usr/bin/env python3
import taxonomy as tax
import viral_genome_classes as obj
import viruses_statistics as stat
import parser_interpro_results as do
import visualisation_alignment as view_aln
import visualisation_protein as view_prot
import multiple_alignment_analysis as analysis

from  os import path, listdir
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


def split_cluster(cds_list):
     # includer when cds has cds obj in sub_prot attr: meaning that includer
     # cds include in their sequence a smaller cds that start or finish at a different position
    parental_cds = {}
    correponding_sub_protein = []
    parental_cds_lvl = {}
    level = {} # lvl of inclusion we have in the cluster. in most of the case it's 1
                #but some Retro virus have 2 lvl (so a stack of 3 proteins sharing a common part)

    for cds in cds_list:
        if cds.sub_prot:
            # print(f'cds {cds.protein_id} has {len(cds.sub_prot)} sub prot')

            sub_cds_in_cluster = [sub for sub in cds.sub_prot if sub in cds_list]

            if sub_cds_in_cluster:
                # print(f'   {len(sub_cds_in_cluster)} are in the cluster')
                parental_cds[cds] = sub_cds_in_cluster
                parental_cds_lvl.setdefault(len(sub_cds_in_cluster),[]).append(cds)
    if not parental_cds:
        return
    # print('---------------')
    # print(parental_cds_lvl)
    # print(parental_cds)
    level_max = max(parental_cds_lvl)
    # print('level max', level_max)
    # print('prot lvl max is ', parental_cds_lvl[level_max])
    # print('nb prot lvl max', len(parental_cds_lvl[level_max]))
    group_level_cds= {level_max: list(parental_cds_lvl[level_max])}
    for cds in parental_cds_lvl[level_max]:
        # print(parental_cds_lvl)
        # print([cds])
        for sub_prot in parental_cds[cds]:
            # print(f'sub prot has {len(sub_prot.sub_prot)}')
            sub_prot_lvl = 0 if sub_prot not in parental_cds else len(parental_cds[sub_prot])
            group_level_cds.setdefault(sub_prot_lvl,[]).append(sub_prot)



    ## Check that length min max of each categoty if not overlapping
    length_distribution_prev = 0
    size_threshold_for_splitting = []
    for lvl in range(level_max+1):
        cds_lengths = [len(cds) for cds in group_level_cds[lvl]]
        # print('previous', length_distribution_prev)
        # print(cds_lengths)

        if length_distribution_prev < min(cds_lengths)-percentage_len*min(cds_lengths):
            # print('OK')
            if length_distribution_prev >0:
                size_threshold_for_splitting.append((length_distribution_prev+ min(cds_lengths))/2)

        else:
            logging.warning('problem: size of group is overlapping. Cluster is not splitted')
            return

        length_distribution_prev = max(cds_lengths)

    # print(size_threshold_for_splitting)

    cds_group = {}
    for cds in cds_list:
        for i,len_threshold in enumerate(size_threshold_for_splitting):
            if len(cds) < len_threshold:
                cds_group.setdefault(i, []).append(cds)
                sequence_stored = True
                break
            else:
                sequence_stored = False

        if sequence_stored is False:
            #store the sequence that are bigger than the last threshold
            cds_group.setdefault(len(size_threshold_for_splitting), []).append(cds)


    singleton_group_index_keys = []
    for i, group in cds_group.items():
        # print(i)
        # print([len(cds) for cds in group])
        if len(group) == 1:
            logging.warning('The splitting would make a singleton cluster.. the sequence is merged in the group of long sequence if possible')
            singleton_group_index_keys.append(i)

    #MERGE POTENTIAL SINGLETON GROUP CREATED BY THE SPLITTING
    for singleton_group_index_key in singleton_group_index_keys:
        if singleton_group_index_key < len(cds_group) -1:
            # print('MERGE HIGH LEN')
            cds_group[singleton_group_index_key+1].append(cds_group[singleton_group_index_key][0])
        elif singleton_group_index_key == len(cds_group) -1 and len(cds_group) > 1:
            # print('MERGE LOW LEN')
            cds_group[singleton_group_index_key-1].append(cds_group[singleton_group_index_key][0])
        else:
            logging.warning('problem: the merging is not possible.  Cluster is not splitted')
            return
        del cds_group[singleton_group_index_key]

    #CREATION OF THE LINES OF THE NEW CLUSTERS.
    if len(cds_group) < 2:
        return
    else:
        return create_cluster_lines(cds_group)


def create_cluster_lines(cds_group):
    lines={}
    for i, group in cds_group.items():
        # print()
        line_list = [ f'{cds.segment.taxon_id}|{cds.protein_id}|{int(len(cds)/3)}' for cds in group]
        line = '\t'.join(line_list)
        lines[i] = line

    return lines


def get_taxon_protein_ids(line):
    taxon_protein_ids = {}
    line = line.rstrip()
    line_list = line.split('\t')
    for sequence_header in line_list:
        ## 	138950|NP_041277.1|2210	138951|NP_040760.1|2195	1826059|YP_009246449.1|2210
        taxon_id, protein_id, len = sequence_header.split('|')
        taxon_protein_ids.setdefault( taxon_id, set()).add(protein_id)
    return taxon_protein_ids


def main():
    count_new_cluster = 0
    old_cluster_splitted = 0
    fl_out = open(cluster_filter_out, 'w')
    with open(cluster_file, 'r') as fl:
        for i, l in enumerate(fl):
            taxon_protein_ids = get_taxon_protein_ids(l)
            cds_list = analysis.getCdsObject(taxon_protein_ids, taxonomy_file, gff_file, sp_treshold)
            # print("=="*20)
            # print("CLUSTER ", i)
            new_lines = None
            if any((True for cds in cds_list if cds.sub_prot)):
                new_lines = split_cluster(cds_list)
            # else:
            #     print(f"No sub prot in cluster {i}")

            if new_lines:
                old_cluster_splitted += 1
                assert sum((len(n_l.split('\t')) for n_l in new_lines.values())) == len(l.split('\t')), 'sum elments in new cluster is different of nb element in old cluster!'
                for new_l in new_lines.values():
                    count_new_cluster += 1
                    fl_out.write(f'{new_l}\n')
            else:
                fl_out.write(l)

    print(f'{old_cluster_splitted} clusters have been splitted into {count_new_cluster} new clusters')
    fl_out.close()

if __name__ == '__main__':


    try:
        cluster_file=  sys.argv[1]
        cluster_filter_out =  sys.argv[2]
    except:
        cluster_file = "data/alignment/ssRNA_viruses/RefSeq_download_date_2018-07-21/ssRNA_viruses_evalue_1e-20coverage60_I1_8/clusters_with_identified_polyprotein.out"
        # cluster_file = 'test.csv'
        cluster_filter_out = "test/filter_cluster"

    taxonomy_file = "data/taxonomy/taxonomy_virus.txt"
    gff_file = None
    sp_treshold = 90
    percentage_len = 0.1

    main()

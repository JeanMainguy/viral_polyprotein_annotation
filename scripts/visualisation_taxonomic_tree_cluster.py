import taxonomic_homogeneity_checking as homo
from ete3 import Tree
import logging, sys
from os import path
import re, csv
from time import clock

def constructETEtree(node_dict, count_leaves, genome_id_to_name, color, current_node, children, cluster_info={}):

    for child in children:
        count = count_leaves.setdefault(child, child)
        if count == child:
            name = color['leaf']+genome_id_to_name[child]+color['end']
            if cluster_info and child in cluster_info["genome_ids"]:
                    name = color['cluster_leaf']+genome_id_to_name[child]+color['end']
        else:
            name=str(count)+':'+color['node']+child+color['end']
        child_node = current_node.add_child(name = name )
        if child in node_dict:
            if cluster_info and child in cluster_info['unvalid_branches']:
                continue
            next_children = node_dict[child]
            constructETEtree(node_dict, count_leaves, genome_id_to_name, color, child_node, next_children, cluster_info)


def displayTree(node_dict, root_name, count_leaves, genome_id_to_name):
    color = {'end':'\033[0m',
            'leaf':'\033[94m',
            'node':'\033[93m',
            'cluster_leaf':'\033[92m' }
    t = Tree()
    root = t.add_child(name=str(count_leaves[root_name])+':'+color['node']+root_name+color['end'])
    children = node_dict[root_name]
    constructETEtree(node_dict,count_leaves, genome_id_to_name, color, root, children)
    print(t.get_ascii(show_internal=True))
    return t

def displayClusterTree(node_dict, cluster_info, count_leaves, genome_id_to_name, leaves_taxonomy, show_unvalid_branches=False):
    #cluster info is a dict with information on the cluster
    # with inject into the constructETEtree fct: the shared_taxonomy of the cluster, the genome ids, valid branches

    color = {'end':'\033[0m',
            'leaf':'\033[94m',
            'node':'\033[93m',
            'cluster_leaf':'\033[91m' }
    t = Tree()
    root_name = cluster_info["shared_taxonomy"].split(';')[-1]

    root = t.add_child(name=str(count_leaves[root_name])+':'+color['node']+root_name+color['end'])
    children = node_dict[root_name]
    # print('children node', {c for c in children if not c.isdigit() })
    # print('valid branches',set(cluster_info['valid_branches'].split('|') ))

    cluster_info['unvalid_branches'] = set({c for c in children if not c.isdigit() }) - set(cluster_info['valid_branches'].split('|'))
    # print('unvalid branch', cluster_info['unvalid_branches'] )
    if show_unvalid_branches:
        cluster_info['unvalid_branches']  = set()

    
    constructETEtree(node_dict,count_leaves,genome_id_to_name, color, root, children, cluster_info)
    print(t.get_ascii(show_internal=True))
    for k in cluster_info:
        print(k, cluster_info[k])
    return t

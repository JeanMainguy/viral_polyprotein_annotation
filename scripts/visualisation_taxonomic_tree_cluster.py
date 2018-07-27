import taxonomic_homogeneity_checking as homo
from ete3 import Tree, TreeStyle, Tree, TextFace, add_face_to_node, NodeStyle, SeqMotifFace
import logging, sys
from os import path
import re, csv
from time import clock

def constructETEtree(node_dict, count_leaves, genome_id_to_name, color, current_node, children, cluster_info, node_styles):

    for child in children:
        count = count_leaves.setdefault(child, child)
        if count == child:
            style = node_styles['leaf']
            type = 'leaf'
            name = color['leaf']+genome_id_to_name[child]+color['end']
            if cluster_info and child in cluster_info["genome_ids"]:
                    style = node_styles['cluster_leaf']
                    type = 'cluster_leaf'
                    name = color['cluster_leaf']+genome_id_to_name[child]+color['end']
        else:
            style = node_styles['node']
            type = 'node'
            name=str(count)+':'+color['node']+child+color['end']
        child_node = current_node.add_child(name = name )
        child_node.add_features(nodetype=type)
        child_node.set_style(style)
        if child in node_dict:
            if cluster_info and child in cluster_info['unvalid_branches']:
                continue
            next_children = node_dict[child]
            constructETEtree(node_dict, count_leaves, genome_id_to_name, color, child_node, next_children, cluster_info, node_styles)


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

def displayClusterTree(node_dict, cluster_info, count_leaves, genome_id_to_name, leaves_taxonomy, show_unvalid_branches=False, color={}, node_styles={}):
    #cluster info is a dict with information on the cluster
    # with inject into the constructETEtree fct: the shared_taxonomy of the cluster, the genome ids, valid branches
    t = Tree()
    root_name = cluster_info["shared_taxonomy"].split(';')[-1]

    root = t.add_child(name=str(count_leaves[root_name])+':'+color['node']+root_name+color['end'])
    children = node_dict[root_name]
    # print('children node', {c for c in children if not c.isdigit() })
    # print('valid branches',set(cluster_info['valid_branches'].split('|') ))

    # print('unvalid branch', cluster_info['unvalid_branches'] )
    cluster_info['unvalid_branches'] = set({c for c in children if not c.isdigit() }) - set(cluster_info['valid_branches'].split('|'))
    if show_unvalid_branches:
        cluster_info['unvalid_branches']  = set()


    constructETEtree(node_dict,count_leaves,genome_id_to_name, color, root, children, cluster_info, node_styles)
    print(t.get_ascii(show_internal=True))
    for k in cluster_info:
        print(k, cluster_info[k])
    return t

def my_layout(node):
    # print(dir(node))
    # if node.nodetype == 'leaf':
    #     node_face = TextFace(node.name, tight_text=True)
    #     node_face.background.color = "LightGreen"
    #
    # elif node.nodetype == 'node':
    if not node.is_leaf():
        node_face = TextFace(node.name, tight_text=True)
        add_face_to_node(node_face, node, column=0)

if __name__ == '__main__':

    cluster_file = sys.argv[1]
    taxonomy_file = "data/taxonomy/taxonomy_virus.txt"
    alternative_taxon_id_file =  "data/taxonomy/heterogeneous_taxon_id_taxonomy_virus.txt"
    unclassified_term_file = "data/taxonomy/unclassified_terms.txt"

    display_ETE3 = True
    # we retrieve the coverage value and the inflation from the file name
    re_result = re.search("coverage(\d+)_I([\d_]+).out", cluster_file)
    #if the regex doesn't work an expception will be raised
    coverage = re_result.group(1)
    inflation = re_result.group(2)


    unclassified_term_list = homo.getListFromFile(unclassified_term_file)
    node_dict, count_leaves, leaves_taxonomy, genome_id_to_name = homo.constructTaxonTree(taxonomy_file, unclassified_term_list)

    iter_clusters = homo.clusterFileParser(cluster_file)
    alternative_taxon_id =homo.getAlternativeTaxonId(alternative_taxon_id_file)

    # csv_writer = initiateOuputFile(cluster_homogeneity_file)
    # print(alternative_taxon_id)
    #color for the terminal ...
    color = {'end':'\033[0m',
            'leaf':'\033[94m',
            'node':'\033[93m',
            'cluster_leaf':'\033[91m' }

    node_styles={}
    #LEAF
    nstyle = NodeStyle()
    nstyle["fgcolor"] = "blue"
    nstyle["size"] = 5

    node_styles['leaf'] = nstyle

    #NODE
    nstyle = NodeStyle()
    nstyle["fgcolor"] = "yellow"
    nstyle["size"] = 5
    node_styles['node'] = nstyle

    #Cluster leaf
    nstyle = NodeStyle()
    nstyle["fgcolor"] = "red"
    nstyle["size"] = 5
    # nstyle['bgcolor'] = 'FireBrick'
    node_styles['cluster_leaf'] = nstyle

    if display_ETE3:
        for k in color:
            color[k] = ''

    for i, cluster_info in enumerate(iter_clusters):
        #cat data/taxonomy/taxonomy_virus.txt | cut -f3 | cut -d';' -f2 | sort | uniq | more
        dict_info = homo.computeHomogeneity(cluster_info, unclassified_term_list, leaves_taxonomy, alternative_taxon_id, count_leaves)

        dict_info["inflation"] = inflation
        dict_info["coverage"] = coverage
        cluster_info.update(dict_info)
        # if len(cluster_info['genome_ids']) < 80:
        #     print("HOMOGENEITY", dict_info["homogeneity"])
        if dict_info['number_of_genome'] > 30:
            continue
        t = displayClusterTree(node_dict, cluster_info, count_leaves, genome_id_to_name, leaves_taxonomy, False, color, node_styles )
        # tree style
        ts = TreeStyle()
        ts.show_leaf_name = False
        ts.layout_fn = my_layout
        # ts.mode = "c"
        # ts.arc_start = 0 # 0 degrees = 3 o'clock
        # ts.arc_span = 180
        for node in t.search_nodes(nodetype="cluster_leaf"):
            t_face = TextFace(node.name, ftype='Verdana', fsize=10, fgcolor='FireBrick', penwidth=0, fstyle='normal', tight_text=True, bold=False)
            print(node)
            node.add_face(t_face, column=0)
            print(node.nodetype)

        for node in t.search_nodes(nodetype="leaf"):
            t_face = TextFace(node.name, ftype='Verdana', fsize=10, fgcolor='Black', penwidth=0, fstyle='normal', tight_text=True, bold=False)
            print(node)
            node.add_face(t_face, column=0)

        for node in t.search_nodes(nodetype="node"):
            t_face = TextFace(node.name, ftype='Verdana', fsize=10, fgcolor='Black', penwidth=0, fstyle='normal', tight_text=True, bold=False)
            print(node)
            node.add_face(t_face, column=0)
        print(t.search_nodes(nodetype="cluster_leaf"))
        t.show(tree_style=ts)
        input()

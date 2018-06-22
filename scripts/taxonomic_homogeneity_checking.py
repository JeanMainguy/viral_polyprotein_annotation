from ete3 import Tree
import logging, sys
from os import path
import re, csv
from time import clock
import visualisation_taxonomic_tree_cluster as view

# PythonDecorators/decorator_without_arguments.py
class timer(object):

    def __init__(self, f):
        """
        If there are no decorator arguments, the function
        to be decorated is passed to the constructor.
        """
        self.f = f

    def __call__(self, *args, **kwargs):
        """
        The __call__ method is not called until the
        decorated function is called.
        """
        START_TIME = clock()
        result = self.f(*args, **kwargs)
        print(self.f.__name__, 'run time:%6.2f'%(clock()-START_TIME), 'seconds.');
        print('-'*20)
        return result
@timer
def constructTaxonTree(taxonomy_file, unclassified_term_list):


    node_dict = {} #name:[child1, child2....]
    count_leaves = {} #"node name":int(number of leaves)
    leaves_taxonomy = {} # genome_id:string taxonomy
    genome_id_to_name = {}
    #438782	Abaca bunchy top virus	Viruses;ssDNA viruses;Nanoviridae;Babuvirus	1	/mirr[...]gbff.gz
    with open(taxonomy_file, 'r') as taxon_handle:
        for l in taxon_handle:
            genome_id, genome_name, taxonomy, genetic_code, gbff_file = l.split("\t")

            taxonomy = taxonomy.split(';')

            genome_id_to_name[genome_id] = genome_name
            #We don't want to count the unclassified taxonomy node
            #For example the virus "Tai Forest alphavirus" has the folowing taxonomy : Viruses;[...];Alphavirus;unclassified Alphavirus
            # we remove the last node and associate it with the name of the virus in the dictionnary genome_id_to_name
            for unclassified_term in unclassified_term_list:
                if len(taxonomy)>2 and taxonomy[-1].startswith(unclassified_term):
                    genome_id_to_name[genome_id] = "{} ({})".format(genome_name, taxonomy.pop())
            leaves_taxonomy[genome_id] = tuple(taxonomy)
            root = taxonomy[0]
            node_dict.setdefault(root, set()) # add the root to the dict if not there


            updateTree(taxonomy,node_dict, count_leaves, genome_id)

    return node_dict, count_leaves, leaves_taxonomy, genome_id_to_name


def updateTree(taxonomy, node_dict, count_leaves, child):
    try:
        last_taxon = taxonomy.pop()
    except IndexError:
        return

    #update the node leaf count. He has one more leaf
    count_leaves[last_taxon] = count_leaves.setdefault(last_taxon, 0) + 1

    #add child to the last taxon if it is not already there
    node_dict.setdefault(last_taxon, set()).add(child)

    updateTree(taxonomy, node_dict, count_leaves, last_taxon)


def clusterFileParser(cluster_file):
    # 1093958|YP_004901696.1|362.0|None	759390|YP_003620408.1|362.0|None	2021667|YP_009408633.1|357.0|None
    with open(cluster_file, 'r') as cluster_handle:
        for i, cluster_line in enumerate(cluster_handle):


            cluster_elements = cluster_line.rstrip().split("\t")
            cluster_info = {'genome_ids':[], 'nb_polyprotein':0, 'nb_protein':len(cluster_elements), "cluster_id":i}
            # 1774200|YP_009216586.1|357.0|None
            for element in cluster_elements:
                genome_ids, protein_id, length, peptide_annotation = element.split('|')
                if peptide_annotation == 'Peptide':
                    cluster_info['nb_polyprotein'] += 1
                cluster_info['genome_ids'].append(element.split('|')[0])
            # genome_ids = {element.split('|')[0] for element in cluster}
            yield cluster_info


def getSharedTaxonomy(taxonomy_set):

    # valid_branches = set()
    taxonomy_iter = iter(taxonomy_set)
    current_shared_taxonomy = next(taxonomy_iter)

    for new_taxonomy in taxonomy_iter:
        if  current_shared_taxonomy[-1] in new_taxonomy:
            continue

        current_taxonomy_iter = iter(current_shared_taxonomy)
        #We want to get the last smallest share taxonomy between the two taxonomy
        # we take the smallest one as itererator and we loop on the other one
        smallest_taxonomy_iter, biggest_taxonomy_iter = (iter(new_taxonomy), iter(current_shared_taxonomy)) if len(new_taxonomy) <= len(current_shared_taxonomy) else (iter(current_shared_taxonomy), iter(new_taxonomy))


        for i, taxon in enumerate(biggest_taxonomy_iter):
            # print("taxon tested to be in current_shared_taxonomy ",taxon )
            if taxon != next(smallest_taxonomy_iter, False):
                # print(taxon, 'not in current ..')
                current_shared_taxonomy = new_taxonomy[:i]
                # print("new shared tax is set", ";".join(current_shared_taxonomy))
                break
    # last_common_node = current_shared_taxonomy[-1]
    # # print(last_common_node)
    # # print(current_shared_taxonomy)
    # for tax in taxonomy_set:
    #     i_last = tax.index(last_common_node)
    #     try:
    #         valid_branches.add(tax[i_last+1])
    #     except:
    #         valid_branches = {last_common_node}
    #         break

    # print("LAST TAX")
    # print(current_shared_taxonomy)
    # print('valid branch', valid_branches )
    # # print('NB GENOmE in cluster', len(cluster))
    # # print('NB GENOmE in subtree', sum( [count_leaves[t] for t in valid_branches]))
    # print("--------------")
    # subtree_nb_leaves =  sum( [count_leaves[t] for t in valid_branches])

    return current_shared_taxonomy #, valid_branches

def getValidBranches(current_shared_taxonomy, taxonomy_set):
    last_common_node = current_shared_taxonomy[-1]

    valid_branches = set()
    for tax in taxonomy_set:
        i_last = tax.index(last_common_node)
        if len(tax) > i_last+1:
            valid_branches.add(tax[i_last+1])
    if not valid_branches:
        valid_branches.add(last_common_node)
    return valid_branches



def getAlternativeTaxonId(alternative_taxon_id_file):
    alternative_taxon_id={} # key: taxon id value : set of alternative taxon id
    with open(alternative_taxon_id_file, 'r') as handle:
        for l in handle:
            #429564;77811	/mirror/ncbi/current/g...
            set_taxon_id = set(l.split("\t")[0].split(';'))
            for i, taxon in enumerate(set_taxon_id):
                alternative_taxon_id[taxon] = set_taxon_id - {taxon}
    return alternative_taxon_id


def get_filtered_taxonomy_set(genomes_ids, unclassified_term_list, leaves_taxonomy, alternative_taxon_id):
    '''
    Return set of different taxonomy found in the cluster

    '''
    taxonomy_set = set()
    unclassified_taxonomy_set = set()
    unclassified_genome = 0
    # unfound_genome = 0
    for i, taxon_id in enumerate(genomes_ids):
        taxonomy = None
        try:
            taxonomy = leaves_taxonomy[taxon_id]
        except:
            if taxon_id in alternative_taxon_id:
                for alternative in alternative_taxon_id[taxon_id]:
                    if alternative in leaves_taxonomy:
                        taxonomy = leaves_taxonomy[alternative]
                        genomes_ids[i] = alternative

            if not taxonomy:

                # logging.warning("taxon id not found in taxonomy file or in alternative taxon id file: "+taxon_id)
                raise NameError("taxon id not found in taxonomy file or in alternative taxon id file: ",taxon_id)
                # unfound_genome += 1
                continue
        # Check if the first node after the root is an unclassified node
        # if it is this genome is not take into account in the homogeneity
        # print(taxon_id, taxonomy)
        for unclassified_term in unclassified_term_list:
            if taxonomy[1].startswith(unclassified_term):
                unclassified_genome += 1
                unclassified_taxonomy_set.add(taxonomy)
        taxonomy_set.add(taxonomy)


    return taxonomy_set, unclassified_genome, unclassified_taxonomy_set


def initiateOuputFile(file_output):

    handle_stat = open(file_output, "w")
    #segment.taxon_id, cds.protein_id, has_peptide, len(cds.peptides), len(cds.cleavage_sites), is_sub_protein
    header = ['cluster_id',
                    'nb_of_protein',
                    'number_of_genome',
                    'number_of_genome_in_valid_branches',
                    'unclassified_genome',
                    'nb_polyprotein',
                    "number_leaves_in_subtree",
                    "homogeneity",
                    "shared_taxonomy",
                    'valid_branches',
                    'inflation',
                    'coverage',
                    "singleton",
                    "unclassified_cluster",
                    "first_node"]

    csv_writer = csv.DictWriter(handle_stat, fieldnames=header, delimiter='\t')
    csv_writer.writeheader()

    return csv_writer


def computeHomogeneity(cluster_info, unclassified_term_list, leaves_taxonomy, alternative_taxon_id, count_leaves):

    taxonomy_set, unclassified_genome,  unclassified_taxonomy_set= get_filtered_taxonomy_set(cluster_info['genome_ids'], unclassified_term_list, leaves_taxonomy,alternative_taxon_id)


    if unclassified_genome == len(cluster_info['genome_ids']): # the all proteins of the cluster are coming from unclassified genomes
        # assert  all([True for t in taxonomy_set if t[1] in unclassified_term_list]), "Some taxonomy are not unclassified while there are predicted as"
        unclassified_cluster = True
        effective_taxonomy_set = taxonomy_set
        assert taxonomy_set == unclassified_taxonomy_set
    else:
        unclassified_cluster = False
        effective_taxonomy_set = taxonomy_set - unclassified_taxonomy_set
    current_shared_taxonomy = getSharedTaxonomy(effective_taxonomy_set)
    valid_branches = getValidBranches(current_shared_taxonomy, effective_taxonomy_set)
    ##Get nb genome from the cluster that are included in the valid branches
    nb_genome_in_valid_branch = 0
    nb_genome = len(set(cluster_info['genome_ids']))
    for genome_id in set(cluster_info['genome_ids']):
        if any([True for node in leaves_taxonomy[genome_id] if node in valid_branches]):
            nb_genome_in_valid_branch += 1



    subtree_nb_leaves = sum([count_leaves[t] for t in valid_branches])
    homogeneity = (nb_genome_in_valid_branch)/subtree_nb_leaves
    # print(nb_genome,"/",subtree_nb_leaves)
    try:
        first_node = current_shared_taxonomy[1]
    except IndexError: # the current shared tax is only Viruses
        first_node = current_shared_taxonomy[0]
    if unclassified_cluster:
        first_node = 'unclassified'

    current_shared_taxonomy = ';'.join(current_shared_taxonomy)
    valid_branches = '|'.join(valid_branches)

    assert (0 < homogeneity <= 1),"homogeneity value is over 1: %f"  % homogeneity

    dict_info = {'cluster_id':cluster_info['cluster_id'],
                    'nb_of_protein':cluster_info['nb_protein'],
                    'number_of_genome_in_valid_branches':nb_genome_in_valid_branch,
                    'number_of_genome':nb_genome,
                    'unclassified_genome':unclassified_genome,
                    'nb_polyprotein':cluster_info['nb_polyprotein'],
                    "number_leaves_in_subtree":subtree_nb_leaves,
                    "homogeneity":homogeneity,
                    "shared_taxonomy":current_shared_taxonomy,
                    'valid_branches':valid_branches,
                    "unclassified_cluster":unclassified_cluster,
                    "first_node":first_node
                    }
    return dict_info

@timer
def writeHomogeneity(iter_clusters, inflation, coverage, unclassified_term_list, leaves_taxonomy, alternative_taxon_id, count_leaves, genome_id_to_name, node_dict):
    for i, cluster_info in enumerate(iter_clusters):
        #cat data/taxonomy/taxonomy_virus.txt | cut -f3 | cut -d';' -f2 | sort | uniq | more
        dict_info = computeHomogeneity(cluster_info, unclassified_term_list, leaves_taxonomy, alternative_taxon_id, count_leaves)

        dict_info["inflation"] = inflation
        dict_info["coverage"] = coverage
        csv_writer.writerow(dict_info)
        # cluster_info.update(dict_info)
        # if len(cluster_info['genome_ids']) < 80:
        #     print("HOMOGENEITY", dict_info["homogeneity"])
        # displayClusterTree(node_dict, cluster_info, count_leaves, genome_id_to_name, leaves_taxonomy)

def getListFromFile(file):
    list_line = []
    with open(file, "r") as hd:
        for l in hd:
            list_line.append(l.rstrip())
    return list_line


if __name__ == '__main__':


    cluster_file = sys.argv[1]
    cluster_homogeneity_file = sys.argv[2]
    taxonomy_file = sys.argv[3]  #"data/taxonomy/taxonomy_virus.txt"
    alternative_taxon_id_file = sys.argv[4] # "data/taxonomy/heterogeneous_taxon_id_taxonomy_virus.txt"
    unclassified_term_file = sys.argv[5] # "data/taxonomy/unclassified_terms.txt"
    try:
        visualisation = True if sys.argv[6] == 'visualisation' else False
    except:
        visualisation = False


    # unclassified_nodes =['unassigned viruses',
    #                 'unclassified archaeal viruses',
    #                 'unclassified bacterial viruses',
    #                 'unclassified RNA viruses',
    #                 'unclassified viroids',
    #                 'unclassified virophages',
    #                 'unclassified viruses',
    #                 'environmental samples',
    #                 "Virus families not assigned to an order"]

    # we retrieve the coverage value and the inflation from the file name
    re_result = re.search("coverage(\d+)_I([\d_]+).out", cluster_file)
    #if the regex doesn't work an expception will be raised
    coverage = re_result.group(1)
    inflation = re_result.group(2)
    # base_cluster_file = path.basename(cluster_file)
    # file_name =  base_cluster_file[:base_cluster_file.index(".")] +'_tax_homogeneity.csv'

    # cluster_homogeneity_file = path.join(output_dir, file_name)
    # cluster_homogeneity_file =

    # print(file_name)
    print(inflation)
    print(cluster_homogeneity_file)

    unclassified_term_list = getListFromFile(unclassified_term_file)
    node_dict, count_leaves, leaves_taxonomy, genome_id_to_name = constructTaxonTree(taxonomy_file, unclassified_term_list)

    iter_clusters = clusterFileParser(cluster_file)
    alternative_taxon_id =getAlternativeTaxonId(alternative_taxon_id_file)

    csv_writer = initiateOuputFile(cluster_homogeneity_file)
    # print(alternative_taxon_id)
    if visualisation:
        for i, cluster_info in enumerate(iter_clusters):
            #cat data/taxonomy/taxonomy_virus.txt | cut -f3 | cut -d';' -f2 | sort | uniq | more
            dict_info = computeHomogeneity(cluster_info, unclassified_term_list, leaves_taxonomy, alternative_taxon_id, count_leaves)

            dict_info["inflation"] = inflation
            dict_info["coverage"] = coverage
            cluster_info.update(dict_info)
            # if len(cluster_info['genome_ids']) < 80:
            #     print("HOMOGENEITY", dict_info["homogeneity"])
            view.displayClusterTree(node_dict, cluster_info, count_leaves, genome_id_to_name, leaves_taxonomy)
    else:
        writeHomogeneity(iter_clusters, inflation, coverage, unclassified_term_list, leaves_taxonomy, alternative_taxon_id, count_leaves, genome_id_to_name, node_dict)

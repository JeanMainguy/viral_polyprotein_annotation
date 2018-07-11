# from ete3 import Tree
import logging, sys
from os import path
import re, csv, os
from time import clock
# import visualisation_taxonomic_tree_cluster as view
import numpy as np
# PythonDecorators/decorator_without_arguments.py
# class timer(object):
#
#     def __init__(self, f):
#         """
#         If there are no decorator arguments, the function
#         to be decorated is passed to the constructor.
#         """
#         self.f = f
#
#     def __call__(self, *args, **kwargs):
#         """
#         The __call__ method is not called until the
#         decorated function is called.
#         """
#         START_TIME = clock()
#         result = self.f(*args, **kwargs)
#         # print(self.f.__name__, 'run time:%6.2f'%(clock()-START_TIME), 'seconds.');
#         # print('-'*20)
#         return result
# @timer

def getAnnotatedProteinTaxId(stat_protein_file):
    annotated_genomes = {}
    with open(stat_protein_file) as csvfile:
        reader = csv.DictReader(csvfile, delimiter='\t')
        for row in reader:
            # print(row)
            if row['polyprotein_outline'] == "True":
                annotated_genomes.setdefault(row['taxon_id'], []).append(row['protein_id'])
                # print(row['taxon_id'], row['protein_id'])
    return annotated_genomes


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
            unclassified_node_name_to_merge = []
            taxon_to_remove = []

            for i, taxon in enumerate(taxonomy):
                if any([taxon.startswith(unclassified_term) for unclassified_term in unclassified_term_list]):
                    taxon_to_remove.append(taxon)
                    unclassified_node_name_to_merge.append(taxon)
                    continue
                elif unclassified_node_name_to_merge:
                    f = True
                    taxonomy[i] = '{} ({})'.format(taxon, ';'.join(unclassified_node_name_to_merge) )
                    unclassified_node_name_to_merge = []
            if unclassified_node_name_to_merge:
                genome_id_to_name[genome_id] = "{} ({})".format(genome_name, ';'.join(unclassified_node_name_to_merge))
            [taxonomy.remove(t) for t in taxon_to_remove]

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


def clusterFileParser(cluster_file, annotated_genomes):
    # 1093958|YP_004901696.1|362.0|None	759390|YP_003620408.1|362.0|None	2021667|YP_009408633.1|357.0|None
    with open(cluster_file, 'r') as cluster_handle:
        for i, cluster_line in enumerate(cluster_handle):
            cluster_elements = cluster_line.rstrip().split("\t")
            cluster_info = {'genome_ids':[], 'nb_polyprotein':0, 'nb_protein':len(cluster_elements), "cluster_id":i}
            cluster_info['polyproteins'] = []

            for element in cluster_elements:
                # print(element)
                #temporary solution should be change at some point
                try:
                    genome_id, protein_id, length, peptide_annotation = element.split('|')
                except ValueError: # Peptide info is not anymore in fasta header. use the dict annotated_genomes to identify annotated protein
                    genome_id, protein_id, length = element.split('|')
                if genome_id in annotated_genomes and protein_id in annotated_genomes[genome_id]:
                    peptide_annotation = "Peptide"
                else:
                    peptide_annotation = "None"
                if peptide_annotation == 'Peptide':
                    cluster_info['nb_polyprotein'] += 1
                    cluster_info['polyproteins'].append(genome_id)
                cluster_info['genome_ids'].append(element.split('|')[0])
            # genome_ids = {element.split('|')[0] for element in cluster}
            yield cluster_info


def getSharedTaxonomy(taxonomy_set):

    # valid_branches = set()
    taxonomy_iter = iter(taxonomy_set)
    current_shared_taxonomy = next(taxonomy_iter)

    for new_taxonomy in taxonomy_iter:
        if  current_shared_taxonomy[-1] in new_taxonomy:
            if new_taxonomy.index(current_shared_taxonomy[-1]) != current_shared_taxonomy.index(current_shared_taxonomy[-1]):
                logging.warning("The same node name found with different index in taxonomies:"+str(current_shared_taxonomy)+'_and_'+str(new_taxonomy))
            else:
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

def getValidBranches(last_common_node, current_shared_taxonomy, taxonomy_set):

    valid_branches = set()
    for tax in taxonomy_set:
        if last_common_node in tax:
            i_last = tax.index(last_common_node)

            if len(tax) > i_last+1:
                valid_branches.add(tax[i_last+1])



    if not valid_branches: # there is no better defined node in the cluster than the last_common_node
        return {last_common_node} # last common node is the node use to compute homogeneity

    if len(valid_branches) > 1:
        return valid_branches

    if len(valid_branches) == 1:
        # only one branch is in valid_branches which mean some genome are attached to last_common_node and all the other
        # ones have the better defined node of valid branch. which means we find potential valid branches in the valid ranch node
        last_common_node_from_valid_branch = valid_branches.pop()

        return getValidBranches(last_common_node_from_valid_branch, current_shared_taxonomy, taxonomy_set)


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

                logging.warning("taxon id not found in taxonomy file or in alternative taxon id file: "+taxon_id)
                raise NameError("taxon id not found in taxonomy file or in alternative taxon id file: ",taxon_id)
                # unfound_genome += 1
                continue
        # Check if the first node after the root is an unclassified node
        # if it is this genome is not take into account in the homogeneity
        # print(taxon_id, taxonomy)
        # for unclassified_term in unclassified_term_list:
        #     if taxonomy[1].startswith(unclassified_term):
        #         unclassified_genome += 1
        #         unclassified_taxonomy_set.add(taxonomy)
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
                    "evalue",
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
    valid_branches = getValidBranches(current_shared_taxonomy[-1], current_shared_taxonomy, effective_taxonomy_set)
    # input()
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
    # singleton = False if
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

# @timer
def writeHomogeneity(iter_clusters, inflation,evalue, coverage, unclassified_term_list, leaves_taxonomy, alternative_taxon_id, count_leaves, genome_id_to_name, node_dict):
    homogeneities = []
    homogeneities_poly = []
    homogeneities_no_poly = []
    homogeneities_dsDNA = []
    homogeneities_no_dsDNA = []


    for i, cluster_info in enumerate(iter_clusters):
        #cat data/taxonomy/taxonomy_virus.txt | cut -f3 | cut -d';' -f2 | sort | uniq | more
        dict_info = computeHomogeneity(cluster_info, unclassified_term_list, leaves_taxonomy, alternative_taxon_id, count_leaves)

        dict_info["inflation"] = inflation
        dict_info["coverage"] = coverage
        dict_info["evalue"] = evalue

        if csv_writer:
            csv_writer.writerow(dict_info)
            # continue
        if not dict_info['unclassified_cluster']:
            homogeneities.append(dict_info['homogeneity'])

            if dict_info["nb_polyprotein"] == 0:
                homogeneities_no_poly.append(dict_info['homogeneity'])
            else:
                homogeneities_poly.append(dict_info['homogeneity'])

            if 'dsDNA' in dict_info["shared_taxonomy"] or 'dsDNA' in dict_info["valid_branches"]:
                homogeneities_dsDNA.append(dict_info['homogeneity'])
            else:
                homogeneities_no_dsDNA.append(dict_info['homogeneity'])

    assert (len(homogeneities_dsDNA) + len(homogeneities_no_dsDNA)== len(homogeneities) ), 'homogeneity array are no symetric'

    print(build_summary_line(homogeneities, 'all', inflation, evalue, coverage))
    print(build_summary_line(homogeneities_no_poly, 'no poly', inflation, evalue, coverage))
    print(build_summary_line(homogeneities_poly, 'poly', inflation, evalue, coverage))
    print(build_summary_line(homogeneities_dsDNA, 'dsDNA', inflation, evalue, coverage))
    print(build_summary_line(homogeneities_no_dsDNA,'no dsDNA', inflation, evalue, coverage))

    summary_fl.write(build_summary_line(homogeneities, 'all', inflation, evalue, coverage)+'\n')
    summary_fl.write(build_summary_line(homogeneities_no_poly, 'no poly', inflation, evalue, coverage)+'\n')
    summary_fl.write(build_summary_line(homogeneities_poly, 'poly', inflation, evalue, coverage)+'\n')
    summary_fl.write(build_summary_line(homogeneities_dsDNA, 'dsDNA', inflation, evalue, coverage)+'\n')
    summary_fl.write(build_summary_line(homogeneities_no_dsDNA,'no dsDNA', inflation, evalue, coverage)+'\n')


    return homogeneities

def getListFromFile(file):
    list_line = []
    with open(file, "r") as hd:
        for l in hd:
            list_line.append(l.rstrip())
    return list_line

def build_summary_line(homogeneities, type, inflation, evalue, coverage):
    try:
        median=np.median(homogeneities)
        Q1 = np.percentile(homogeneities, 25)
        Q3 = np.percentile(homogeneities, 75)
        ymax = max(homogeneities)
        ymin = min(homogeneities)
        mean=np.mean(homogeneities)
    except IndexError:
        median, Q1, Q3, ymax, ymin, mean = None, None, None, None, None, None
    length = len(homogeneities)
    line = [median,Q1, Q3, ymax,  ymin, length, mean, coverage,  evalue, inflation, type]
    line = [str(i) for i in line]
    line = ','.join(line)
    return line


if __name__ == '__main__':

    logging.basicConfig(filename='log/homogeneity_cluster.log', level=logging.WARNING)
    clusters = sys.argv[1]
    cluster_homogeneity_dir = sys.argv[2]
    taxonomy_file = sys.argv[3]  #"data/taxonomy/taxonomy_virus.txt"
    alternative_taxon_id_file = sys.argv[4] # "data/taxonomy/heterogeneous_taxon_id_taxonomy_virus.txt"
    unclassified_term_file = sys.argv[5] # "data/taxonomy/unclassified_terms.txt"
    stat_protein_file= sys.argv[6] #'results/stat_viral_protein/stat_proteins_Viruses.csv'
    summary_file =  sys.argv[7]

    #initiate summary file
    summary_fl = open(summary_file, 'w')
    header = ["median","Q1","Q3","max","min","length","mean","coverage","evalue","inflation","type"]
    summary_fl.write(','.join(header)+'\n')

    annotated_genomes = getAnnotatedProteinTaxId(stat_protein_file)
    unclassified_term_list = getListFromFile(unclassified_term_file)
    node_dict, count_leaves, leaves_taxonomy, genome_id_to_name = constructTaxonTree(taxonomy_file, unclassified_term_list)
    alternative_taxon_id =getAlternativeTaxonId(alternative_taxon_id_file)

    if os.path.isdir(clusters):
        print('clusters is a directory, so we compute homogeneity for each of the cluster file we found')
        cluster_files = (os.path.join(clusters, f) for f in os.listdir(clusters) if f.endswith('.out')) # in the directory only cluster need to be present
    else:
        cluster_files = [clusters]

    for cluster_file in cluster_files:
        # we retrieve the coverage value and the inflation from the file name
        re_result = re.search("evalue_(1e[-\d]+)coverage(\d+)_I([\d_]+).out", cluster_file)
        #if the regex doesn't work an expception will be raised
        if re_result:
            evalue=  re_result.group(1)
            coverage = re_result.group(2)
            inflation = re_result.group(3)
        else:
            raise NameError('parsing cluster file name failed %s' % cluster_file)

        # re_result = re.search("evalue_(1e[-\d]+)coverage", cluster_file)
        # evalue = None if not re_result else re_result.group(1)
        # if evalue == None:
        #     raise NameError('parsing cluster file name failed %s' % cluster_file)
        iter_clusters = clusterFileParser(cluster_file, annotated_genomes)

        if cluster_homogeneity_dir == 'None':
            csv_writer = False
        else:
            base_name = os.path.basename(cluster_file)
            cluster_homogeneity_file = os.path.join(cluster_homogeneity_dir,base_name.replace('.out', '_homogeneity.csv'))
            csv_writer = initiateOuputFile(cluster_homogeneity_file)


        writeHomogeneity(iter_clusters, inflation, evalue, coverage, unclassified_term_list, leaves_taxonomy, alternative_taxon_id, count_leaves, genome_id_to_name, node_dict)

    summary_fl.close()

#!/usr/bin/env python3

import taxonomy as tax
import viral_genome_classes as obj
import viruses_statistics as stat
import parser_interpro_results as do
import visualisation_alignment as view_aln
import visualisation_protein as view_prot
import viral_protein_extraction as extract

from os import path, listdir, rename
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

        # print(set_retrieved)
        # print(set_requested - set_retrieved)

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

        for site in cds.cleavage_sites:
            site.cds_of_aln.append(cds)
            site_start = site.start_aa(cds) - 1  # to be in base 0 we need the -1
            site_end = site.end_aa(cds) - 1
            site_start_in_aln = cds.aln_list.index(site_start)
            site_end_in_aln = cds.aln_list.index(site_end)
            site.start_in_aln = site_start_in_aln
            site.end_in_aln = site_end_in_aln

            findCloseSites(site, site_start_in_aln, site_end_in_aln, annotated_cds_list, window)

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
            # if not cds.polyprotein:
            #     print(cds)
            #     print(site)
            #     input()
            # -1 to be in base 0 because it'll be use as a list index
            site_start = site.start_aa(cds) - 1
            site_end = site.end_aa(cds) - 1
            protein_in_group_list += site.cds_of_aln
            # pfnt(len([i for i in cds.aln_list if i is not None]))
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
    # [ nb_of_cleavage_site, nb_protein_in_group, nb_cleavage_from_the_same_protein, standard_dev, round(standard_dev)]
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


def get_predicted_mat_peptide(cds, cleavage_site_attrb="predicted_cleavage_sites"):

    end_mat_pep = [cs.start_aa(cds) for cs in getattr(cds, cleavage_site_attrb)]
    end_mat_pep.sort()
    # end of the prot -3 (because of the stop codon) == end of the last mat_pep
    end_mat_pep.append(int(len(cds)/3 - 1))
    # print('end_mat_pep', end_mat_pep)
    # print(len(cds))
    # print(len(cds)/3)
    # input()
    start_position = 1  # cds.start  # start of the first pep
    for end_position in end_mat_pep:
        # Write the mat pep:
        # print(start_position, end_position)
        pep = obj.Predicted_peptide(cds, start_in_prot=start_position, end_in_prot=end_position)
        cds.predicted_mat_peptides.append(pep)

        start_position = end_position + 1


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


def visualisation_of_processed_aln(cds_list, alignment_file, group_info_list, file_handle, display_line_size=180):

    cds_annotated = [cds for cds in cds_list if cds.polyprotein]
    # VISUALISATION OF THE ALIGNMENT
    len_cds_max = max((len(cds) for cds in cds_list))
    # for cds in cds_annotated:
    # for cds in cds_list:
    # print(view_prot.visualisation_protein(cds, 1, len_cds_max))

    for cds in cds_list:
        cds.matchs = []

    view_aln.visualisation(cds_list, alignment_file, display_line_size,
                           group_info_list, file_handle)


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
            # remove valid needed to visualise but useless for the csv
            group_info.pop('valid', None)
            group_site_csv_writer.writerow(group_info)


def write_predicted_annotation(writer, cds, cluster_nb):
    info_dict = {"cluster_nb": cluster_nb,
                 "taxon_id": cds.segment.taxon_id,
                 "virus_name": cds.segment.organism,
                 "assembly_version": "GCF_..",
                 "protein_id": cds.protein_id}

    for predicted_cs in cds.predicted_cleavage_sites:
        info_dict["cs_P1_protein_position"] = predicted_cs.start_in_prot
        info_dict["confidence_score"] = predicted_cs.confidence_score
        writer.writerow(info_dict)


def write_gff_annotation(writer, cds):
    # attribute = f'GCF..;taxon_id_{cds.segment.taxon_id};cluster_{cluster_nb}'
    info_dict = {"attributes": f"Dbxref=taxon:{cds.segment.taxon_id}",
                 "source": '.',
                 "score": '.',
                 "strand": '.',
                 "phase": '.'}
    # Write the cds:
    info_dict["seqid"] = cds.protein_id
    info_dict["type"] = "polypeptide"
    info_dict["start"] = 1
    info_dict["end"] = int(len(cds)/3) - 1  # doesn't take into account stop codon
    writer.writerow(info_dict)

    # start of the cs is the first aa of the cs
    # which is the last aa of mat_peptides
    # +1 because gff starts at 1 and not 0
    cs_positions = [cs.start_in_prot + 1 for cs in cds.predicted_cleavage_sites]
    cs_positions.append(int(len(cds)/3)-1)  # end of the prot == end of the last mat_pep
    start_position = 1  # start of the first pep
    for end_position in cs_positions:
        # Write the mat pep:
        info_dict["seqid"] = cds.protein_id
        info_dict["type"] = "mat_peptide"
        info_dict["start"] = start_position
        info_dict["end"] = end_position
        writer.writerow(info_dict)
        start_position = end_position + 1


def write_reannotated_genbank_file(cds, output_dir):

    # gb_file = cds.segment.gb_file
    # # for cs in cds.
    # Predicted_peptide(cds, start, end)
    reannotated_gb_file = path.join(output_dir, f'genome_{cds.segment.taxon_id}.gbff')
    tmp_file = path.join(output_dir, 'tmp.gbff')
    if path.exists(reannotated_gb_file):
        fl_read = open(reannotated_gb_file, 'r')
    else:
        gb_file = cds.segment.gb_file
        proper_open = gzip.open if gb_file.endswith('.gz') else open
        fl_read = proper_open(gb_file, 'rt')

    fl_write = open(tmp_file, 'w')
    # print(reannotated_gb_file)
    for record in SeqIO.parse(fl_read, "genbank"):
        if record.id == cds.segment.record.id:
            cds_index = [i for i, f in enumerate(record.features) if cds.bp_obj.qualifiers ==
                         f.qualifiers and cds.bp_obj.type == f.type and cds.bp_obj.location == f.location].pop()
            predicted_peps_bp_obj = [pep.bp_obj for pep in cds.predicted_mat_peptides]
            record.features = record.features[:cds_index+1] + \
                predicted_peps_bp_obj + record.features[cds_index+1:]

        SeqIO.write(record, fl_write, 'genbank')

    rename(tmp_file, reannotated_gb_file)  # os import

    fl_read.close()
    fl_write.close()


def initiate_gff_file(genome_file_name):

    genome_header = ["seqid",  # name of the chromosome or scaffold; chromosome names can be given with or without the 'chr' prefix. Important note: the seq ID must be one used within Ensembl, i.e. a standard chromosome name or an Ensembl identifier such as a scaffold ID, without any additional content such as species or assembly. See the example GFF output below.
                     # name of the program that generated this feature, or the data source (database or project name)
                     "source",
                     "type",  # type of feature. Must be a term or accession from the SOFA sequence ontology
                     # Start position of the feature, with sequence numbering starting at 1.
                     "start",
                     "end",  # End position of the feature, with sequence numbering starting at 1.
                     "score",  # A floating point value.
                     "strand",  # defined as + (forward) or", # (reverse).
                     "phase",  # One of '0', '1' or '2'. '0' indicates that the first base of the feature is the first base of a codon, '1' that the second base is the first base of a codon, and so on..
                     "attributes",  # A semicolon#separated list of tag#value pairs, providing additional information about each feature. Some of these tags are predefined, e.g. ID, Name, Alias, Parent", # see the GFF documentation for more details.
                     ]

    handle_genome_out = open(genome_file_name, "w")
    handle_genome_out.write("##gff-version 3\n")
    genome_csv_writer = csv.DictWriter(handle_genome_out, fieldnames=genome_header, delimiter='\t')
    # genome_csv_writer.writeheader()

    return genome_csv_writer, handle_genome_out


def initiate_reannotated_genome_file(genome_file_name):

    genome_header = ["cluster_nb",
                     "taxon_id",
                     "virus_name",
                     "assembly_version",
                     "protein_id",
                     "cs_P1_protein_position",  # P1 is aa juste before the cleavage
                     "cs_P1_genome_position",
                     "confidence_score"]

    handle_genome_out = open(genome_file_name, "w")
    genome_csv_writer = csv.DictWriter(handle_genome_out, fieldnames=genome_header, delimiter='\t')
    genome_csv_writer.writeheader()

    return genome_csv_writer, handle_genome_out


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
        # print(taxon_prot_ids)
        cds_list = getCdsObject(taxon_prot_ids, taxonomy_file, interpro_domains_file, sp_treshold)
        for cds in cds_list:
            if f'{cds.segment.taxon_id}|{cds.protein_id}' in black_list:
                # cds annotation has been black listed
                cds.black_listed_peptides = cds.peptides
                cds.black_listed_cleavage_sites = cds.cleavage_sites
                cds.peptides = set()
                cds.cleavage_sites = []
                cds.polyprotein = None

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
                group_info['valid'] = True
            else:
                group_info['valid'] = False

        iter_cds_non_annotated = (cds for cds in cds_list if not cds.polyprotein)
        for blank_cds in iter_cds_non_annotated:
            get_predicted_mat_peptide(blank_cds)

        ######################
        #   OUTPUT
        ######################

        if not aln_csv_writer or 1:
            # print(f'VISUALISATION WITH WINDOW {window}')

            # for cds in cds_list:
            #     print(cds.segment.record.annotations['taxonomy'])

            # for group in group_info_list:
            #     print(group)
            if write_visualization_path:
                file_handle = open(path.join(write_visualization_path,
                                             f"visualization_cluster{cluster_nb}.aln"), 'w')
            else:
                file_handle = sys.stdout
            visualisation_of_processed_aln(cds_list, alignment_file,
                                           group_info_list, file_handle, display_line_size)

        if aln_csv_writer:
            parameters["cluster_nb"] = cluster_nb
            parameters['window'] = window
            compute_alignment_stat(group_info_list, parameters, cds_list,
                                   aln_csv_writer, group_site_csv_writer)

        if reannotated_genome_dir:
            iter_cds_non_annotated = (cds for cds in cds_list if not cds.polyprotein)

            for blank_cds in iter_cds_non_annotated:
                taxon_id = blank_cds.segment.taxon_id
                if taxon_id not in genomes_output_writers:
                    genomes_output_writers[taxon_id], fl_handle = initiate_reannotated_genome_file(
                        path.join(reannotated_genome_dir, f'cleavage_sites_{taxon_id}.csv'))
                    file_handles.append(fl_handle)  # to close all files properly at the end

                    gff_writers[taxon_id], fl_handle = initiate_gff_file(
                        path.join(reannotated_genome_dir, f'mature_peptides_{taxon_id}.gff'))
                    file_handles.append(fl_handle)  # to close all files properly at the end

                genome_csv_writer = genomes_output_writers[taxon_id]
                write_predicted_annotation(genome_csv_writer, blank_cds, cluster_nb)

                gff_writer = gff_writers[taxon_id]
                write_gff_annotation(gff_writer, blank_cds)

                write_reannotated_genbank_file(blank_cds, reannotated_genome_dir)


if __name__ == '__main__':

    # logging.basicConfig(filename='log/alignment_analysis.log', level=logging.INFO)
    # alignment_file= "data/alignment/Retro-transcribing_viruses_1e-5_coverage90_I3/seq_cluster16.aln"

    try:
        # 'data/alignment/Viruses_1e-5_coverage90_I2/seq_cluster1037.aln'
        alignment_file_or_dir = sys.argv[1]
    except:
        alignment_file_or_dir = 'results_db_viral_2018-10-19/Picornavirales_evalue_1e-60coverage60_I1_8/alignment/seq_cluster6.aln'
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
        taxonomy_file = "results_db_viral_2018-10-19/genomes_index/taxonomy_virus.txt"

    try:
        reannotated_genome_dir = sys.argv[6]
    except:
        reannotated_genome_dir = "test/"
    try:
        interpro_domains_file = sys.argv[7]
    except:
        interpro_domains_file = None

    black_list_file = "results_db_viral_2018-10-19/viral_protein_stat/black_list_annotation_Picornavirales.txt"

    sp_treshold = 90
    display_line_size = 75  # 140/2
    confidence_score_threshold = 4
    if output_stat_aln:
        write_visualization_path = path.dirname(output_stat_aln)
    else:
        write_visualization_path = False

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

    black_list = []
    # Black list preparation
    if black_list_file:
        with open(black_list_file, "r") as fl:
            black_list = fl.read().split('\n')

    genomes_output_writers = {}
    gff_writers = {}

    for alignment_file in alignment_files:
        print(alignment_file)
        main()  # variable needed in main() are global

    for fl in file_handles:
        fl.close()

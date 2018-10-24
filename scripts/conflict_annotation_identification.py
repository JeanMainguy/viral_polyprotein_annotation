#!/usr/bin/env python3
import csv
import sys
from os import path


def conflicting_overlapping_events_identification(domain_dict, dict_seq,
                                                  threshold_overlap_prct,
                                                  threshold_overlap_aa,
                                                  ignoring_threshold_ratio,
                                                  notify_file=None):
    """
    Find domain annotations that seem to always overlap (use of ignore_threshold).
    In order to ignored them in the rest of the analysis (ignore)
    All the other domain annotations that are found overlapping more than the acceptable_overlapping_threshold are consider
    as conflicting domain annotation and are written in a file
    Return list of taxon_id|protein_id where conflict have been identified, that can be blacklisted
    INPUT:
        domain_dict: dictionnary with as key: domain name and as value a dictionnary
                    containing the folowing content:
                    'overlapping_dist':[],
                    'overlapping_dist_prct':[],
                    'not_overlapping_count':0}

        acceptable_overlap_threshold : threshold of the overlapping length domain
                                           (in aa or in percentage) which is consider
                                           as accepatble
        ignoring_threshold_ratio : percentage of time
                                  a domain annotation has
                                  to be found overlapping to be
                                  ignored by the program
    """
    if not threshold_overlap_aa and not threshold_overlap_prct:
        raise ValueError  # 'accepetable overlap thresholds are not provided'

    threshold_overlap_aa_param = threshold_overlap_aa
    domain_total_count = 0
    conflict_cds_count = 0
    ignored_domains = {}
    conflicting_domains = {}
    safe_domains = {}
    black_listed_cds = []
    if notify_file:
        intro_param = ''
        fl_handle = open(notify_file, 'w')
        intro_param += 'Identification of conflict between domain annotations and cleavage site annotations:\n'
        intro_param += 'Parameters:\n'
        intro_param += "Maximum overlapping distance (in aa or in percentage of the domain length) which is consider as accepatble:\n"
        intro_param += f"acceptable_overlap_threshold = {threshold_overlap_aa}aa "
        intro_param += f"or {threshold_overlap_prct}% of the domain length\n"

        intro_param += f'A domain annotation need to be found overlapping '
        intro_param += f'in {ignoring_threshold_ratio*100}% of cases to be ignored by the program\n\n'
        intro_param += "The distance of overlapping corresponds to the smallest distance of the domain "
        intro_param += "on the right or on the left of the cleavage site "
        fl_handle.write(intro_param)

    for domain, info in domain_dict.items():
        domain_total_count += info['total']
        # SET UP nb of acceptable overlapping event
        nb_acceptable_overlap = 0
        nb_overlap = 0
        overlapping_info = []
        # print(info)
        # input()
        for dist, cds, len_do in zip(info['overlapping_dist'], info['overlapped_cds'], info['len_do']):
            # transformation of parameters to get only one value to test:
            if not threshold_overlap_aa_param:
                threshold_overlap_aa = (threshold_overlap_prct/100)*len_do
            if threshold_overlap_aa_param and threshold_overlap_prct:
                threshold_overlap_aa = min(
                    threshold_overlap_aa_param, (threshold_overlap_prct/100)*len_do)

            if int(dist) < threshold_overlap_aa:
                nb_acceptable_overlap += 1

            else:
                overlapping_info.append((cds, dist, len_do))
                nb_overlap += 1

        overlapping_time_ratio = nb_overlap/info['total']
        if overlapping_time_ratio >= ignoring_threshold_ratio:
            print('//'*20)
            print('ignored DOMAIN')
            print('acceptable overlapping', nb_acceptable_overlap)
            print('OVERLAPPING', nb_overlap)
            print(domain)
            [print(k, v) for k, v in info.items()]
            ignored_domains[domain] = info

        elif nb_overlap != 0:
            print('*!!*'*20)
            print('conflicting DOMAIN')
            print("overlapping_time_ratio", overlapping_time_ratio)
            print('acceptable overlapping', nb_acceptable_overlap)
            print('OVERLAPPING', nb_overlap)
            print(domain)
            [print(k, v) for k, v in info.items()]
            string = '\n' + '--'*50 + '\n'
            string += f'Conflict with the domain annotation {domain}\n'
            string += f'The domain annotation is found in {info["total"]} sequences\n'
            string += f'It is overlapping cleavage sites in {info["not_overlapping_count"]} sequences\n'
            string += f'In {nb_overlap} sequences the overlapping distance is greater than the given treshold\n'
            for cds, dist, len_do in overlapping_info:
                string += f'   Sequence {cds} : overlapping distance {dist}aa ({dist/len_do*100:.1f}% of the domain lenght)\n'
                black_listed_cds.append(cds)
                conflict_cds_count += 1
            conflicting_domains[domain] = info
            if notify_file:
                fl_handle.write(string)
            else:
                print(string)
        else:
            safe_domains[domain] = info

    print('ignored', len(ignored_domains))
    print('safe_domain', len(safe_domains))
    print('conflicting domain', len(conflicting_domains))

    if notify_file:
        summary = '\n' + '=='*50 + '\n'
        summary += f'{len(domain_dict)} different domain annotations have been '
        summary += f'found {domain_total_count} times in the {len(dict_seq)} analysed sequences\n'
        summary += f'   Domain annotations that do not overlap cleavage sites: {len(safe_domains)}\n'
        summary += '   Domain annotations overlapping but ignored because they overlap '
        summary += f'in more than {ignoring_threshold_ratio*100}% of the cases: {len(ignored_domains)}\n'
        summary += f'   Domain annotations in conflict with cleavage sites; {len(conflicting_domains)}\n'
        summary += f"   Sequence annotation in conflict with domain annotation: {conflict_cds_count}"
        fl_handle.write(summary)
        fl_handle.close()
    return black_listed_cds


def get_domain_overlapping_dict(domain_stat_file):
    """
    From the stat file of domain, compute a dictionary with overlapping info for each domain annotation.
    """
    domain_dict = {}
    seq_dict = {}
    set_seq = set()
    with open(domain_stat_file) as fl:

        reader = csv.DictReader(fl, delimiter='\t')

        for row in reader:
            cds = f"{row['taxon_id']} | {row['protein_id']}"
            if row['name'] not in domain_dict:
                domain_dict[row['name']] = {'overlapping_dist': [],
                                            'overlapping_dist_prct': [],
                                            "overlapped_cds": [],
                                            "len_do": [],
                                            'not_overlapping_count': 0,
                                            "total": 0}
            if cds not in seq_dict:
                seq_dict[cds] = {'overlapping_domains': {},
                                 'overlapping_domain_count': 0,
                                 'not_overlapping_count': 0,
                                 "total": 0}

            if row['overlaps_cleavage_site'] == 'True':
                domain_dict[row['name']]["overlapping_dist"].append(
                    int(row['overlapping_distance']))
                domain_dict[row['name']]["overlapping_dist_prct"].append(
                    float(row['overlapping_size_percentage']))
                domain_dict[row['name']]["overlapped_cds"].append(cds)
                len_do = int(row['end_in_prot']) - int(row['start_in_prot']) + 1
                domain_dict[row['name']]['len_do'].append(len_do)

                if row['name'] not in seq_dict[cds]["overlapping_domains"]:
                    seq_dict[cds]["overlapping_domains"][row['name']] = [domain_dict]
                else:
                    seq_dict[cds]["overlapping_domains"][row['name']].append(domain_dict)
                seq_dict[cds]['overlapping_domain_count'] += 1
            else:
                domain_dict[row['name']]["not_overlapping_count"] += 1
                seq_dict[cds]["not_overlapping_count"] += 1

            domain_dict[row['name']]['total'] += 1
            seq_dict[cds]["total"] += 1
    return domain_dict, seq_dict


if __name__ == '__main__':

    domain_stat_file = sys.argv[1]

    black_liste_cds_file = sys.argv[2]

    threshold_overlap_prct = sys.argv[3]  # 10
    threshold_overlap_aa = sys.argv[4]  # 10
    ignoring_threshold_ratio = float(sys.argv[5])  # 0.5

    notify_file = sys.argv[6]

    try:
        threshold_overlap_prct = float(threshold_overlap_prct)
    except ValueError:
        threshold_overlap_prct = None

    try:
        threshold_overlap_aa = int(threshold_overlap_aa)
    except ValueError:
        threshold_overlap_aa = None

    dict_domain, dict_seq = get_domain_overlapping_dict(domain_stat_file)
    black_listed_cds = conflicting_overlapping_events_identification(dict_domain, dict_seq,
                                                                     threshold_overlap_prct,  # 5% of the domain length
                                                                     threshold_overlap_aa,  # 5aa
                                                                     # if more than 0.5 percent of the domain annotation is overlapping.. then it is blacklister
                                                                     ignoring_threshold_ratio,
                                                                     notify_file)
    with open(black_liste_cds_file, 'w') as fl:
        fl.write('\n'.join(black_listed_cds))

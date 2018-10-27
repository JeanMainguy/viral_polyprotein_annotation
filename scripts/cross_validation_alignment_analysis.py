import multiple_alignment_analysis as analysis
from os import path, listdir
import csv
import re
import logging
from numpy import std, mean
import sys


def get_positive_negative_distribution(annotated_cds, cds_unannotated_version, window, positive, negative):

    for predicted_site in cds_unannotated_version.predicted_cleavage_sites:
        state = "False_positive"
        for annotation_site in annotated_cds.cleavage_sites:
            if predicted_site.start_aa(cds_unannotated_version)-window/2 <= annotation_site.start_aa(annotated_cds) <= predicted_site.start_aa(cds_unannotated_version)+window/2:
                # print('prediction match annotation')
                state = "True_positive"
        if state == "True_positive":
            positive.append(predicted_site.confidence_score)
        else:
            negative.append(predicted_site.confidence_score)


def annotation_vs_prediction(annotated_cds, cds_unannotated_version, window, score_cutoff):
    state_summary = {"False_positive": 0, "True_positive": 0, "False_negative": 0}
    # FALSE POSITIVE DETECTION. predicted cleavage_site not found in the original annotation
    prediction_vs_annotation = []
    predicted_cleavage_site = [
        site for site in cds_unannotated_version.predicted_cleavage_sites if site.confidence_score < score_cutoff]
    for predicted_site in predicted_cleavage_site:
        state = "False_positive"
        for annotation_site in annotated_cds.cleavage_sites:
            if predicted_site.start_aa(cds_unannotated_version)-window/2 <= annotation_site.start_aa(annotated_cds) <= predicted_site.start_aa(cds_unannotated_version)+window/2:
                # print('prediction match annotation')
                state = "True_positive"
        prediction_vs_annotation.append(state)
    # False NEGATIVE DETECTION:When there is no prediction corresponding to an annotated cleavage site
    annotation_vs_prediction = []
    for annotation_site in annotated_cds.cleavage_sites:
        state = "False_negative"
        for predicted_site in predicted_cleavage_site:
            if annotation_site.start_aa(annotated_cds)-window/2 <= predicted_site.start_aa(cds_unannotated_version) <= annotation_site.start_aa(annotated_cds)+window/2:
                state = "True_positive"
        annotation_vs_prediction.append(state)

    # print("prediction_vs_annotation", prediction_vs_annotation)
    # print('annotation_vs_prediction', annotation_vs_prediction)
    state_summary["False_positive"] += prediction_vs_annotation.count("False_positive")
    state_summary["True_positive"] += prediction_vs_annotation.count("True_positive")
    state_summary["False_negative"] += annotation_vs_prediction.count("False_negative")

    # assert prediction_vs_annotation.count('True_positive') == annotation_vs_prediction.count('True_positive'), 'Nb of True Positive is different in prediction vs annotation and annotation vs prediction'
    # True positive from annotation can be more abundant because if more than one annotated cleavage site  from a cds
    # is found in the same group (very close cleavage site) then each of them would be count as True positive because they would be close to a
    # predicted one
    # then we take only into account True Positive computed by the comparison of prediction compare to annotation

    # print(state_summary)

    relevant_annotation = (state_summary["False_negative"]+state_summary["True_positive"])
    positive_element = (state_summary["False_positive"]+state_summary["True_positive"])

    try:
        precision = state_summary["True_positive"] / positive_element
    except ZeroDivisionError:
        precision = 0

    try:
        recall = state_summary["True_positive"] / relevant_annotation
    except ZeroDivisionError:
        recall = 0

    # print('recall', recall)
    # print('precision', precision)
    return precision, recall


def compute_cross_validation_stat(group_info_list, parameters, cds_list, aln_csv_writer, precision, recall):

    cds_annotated = [cds for cds in cds_list if cds.polyprotein]
    cleavages_sites_per_seq = [len(cds.cleavage_sites) for cds in cds_annotated]

    general_info = {'nb_sequences': len(cds_list),
                    'nb_seq_with_annotation': len(cds_annotated),
                    'cleavages_sites_per_seq': cleavages_sites_per_seq}
    general_info.update(parameters)

    try:
        # Need to float the value because otherwise the denominator is of type numpy float and don't trow any
        # exception but an obscur warning when the float == 0.0
        F1_score = 2*(float(precision*recall)/float((precision + recall)))
    except ZeroDivisionError:
        F1_score = 0.0

    # print('number of group', len(group_info_list))
    nb_group_with_valid_score = sum(
        [1 for grp in group_info_list if grp['confidence_score'] < parameters["confidence_score_threshold"]])
    # print( f'nb group that have a completeness > of {completeness_threshold} is {nb_group_with_valid_score}')
    # print(f'cleavage site per seq: {cleavages_sites_per_seq}')

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
                     "nb_group_with_valid_score": nb_group_with_valid_score,
                     "precision": precision,
                     "recall": recall,
                     "F-score": F1_score
                     }
    aln_dict_stat.update(general_info)

    if aln_csv_writer:

        aln_csv_writer.writerow(aln_dict_stat)

    # if group_site_csv_writer:
    #     for group_info in group_info_list:
    #         group_info.update(general_info)
    #         group_info['aln_validity'] = aln_validity
    #
    #         group_info['cleavages_sites_per_seq'] = None
    #         print(f'group {group_info["group_index"]} score: {group_info["confidence_score"]}')
    #         group_site_csv_writer.writerow(group_info)
    #


def initiate_ouput(output_stat_aln):

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
                  "precision",
                  "recall",
                  "F-score",
                  "cleavages_sites_per_seq"]

    handle_aln_out = open(output_stat_aln, "w")
    aln_csv_writer = csv.DictWriter(handle_aln_out, fieldnames=aln_header, delimiter='\t')
    aln_csv_writer.writeheader()
    file_handles = [handle_aln_out]
    return aln_csv_writer, file_handles


def get_cds_lists_with_artificial_unannotated_cds(cds_list, annotated_cds_list):
    cds_list_original = list(cds_list)

    if len(annotated_cds_list) == 1:  # cross validation is not possible
        return []
    for annotated_cds in annotated_cds_list:
        # print('annotated_cds')
        # print(annotated_cds.protein_id)
        cds_list = list(cds_list_original)
        cds_unannotated_version = annotated_cds.get_unnannotated_version()
        cds_unannotated_version.info_to_display = 'validation'
        # Give attribute relative to the aln to the cds_unannotated_version:
        cds_unannotated_version.aligned_sequence = annotated_cds.aligned_sequence
        cds_unannotated_version.aln_list = annotated_cds.aln_list

        cds_list.remove(annotated_cds)
        cds_list.append(cds_unannotated_version)
        # cds_list[cds_list_original.index(annotated_cds)] = cds_unannotated_version
        # print('original:', [(cds.protein_id, cds.polyprotein) for cds in cds_list_original])
        # print("changed :", [(cds.protein_id, cds.polyprotein) for cds in cds_list])
        # print()
        # print("annotated_cds", annotated_cds.protein_id)
        # print("cds_unannotated_version",cds_unannotated_version.number)

        yield annotated_cds, cds_unannotated_version, cds_list


def main():
    re_result = re.search("cluster(\d+).aln", alignment_file)
    cluster_nb = "NA" if not re_result else re_result.group(1)
    print("PROCESS of CLUSTER ", cluster_nb)
    taxon_prot_ids, seq_aln_dict = analysis.parse_alignment_file(alignment_file)

    for window in windows:
        cds_list = analysis.getCdsObject(taxon_prot_ids, taxonomy_file, gff_file, sp_treshold)

        if len(cds_list) != len(seq_aln_dict):
            raise Error('Not all cds from the alignment have been retrieved\n')

        for cds in cds_list:
            cds.aligned_sequence = seq_aln_dict[cds.protein_id]
            cds.aln_list = analysis.convertAlignmentToList(cds.aligned_sequence)
            cds.cleavage_site_positions = {s.start_aa(cds): s for s in cds.cleavage_sites}

        # extract.extract_cleavage_site_sequences(cds_list, 4)
        # original propagation using all the annotated seqeunce
        # use to write in the csv file the general info regarding the cluster
        cds_list_original = list(cds_list)
        group_info_list, cs_group_list = analysis.analyse_cleavage_site_groups(
            cds_list_original, window)
        for group_info in group_info_list:
            group_of_cs = cs_group_list[group_info["group_index"]]
            if group_info['confidence_score'] < confidence_score_threshold:
                analysis.propagate_cleavage_sites(
                    group_of_cs, group_info, cds_list_original, window)

        annotated_cds_list = [cds for cds in cds_list_original if cds.polyprotein]
        nb_cds_annotated = len(annotated_cds_list)
        if nb_cds_annotated <= 1:
            logging.info(
                'Only one annotated cds in the cluster, we cannot compute the cross validation')
            return

        # CROSS VALIDATION
        precisions = []
        recalls = []
        cds_modify_cds_list_iter = get_cds_lists_with_artificial_unannotated_cds(
            cds_list_original, annotated_cds_list)

        i = 0
        for annotated_cds, cds_unannotated_version, cds_list in cds_modify_cds_list_iter:
            i += 1
            print(f'annotated seq {i}/{nb_cds_annotated}')
            group_info_list_cross_val, cs_group_list = analysis.analyse_cleavage_site_groups(
                cds_list, window)

            for group_info in group_info_list_cross_val:
                group_of_cs = cs_group_list[group_info["group_index"]]
                # if group_info['confidence_score'] < confidence_score_threshold:
                analysis.propagate_cleavage_sites(group_of_cs, group_info, cds_list, window)

            # print(f'VISUALISATION WITH WINDOW {window}')
            # cds_list = sorted(cds_list, key=lambda x: x.protein_id)
            # visualisation_of_processed_aln(cds_list, alignment_file, group_info_list, display_line_size)
            if aln_csv_writer:
                precision, recall = annotation_vs_prediction(
                    annotated_cds, cds_unannotated_version, window, confidence_score_threshold)
                precisions.append(precision)
                recalls.append(recall)

            get_positive_negative_distribution(
                annotated_cds, cds_unannotated_version, window, positive, negative)

    # print("positive", len(positive))
    # print('negative',len(negative) )
    #
    # print("positive=", positive)
    # print()
    # print('negative=',negative)

    if aln_csv_writer:
        mean_precision = mean(precisions)
        mean_recall = mean(recalls)
        parameters["cluster_nb"] = cluster_nb
        parameters['window'] = window
        # print(parameters)
        compute_cross_validation_stat(group_info_list, parameters, cds_list_original,
                                      aln_csv_writer, mean_precision, mean_recall)

    # print('mean_precision', mean_precision)
    # print('mean_recall', mean_recall)


if __name__ == '__main__':

    try:
        # 'data/alignment/Viruses_1e-5_coverage90_I2/seq_cluster1037.aln'
        alignment_file_or_dir = sys.argv[1]
    except IndexError:
        alignment_file_or_dir = 'data/alignment/Viruses/RefSeq_download_date_2018-07-21/Viruses_evalue_1e-40coverage40_I2'
        # alignment_file_or_dir = 'data/alignment/Viruses/RefSeq_download_date_2018-07-21/Viruses_evalue_1e-40coverage40_I2/seq_cluster22.aln'
    try:
        # 'data/alignment/Viruses_1e-5_coverage90_I2/seq_cluster1037_stat.csv'
        windows_input = sys.argv[2]
    except IndexError:
        windows_input = "30"

    try:
        output_stat_aln = sys.argv[3]
        total_positive_negative_file = sys.argv[4]
    except IndexError:
        output_stat_aln = "test_aln.csv"
        output_stat_aln = False
        total_positive_negative_file = 'test_total_positive_negative_file.txt'
    try:
        taxonomy_file = sys.argv[5]  # "data/taxonomy/taxonomy_virus.txt"
        gff_file = sys.argv[6]
    except IndexError:
        gff_file = 'data/interpro_results_OLD/interproscan-5.30-69.0/domains_viral_sequences.gff3'
        taxonomy_file = "data/taxonomy/taxonomy_virus.txt"

    sp_treshold = 90
    display_line_size = 180
    confidence_score_threshold = 5

    parameters = {"confidence_score_threshold": confidence_score_threshold}

    if path.isdir(alignment_file_or_dir):
        alignment_files = (path.join(alignment_file_or_dir, f)
                           for f in listdir(alignment_file_or_dir) if f.endswith('.aln'))
        if alignment_file_or_dir.endswith('splitted'):
            parameters['cluster_with_isoforms_splitted'] = True
        else:
            parameters['cluster_with_isoforms_splitted'] = False

        re_result = re.search(
            "_evalue_([\de-]+)coverage(\d+)_I(\d_{0,1}\d*)", alignment_file_or_dir)
        parameters['evalue'] = "NA" if not re_result else re_result.group(1)
        parameters['coverage'] = "NA" if not re_result else re_result.group(2)
        parameters['inflation'] = "NA" if not re_result else re_result.group(3)

    elif path.isfile(alignment_file_or_dir) and alignment_file_or_dir.endswith('.aln'):
        alignment_files = [alignment_file_or_dir]
    else:
        raise ValueError('file or directory provided is not correct', alignment_file_or_dir)

    print('INPUT')
    print(alignment_file_or_dir)
    print("OUTPUT:")
    print(output_stat_aln)
    # window includ the cleavage in the middle
    # cleavage sites have a length of 2
    # then a window will always be even and >= 2
    # consequently we add 1 to odd window
    windows = {int((int(w)+int(w) % 2)) for w in windows_input.split(' ')}  # 10,20,30
    assert min(windows) >= 0
    windows = list(windows)
    windows.sort()
    print('window used to analyse cleavage sites', windows)

    if output_stat_aln:
        print("INITIATE OUTPUT FILES")
        aln_csv_writer, file_handles = initiate_ouput(output_stat_aln)
    else:
        aln_csv_writer = False
        file_handles = []

    positive = []
    negative = []

    for alignment_file in alignment_files:

        main()

    with open(total_positive_negative_file, "w") as fl:

        fl.write(f"positive={','.join([str(i) for i in positive])}\n")
        fl.write(f"negative={','.join([str(i) for i in negative])}\n")

    for fl in file_handles:
        fl.close()

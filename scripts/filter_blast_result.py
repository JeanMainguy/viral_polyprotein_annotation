import sys, os


def filter_by_coverage(blast_result, out_files):
    accepeted_pair = 0
    rejected_pair = 0
    length_rejected = 0

    with open(blast_result, 'r') as result_reader:
        for i, l in enumerate(result_reader):
            qseqid, sseqid, qcovs, qcovhsp, evalue, bitscore = l.split('\t')

            len_query = float(qseqid.split('|')[2])
            len_subject = float(sseqid.split('|')[2])

            if len_query < len_subject:
                length_rejected += 1
                continue
            accepeted =False
            for threshold, file in out_files.items():
                if int(qcovs) >= threshold:
                    file.write(l)
                    accepeted = True

            if accepeted:
                accepeted_pair += 1
            else:
                rejected_pair += 1

    print('number of pair', i)
    print('length pair rejected', length_rejected, (length_rejected/(i+1))*100, '%' )
    print('accepeted_pair', accepeted_pair, (accepeted_pair/(i+1))*100, '%' )
    print('rejected pair', rejected_pair,  (rejected_pair/(i+1))*100, '%' )




if __name__ == '__main__':
    """
    We remove all line with length query < length qubject
    Additionnally we remove  % query coverage per subject < threshold
    to be sure that the match is covering almost all the seq
    """

    blast_result = sys.argv[1]

    output_dir =  sys.argv[2]#'data/blast_result/Alphavirus_blast_result_all_vs_all_evalue_1e-5_coverage%i.out' % threshold

    threshold_min = int(sys.argv[3])
    threshold_max = int(sys.argv[4]) +1
    threshold_int = int(sys.argv[5])

    thresholds = range(threshold_min, threshold_max, threshold_int)
    out_files={}
    for t in thresholds:
        file_name = os.path.join(output_dir, 'coverage%i.out' % t)
        out_files[t] = open(file_name, 'w')

    filter_by_coverage(blast_result, out_files)

    for f in out_files.values():
        f.close()

import sys, os



def filter_by_evalue_and_coverage(blast_result, out_files, evalues, coverages):
    accepeted_pair = 0
    rejected_pair = 0
    length_rejected = 0

    with open(blast_result, 'r') as result_reader:
        for i, l in enumerate(result_reader):
            if i%1000000 == 0:
                print(i)
            qseqid, sseqid, qcovs, qcovhsp, evalue, bitscore = l.split('\t')

            len_query = float(qseqid.split('|')[2])
            len_subject = float(sseqid.split('|')[2])

            evalue = float(evalue)

            if len_query < len_subject:
                length_rejected += 1
                continue
            accepeted =False
            accepeted_evalues = []
            accepeted_coverage = []
            for  evalue_threshold in evalues:
                if evalue <= evalue_threshold:
                    accepeted_evalues.append(evalue_threshold)

            for coverage_threshold in coverages:
                if int(qcovs) >= coverage_threshold:
                    accepeted_coverage.append(coverage_threshold)

            for coverage in accepeted_coverage:
                for acc_evalue in accepeted_evalues:
                    out_files[(acc_evalue, coverage)].write(l)
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
    Blast all vs all : each of sequence are aligned twice one as a subject-query the other one as query-subject.
    We remove all line with length query < length qubject
    Additionnally we remove  % query coverage per subject < threshold
    to be sure that the match is covering almost all the seq
    """

    blast_result = sys.argv[1]

    output_dir =  sys.argv[2]#'data/blast_result/Alphavirus_blast_result_all_vs_all_evalue_1e-5_coverage%i.out' % threshold

    coverages =  sys.argv[3]

    print('coverage', coverages) # 20 30 40 50
    coverages = [int(c) for c in coverages.split(' ')]
    print(coverages)


    evalues = sys.argv[4]

    print('evalues', evalues)
    evalues = [float(e) for e in evalues.split(' ')]
    print(evalues)

    out_files= {}
    for c in coverages:
        for e in evalues:
            file_name = os.path.join(output_dir, 'evalue_{}coverage{}.out'.format(e, c))
            out_files[(e, c)] = open(file_name, 'w')

    filter_by_evalue_and_coverage(blast_result, out_files, evalues, coverages)

    for f in out_files.values():
        f.close()

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.metrics import roc_auc_score, roc_curve

from os import path,listdir

from multiple_alignment_analysis import parse_aln_directory_name


def getTPR_FPR_lists(positive, negative):
    score_max = max(positive+negative)
    thresholds = np.linspace(0, score_max+1, num=score_max*10)

    distribution_positive = get_distribution(thresholds, positive)
    distribution_negative = get_distribution(thresholds, negative)

    #Total
    total_negative = np.sum(distribution_negative)
    total_positive = np.sum(distribution_positive)
    #Cumulative sum
    cum_TP = 0
    cum_FP = 0
    #TPR and FPR list initialization
    TPR_list=[]
    FPR_list=[]
    #Iteratre through all values of x
    for i in range(len(thresholds)):
        cum_TP+=distribution_positive[i]
        cum_FP+=distribution_negative[i]

        FPR=cum_FP/total_negative
        TPR=cum_TP/total_positive
        TPR_list.append(TPR)
        FPR_list.append(FPR)
    # print(FPR_list)
    # print(TPR_list)
    #Calculating AUC, taking the 100 timesteps into account
    # auc=np.sum(TPR_list)/len(TPR_list)
    #Plotting final ROC curve
    # ax.plot(FPR_list, TPR_list, 'ro')
    # ax.plot(FPR_list, TPR_list)

    # ax.legend(["AUC=%.3f"%auc])
    return FPR_list,TPR_list, thresholds

def get_distribution(x, list_score):
    list_score.sort()
    iter_score = iter(list_score)
    score = next(iter_score)
    distribution = []
    for v in x:
        count = 0
        while score <= v and iter_score:
            count += 1
            try:
                score = next(iter_score)
            except StopIteration:
                iter_score = None
        distribution.append(count)

    return distribution

def parse_negative_positive_file(file):
    negative, positive = False, False
    with open(file, "r") as fl:
        for l in fl:
            if l.startswith('#'):
                continue
            case, score_str = l.split('=')
            if case == 'negative':
                negative = [float(s) for s in score_str.split(',')]
            elif case == 'positive':
                positive = [float(s) for s in score_str.split(',')]
    if not negative or not positive:
        raise ValueError(f'File format is not correct: {file}')
    return positive, negative

def plot_roc_curve(display_list):
    fig, ax = plt.subplots(1,1, figsize=(10,5))


    ax.set_xlim([-0.01,1.01])
    ax.set_ylim([-0.01,1.01])
    # ax.set_title("ROC Curve", fontsize=14)
    ax.set_ylabel('True Positive Rate', fontsize=12)
    ax.set_xlabel('False Positive Rate', fontsize=12)
    ax.grid()

    curves = []
    legends = []
    for info_dict in display_list:
        TPR_list = info_dict["TPR_list"]
        FPR_list = info_dict["FPR_list"]
        optimal_threshold_i = info_dict["optimal_threshold_i"]
        auc = info_dict['AUC']

        legend = f'{ info_dict["legend"]} AUC={auc:.2f}'
        legend = f'AUC={auc:.2f}'
        # Plotting final ROC curve
        curve = ax.plot(FPR_list, TPR_list)
        ax.plot(FPR_list[optimal_threshold_i], TPR_list[optimal_threshold_i], 'ro')
        # print([FPR_list[optimal_threshold_i]], [TPR_list[optimal_threshold_i]])
        ax.annotate(f'{info_dict["optimal_threshold"]:.2f}', (FPR_list[optimal_threshold_i], TPR_list[optimal_threshold_i]- 0.05))
        curves.append(curve)
        legends.append(legend)
    print(curves)
    print(legends)
    ax.legend(legends)
    ax.plot([0,1],[0,1], "--")
    txt_param = (f'Clustering parameters:\n'
        f"mcl inflation={info_dict['inflation'].replace('_', '.')}\n"
        f"min coverage={info_dict['coverage']}\n"
        f"bast e-value threshold={info_dict['evalue']}")

    bbox_props = dict(boxstyle="square,pad=0.3", fc="white", ec="black",)
    t = ax.text(0.67,0.15, txt_param, va="center",
        size=11,
        bbox=bbox_props)

    return fig, ax

if __name__ == '__main__':
    # file1= 'data/alignment/Viruses/RefSeq_download_date_2018-08-13/Viruses_evalue_1e-20coverage20_I1_4/stat/cross_validation_positive_negative.txt'
    # file2='data/alignment/Viruses/RefSeq_download_date_2018-07-21/Viruses_evalue_1e-50coverage60_I1_4/stat/cross_validation_positive_negative.txt'
    global_alignment_dir = 'data/alignment/Viruses/RefSeq_download_date_2018-08-13/'
    # aln_directories = [aln_directories[0]]
    legend = []
    display_list = []
    display_list_scikit=[]
    aln_dirs = listdir(global_alignment_dir)
    aln_dirs = ['Viruses_evalue_1e-40coverage70_I2/', 'Viruses_evalue_1e-40coverage40_I2/', 'Viruses_evalue_1e-40coverage20_I2/']
    aln_dirs = ['Viruses_evalue_1e-50coverage60_I1_4/']
    aln_dirs = ['Viruses_evalue_1e-20coverage60_I1_8']

    for dir in aln_dirs:
        # print(dir)
        info_dict = parse_aln_directory_name(dir)
        file=path.join(global_alignment_dir, dir, 'stat/cross_validation_positive_negative.txt')
        print(file)

        if not path.exists(file):
            raise ValueError(f'No "stat/cross_validation_positive_negative.txt" in {dir} ')

        positive, negative = parse_negative_positive_file(file)
        info_dict['positive'] = positive
        info_dict["negative"] = negative
        info_dict['score_max'] = max(positive+negative)
        display_list.append(info_dict)


    # positive = [0.1,0.2,0.3,0.4,0.5]
    # negative = [0.4,0.6,0.8,1.0]

    # positive = [s +1 for s in positive ]
    # negative = [s +1 for s in negative ]
    #
    # info_dict['positive'] = positive
    # info_dict["negative"] = negative
    # info_dict['score_max'] = max(positive+negative)

    for info_dict in display_list:
        info_dict['legend']  = (
        f"I={info_dict['inflation'].replace('_', '.')}"
        f"C={info_dict['coverage']}"
        f"E={info_dict['evalue']}"
        f"S={info_dict['cluster_with_isoforms_splitted']}")

        positive = info_dict['positive']
        negative = info_dict["negative"]

        # FPR_list,TPR_list, thresholds = getTPR_FPR_lists(positive, negative)
        # info_dict['TPR_list'] = TPR_list
        # info_dict['FPR_list'] = FPR_list
        # info_dict['AUC'] = np.sum(TPR_list)/len(TPR_list)
        # print('thr tpr fpr')
        # for i in range(len(thresholds)):
        #     print(thresholds[i], TPR_list[i], FPR_list[i])


        # print('## USE of scikit learn')
        scikit_info_dict = {}
        scikit_info_dict.update(info_dict)

        labels = [0]*len(positive) + [1]*len(negative)
        scores = positive + negative
        # print(labels)
        # print(scores)
        tpr, fpr,  thresholds_scikit = roc_curve(labels, scores, pos_label=0, drop_intermediate=False )
        # print('thr tpr fpr')
        # thresholds_scikit = thresholds_scikit[::-1]
        # fpr = tpr[::-1]
        # tpr = fpr[::-1]
        # for i in range(len(thresholds_scikit)):
        #     print(thresholds_scikit[i], tpr[i], fpr[i])
        # print(thresholds_scikit)
        scikit_info_dict['TPR_list'] = tpr
        scikit_info_dict['FPR_list'] = fpr
        # print(tpr)
        # print(fpr)
        youden_index_list = [tpr_i - fpr_i for tpr_i, fpr_i in zip(tpr, fpr)]
        optimal_threshold_i =  youden_index_list.index(max(youden_index_list))
        scikit_info_dict['optimal_threshold_i'] = optimal_threshold_i

        optimal_threshold = thresholds_scikit[optimal_threshold_i]
        scikit_info_dict['optimal_threshold'] = optimal_threshold # print(optimal_threshold)
        scikit_info_dict['AUC'] =  roc_auc_score(labels, scores )
        scikit_info_dict['legend']  = info_dict['legend']+'scikit '
        display_list_scikit.append(scikit_info_dict)

    # input()
    display_list_scikit = sorted(display_list_scikit, key=lambda x: x['AUC'], reverse=False)
    for info in display_list_scikit:
        print(info["AUC"], info['legend'], info['optimal_threshold'] )

    display_list = display_list_scikit
    fig, ax = plot_roc_curve(display_list)

    file_name = f'plot_test/roc_{dir}.png'
    fig.savefig("plot_test/roc.png")
    fig.show()

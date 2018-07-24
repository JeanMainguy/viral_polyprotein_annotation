# library
from pandas import DataFrame, read_csv
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
from numpy import log10

csv_file = "results/clustering_evaluation/expected_taxon_cluster_evaluation.csv"
# First way to call the 2 group Venn diagram:
venn2(subsets = (606, 33, 425), set_labels = ('Expected Genomes', 'Annotated Genomes'))
# plt.show()

inflation = 2

df = read_csv(csv_file, sep='\t')
print(df.loc[1])
# df = df.sort_values(by='sensitivity', ascending=False)
for index, row in df.iterrows():
    print(row['precision'])
    clustered_legend = 'I:{} Evalue:1e{} Cov:{}'.format(row['inflation'], int(log10(row['Evalue'])), row['coverage'])
    # print(row['inflation'], row['Evalue'], row['coverage'])
    print(clustered_legend)
    set_info = (row['complement_expected_taxid'], row['complement_clustered_taxid_with_poly'], row['expected_and_clustered_intersection'])
    venn2(subsets = (set_info), set_labels = ('Expected Genomes', clustered_legend))
    plt.show()
    input()

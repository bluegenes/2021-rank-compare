import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt

sns.set_context("paper")
import matplotlib.ticker as ticker
import matplotlib.patches as mpatches

from tqdm import tqdm # progress bars :)
tqdm.pandas()

import sourmash
from sourmash.lca import lca_utils
from sourmash.tax import tax_utils
from sourmash.logging import notify

#cluster95 = "regtdbclustering/renamed_gtdb_genomic_k31_s1k_kSpider_clusters_95.0%.tsv.gz"
#cluster96 = "regtdbclustering/renamed_gtdb_genomic_k31_s1k_kSpider_clusters_96.0%.tsv.gz"
#cluster97 = "regtdbclustering/renamed_gtdb_genomic_k31_s1k_kSpider_clusters_97.0%.tsv.gz"
#cluster98 = "regtdbclustering/renamed_gtdb_genomic_k31_s1k_kSpider_clusters_98.0%.tsv.gz"
#cluster99 = "regtdbclustering/renamed_gtdb_genomic_k31_s1k_kSpider_clusters_99.0%.tsv.gz"
#cluster100 = "regtdbclustering/renamed_gtdb_genomic_k31_s1k_kSpider_clusters_100.0%.tsv.gz"

#cluster_thresholds = [95,96,97] #,99,100]
cluster_thresholds = [96] #,99,100]

taxonomy_csv = "gtdb-rs207.taxonomy.csv"
tax = pd.read_csv(taxonomy_csv)
tax['lineage'] = tax["superkingdom"] + ';' + tax["phylum"] + ';' + tax["class"] + ';' + tax["order"] + ';' + tax["family"] + ';' + tax["genus"] + ';' + tax["species"]
tax['smash_lin'] = tax['lineage'].apply(lambda x: lca_utils.make_lineage(x))
tax['split_ident'] = tax['ident'].str.split('.', expand=True)[0]
#make lineage dictionary from tax info
taxD = tax.set_index('split_ident').to_dict()['smash_lin']

# count idents in each species
countD = tax[['lineage']].value_counts().reset_index(name='gtdb-lin-count')
gtdb_lineage_count = countD.set_index('lineage').to_dict()['gtdb-lin-count']

def count_and_find_lca(row, lineages=taxD): #, counts=countD):
    all_idents = row['cluster_idents']
    ident_list = all_idents.split(',')
    row['cluster_len'] = len(ident_list)
    all_lineages=[]
    for ident in ident_list:
        idt = ident.rsplit('.')[0]
        # handle bug
        if idt.startswith('CF'):
            idt = "G" + idt
            print(idt)
        lineage = taxD[idt]
        all_lineages.append(lineage)
    lca_tree = lca_utils.build_tree(all_lineages)
    lca = lca_utils.find_lca(lca_tree)
    row['cluster_lca'] = lca
    row['cluster_lca_pretty'] = lca_utils.display_lineage(lca[0])
    row['lca_rank'] = lca[0][-1].rank
    return row

for cluster_percent in cluster_thresholds:
    notify(f"Clustering at {cluster_percent}%:")
    cluster_file = f"regtdbclustering/renamed_gtdb_genomic_k31_s1k_kSpider_clusters_{cluster_percent}.0%.tsv"
    clusters = pd.read_csv(cluster_file, sep='\t', header=None, index_col=False, names = ["cluster_idents"])
    # add taxonomy
    clusters_w_tax = clusters.progress_apply(count_and_find_lca, axis=1)
    # add counts from GTDB
    clusters_w_tax['gtdb_count'] = clusters_w_tax['cluster_lca_pretty'].map(gtdb_lineage_count)

    singletons = clusters_w_tax[clusters_w_tax['cluster_len'] == 1].shape[0] # count number of rows (clusters) with single member
    notify(f"  num singleton clusters: {singletons}")
    # how many have lca rank of species?
    for rank in lca_utils.taxlist(include_strain=False):
        cluster_rank_lca = clusters_w_tax[clusters_w_tax['lca_rank'] == rank]
        num_cluster_rank_lca = cluster_rank_lca.shape[0]
        notify(f"  {rank} lca: {num_cluster_rank_lca}")
        if rank == 'species':
            num_matching_gtdb = cluster_rank_lca[cluster_rank_lca['cluster_len'] == cluster_rank_lca['gtdb_count']]
            notify(f'At species level, num clusters matching gtdb assignments: {num_matching_gtdb.shape[0]}')
            num_matching_gtdb_ns = num_matching_gtdb[num_matching_gtdb["gtdb_count"] >1]
            notify(f'At species level, num clusters matching gtdb assignments (non-singleton): {num_matching_gtdb_ns.shape[0]}')
        else:
            if num_cluster_rank_lca < 10:
                cluster_lineages = "\n  ".join(cluster_rank_lca["cluster_lca_pretty"].to_list())
                notify(f"  {cluster_lineages}")




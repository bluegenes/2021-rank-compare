##
# Script to get the number of unique k-mers at each taxonomic rank.
# modified from lca rankinfo.
##

import sys
import argparse
import csv
from collections import defaultdict
from collections import Counter

import sourmash
from sourmash.lca import lca_utils
from sourmash.tax import tax_utils

#def get_lca_lineage(lineages):
#    for rank in tax_utils.ascending_taxlist(include_strain=False):
#        lin_at_rank = set()
#        for lineage in lineages:
#            lin_at_rank.add(lca_utils.pop_to_rank(lineage, rank))
#        if len(lin_at_rank) == 1:
#            return lin_at_rank.pop()
#
#    return ""

def get_lca_lineage(lineages):
    lin_tree = lca_utils.build_tree(lineages)
    lca, node_count = lca_utils.find_lca(lin_tree)
    return lca, node_count


def hash2lin_from_LCA(dblist, min_num=2):
    """
    Collect counts of all the LCAs in the list of databases.
    CTB this could usefully be converted to a generator function.
    """
    # gather all hashval assignments from across all the databases
    assignments = defaultdict(set)
    for lca_db in dblist:
        for hashval, idx_list in lca_db.hashval_to_idx.items():
            if min_num and len(idx_list) < min_num:
                continue

            for idx in idx_list:
                lid = lca_db.idx_to_lid.get(idx)
                if lid is not None:
                    lineage = lca_db.lid_to_lineage[lid]
                    assignments[hashval].add(lineage)
    return assignments


def make_lca_lineageD(hashD):
   # to get unique k-mers at superkingdom, phylum, etc (generally) --> use lca rankinfo
   # here, get # unique k-mers at their LCA
    lca_lineage_counts = defaultdict(int)

    for n, (hashval, lineages) in enumerate(hashD.items()):
        # find LCA for this hashval
        lca_lineage, node_count = get_lca_lineage(lineages)
        if lca_lineage:
            # increment the count for this lineage in our LCA lineage counts dict
            lca_lineage_counts[lca_lineage] += 1
        else:
            lca_lineage_counts["no_lca"] += 1

    print(f"total unique hashes: {n}")
    return lca_lineage_counts, n


def main(args):
    # load all the databases
    dblist, ksize, scaled = lca_utils.load_databases(args.db, args.ksize, args.scaled)

    # count all the shared kmers across these databases
    #counts, total_kmer_count= count_shared_kmers(dblist)
    hashes_to_lineages = hash2lin_from_LCA(dblist)
    lca_lineageD, total_kmer_counts = make_lca_lineageD(hashes_to_lineages)

    print("LCA lineages found. Now sorting and counting...")
    sorted_lcaD = sorted(lca_lineageD.items(),  key=lambda kv: kv[1], reverse=True)
    with open(args.output_csv, 'w') as fp:
        ll_csv = csv.writer(fp)
        ll_csv.writerow(['rank', 'lineage', 'num_unique_kmers', 'f_unique_kmers']) #header

        for(lineage, count) in sorted_lcaD:
            if lineage == "no_lca":
                rank = "no_rank"
                lin = lineage
            else:
                rank = lineage[-1].rank
                lin = lca_utils.display_lineage(lineage)
            f_unique = float(count)/total_kmer_counts
            ll_csv.writerow([rank, lin, str(count), str(f_unique)])




def cmdline(sys_args):
    "Command line entry point w/argparse action."
    p = argparse.ArgumentParser()
    p.add_argument("--db", required=True, nargs="+")
    p.add_argument("--ksize", type=int)
    p.add_argument("--scaled", type=int)
    p.add_argument("--output-csv")
    args = p.parse_args()
    return main(args)

if __name__ == '__main__':
    returncode = cmdline(sys.argv[1:])
    sys.exit(returncode)



def test_get_lca_lineage():
    lin1 = lca_utils.make_lineage("a;b;c")
    lin2 = lca_utils.make_lineage("a;b;d")
    lin3 = lca_utils.make_lineage("a;d")

    lca_lin, num_nodes = get_lca_lineage([lin1, lin2])
    print(lca_lin)
    assert lca_lin == lca_utils.make_lineage("a;b")

    lca_lin, num_nodes = get_lca_lineage([lin1, lin2, lin3])
    print(lca_lin)
    assert lca_lin == lca_utils.make_lineage("a")


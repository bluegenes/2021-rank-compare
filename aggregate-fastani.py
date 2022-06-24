import os
import sys
import argparse
import glob
import pprint

import numpy as np
import pandas as pd

from collections import defaultdict, namedtuple
from itertools import combinations

# sourmash for tax/lca utils
from sourmash.tax import tax_utils
from sourmash.lca import lca_utils
from sourmash.logging import notify

fastani_lcares = namedtuple('fastani_lcares',
                           'comparison_name, a_name, b_name, lca_rank, lca_tax, avg_fastani_ident, avg_fastani_alignment_fraction')

def get_lineage(name, tax_assign):
    ident = tax_utils.get_ident(name) #, keep_identifier_versions=True)
    try:
        lineage = tax_assign[ident]
    except KeyError:
        if "GCF" in ident:
            ident = ident.replace("GCF", "GCA")
        elif "GCA" in ident:
            ident = ident.replace("GCA", "GCF")
        #print("new ident: ", ident)
        # ARGH, having 202 --> 207 lineage matching issues. just let it be None and notify
        lineage = tax_assign.get(ident, None)
    return lineage

def lineages_lca(linA, linB):
    for rank in tax_utils.ascending_taxlist():
        if lca_utils.is_lineage_match(linA, linB, rank):
            lca = lca_utils.pop_to_rank(linA, rank)
            return lca

def get_lca(name1, name2, lineages):
    lin1 = get_lineage(name1, lineages)
    lin2 = get_lineage(name2, lineages)
    if lin1 is None or lin2 is None:
        if lin1 is None:
            notify(f"{name1} is not in taxonomy files")
        if lin2 is None:
            notify(f"{name2} is not in taxonomy files")
        return None
    return lineages_lca(lin1,lin2)


def main(args):
    # taxonomy for lineages
    tax_assign = tax_utils.MultiLineageDB.load([args.taxonomy], keep_identifier_versions=False)

    # read in and format fastani results
    fastani = pd.read_csv(args.fastani_ranktax, sep = "\t", header=None, names=['a_fn','b_fn','fastani_ident','count_bidirectional_frag_mappings','total_query_frags'])
    fastani["a"] = fastani["a_fn"].str.rsplit("/", 1, expand=True)[1].str.rsplit("_genomic.fna.gz", 1, expand=True)[0]
    fastani["b"] = fastani["b_fn"].str.rsplit("/", 1, expand=True)[1].str.rsplit("_genomic.fna.gz", 1, expand=True)[0]
    fastani["comparison_name"] = fastani["a"] + "_x_" + fastani["b"]
    fastani.set_index("comparison_name",inplace=True)
    all_comparisons = fastani.index.to_list()

    # get all comparisons for this ranktax
    taxInfo = pd.read_csv(args.taxonomy, dtype=str, sep=',', header=0)

    # select comparison info
    taxInfo["rank_no_spaces"] = taxInfo[args.rank].str.replace(' ', '_')
    if args.ranktax_name and args.rank: # ranktax_name won't have spaces
        select_ranktax = args.ranktax_name
        taxInfo = taxInfo[taxInfo["rank_no_spaces"] == select_ranktax] # subset to t
    allranktax = taxInfo[args.rank].unique().tolist()

    results = []
    for rt in allranktax:
        acc_list = taxInfo[taxInfo["rank_no_spaces"] == select_ranktax]["ident"].to_list()

        # get all comparisons from list of accs
        for n, (a, b) in enumerate(combinations(acc_list, 2)):
            # get info for each comparison
            comparison_name = f"{a}_x_{b}"
            rev_name = f"{b}_x_{a}"
            # get lca rank/tax
            lca_lin = get_lca(a, b, tax_assign)
            if lca_lin is None:
                lca_rank = ""
                lca_tax = ""
            else:
                lca_rank = lca_lin[-1].rank
                lca_tax = lca_lin[-1].name

            # check if both comparison directions are available (>= 80% ani)
            if comparison_name in all_comparisons and rev_name in all_comparisons:
                # now grab directional fastani values. calc average ANI
                comp_ani = fastani.at[comparison_name, "fastani_ident"]
                rev_ani = fastani.at[rev_name, "fastani_ident"]
                avg_ani = np.mean([comp_ani, rev_ani])

                comp_fmap = fastani.at[comparison_name, "count_bidirectional_frag_mappings"]
                rev_fmap = fastani.at[rev_name, "count_bidirectional_frag_mappings"]

                comp_qfrag = fastani.at[comparison_name, "total_query_frags"]
                rev_qfrag = fastani.at[rev_name, "total_query_frags"]

                # get directional alignment fractions
                comp_af = float(comp_fmap)/float(comp_qfrag)
                rev_af = float(rev_fmap)/float(rev_qfrag)
                avg_af = np.mean([comp_af, rev_af])
                this_info = fastani_lcares(comparison_name, a, b, lca_rank, lca_tax, avg_ani, avg_af)
                results.append(this_info)

    # convert path ANI comparison info to pandas dataframe. to do: just write this with csv dictwriter to save on dict conversion
    aniDF = pd.DataFrame.from_records(results, columns = fastani_lcares._fields)

    # print to csv
    aniDF.to_csv(args.output_csv, index=False)
    print(f"done! ranktax {args.ranktax_name} comparison info written to {args.output_csv}")

def cmdline(sys_args):
    "Command line entry point w/argparse action."
    p = argparse.ArgumentParser()
    p.add_argument("--fastani-ranktax", required=True)
    p.add_argument("--taxonomy", default="gtdb-rs207.taxonomy.csv")
    p.add_argument("--rank") # comparison rank
    p.add_argument("--ranktax-name") # just do the one ranktax
    p.add_argument("--output-csv", required=True)
    args = p.parse_args()
    return main(args)

if __name__ == '__main__':
    returncode = cmdline(sys.argv[1:])
    sys.exit(returncode)


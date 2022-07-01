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

pyani_lcares = namedtuple('pyani_lcares',
                           'comparison_name, a_name, b_name, lca_rank, lca_tax, pyani_ident, pyani_coverage, pyani_aln_length, pyani_sim_errors, pyani_hadamard')

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

    # get basename for these sequences
    results = []
    taxInfo = pd.read_csv(args.taxonomy, dtype=str, sep=',', header=0)

    if "3" in str(args.pyani_version):
        len_fn= "matrix_aln_lengths_1.tab"
        cov_fn = "matrix_coverage_1.tab"
        had_fn = "matrix_hadamard_1.tab"
        id_fn = "matrix_identity_1.tab"
        se_fn = "matrix_sim_errors_1.tab"
    elif "2" in str(args.pyani_version):
        len_fn= "ANIb_alignment_lengths.tab"
        cov_fn = "ANIb_alignment_coverage.tab"
        had_fn = "ANIb_hadamard.tab"
        id_fn = "ANIb_percentage_identity.tab"
        se_fn = "ANIb_similarity_errors.tab"
    # select comparison info
    taxInfo["rank_no_spaces"] = taxInfo[args.rank].str.replace(' ', '_')
    if args.ranktax_name and args.rank: # ranktax_name won't have spaces
        select_ranktax = args.ranktax_name
        taxInfo = taxInfo[taxInfo["rank_no_spaces"] == select_ranktax] # subset to t
    allranktax = taxInfo[args.rank].unique().tolist()

    for rt in allranktax:
        acc_list = taxInfo[taxInfo["rank_no_spaces"] == select_ranktax]["ident"].to_list()

        results_dir = args.pyani_results_dir
        lenF = os.path.join(results_dir, len_fn)
        covF = os.path.join(results_dir, cov_fn)
        hadamardF = os.path.join(results_dir, had_fn)
        identF = os.path.join(results_dir, id_fn)
        simerrF = os.path.join(results_dir, se_fn)

        # read in all matrices
        lenD = pd.read_csv(lenF, sep="\t", header=0, index_col=0)
        covD = pd.read_csv(covF, sep="\t", header=0, index_col=0)
        hadD = pd.read_csv(hadamardF, sep="\t", header=0, index_col=0)
        idD = pd.read_csv(identF, sep="\t", header=0, index_col=0)
        seD = pd.read_csv(simerrF, sep="\t", header=0, index_col=0)

        # use headers on one file to get full column names, in case they're not just the accessions:
        names = lenD.columns.tolist()

        # get all comparisons from list of accs
        for n, (a, b) in enumerate(combinations(acc_list, 2)):
            # get info for each comparison
            comparison_name = f"{a}_x_{b}"
            # get lca rank/tax
            lca_lin = get_lca(a, b, tax_assign)
            if lca_lin is None:
                lca_rank = ""
                lca_tax = ""
            else:
                lca_rank = lca_lin[-1].rank
                lca_tax = lca_lin[-1].name

            # now grab pyani values
            pyani_ident, pyani_coverage, pyani_aln_length, pyani_sim_errors, pyani_hadamard = np.nan, np.nan, np.nan, np.nan, np.nan

            try:
                a_label = [x for x in names if x.startswith(a)][0]
                b_label = [x for x in names if x.startswith(b)][0]
            except IndexError:
                print(f"skipping comparison {a}_x_{b}")
                continue # skip this comparison


            # matrix will not be symmetric. Average values for each direction.
            pyani_identA = idD.at[a_label, b_label]
            pyani_identB = idD.at[b_label, a_label]
            pyani_ident = np.mean([pyani_identA, pyani_identB])

            pyani_coverageA = covD.at[a_label, b_label]
            pyani_coverageB = covD.at[b_label, a_label]
            pyani_coverage = np.mean([pyani_coverageA, pyani_coverageB])

            pyani_aln_lengthA = lenD.at[a_label, b_label]
            pyani_aln_lengthB = lenD.at[b_label, a_label]
            pyani_aln_length = np.mean([pyani_aln_lengthA, pyani_aln_lengthB])

            pyani_sim_errorsA = seD.at[a_label, b_label]
            pyani_sim_errorsB = seD.at[b_label, a_label]
            pyani_sim_errors = np.mean([pyani_sim_errorsA, pyani_sim_errorsB])

            pyani_hadamardA = hadD.at[a_label, b_label]
            pyani_hadamardB = hadD.at[b_label, a_label]
            pyani_hadamard = np.mean([pyani_hadamardA, pyani_hadamardB])

            this_info = pyani_lcares(comparison_name, a, b, lca_rank, lca_tax, pyani_ident, pyani_coverage, pyani_aln_length, pyani_sim_errors, pyani_hadamard)
            results.append(this_info)

    # convert path ANI comparison info to pandas dataframe. to do: just write this with csv dictwriter to save on dict conversion
    aniDF = pd.DataFrame.from_records(results, columns = pyani_lcares._fields)

    # print to csv
    aniDF.to_csv(args.output_csv, index=False)
    print(f"done! path comparison info written to {args.output_csv}")

def cmdline(sys_args):
    "Command line entry point w/argparse action."
    p = argparse.ArgumentParser()
    p.add_argument("pyani_results_dir")
    p.add_argument("--pyani-version", default = 'v0.3')
    p.add_argument("--taxonomy", default="gtdb-rs207.taxonomy.csv")
    p.add_argument("--rank") # comparison rank
    p.add_argument("--ranktax-name") # just do the one ranktax
    p.add_argument("--labels")
    p.add_argument("--output-csv", required=True)
    args = p.parse_args()
    return main(args)

if __name__ == '__main__':
    returncode = cmdline(sys.argv[1:])
    sys.exit(returncode)


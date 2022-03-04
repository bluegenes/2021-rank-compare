import os
import sys
import csv
import argparse
import glob
import pprint

import pandas as pd
from tqdm import tqdm
tqdm.pandas()

import screed
import sourmash
from sourmash.sourmash_args import load_file_as_signatures
import sourmash.lca.lca_utils as lca_utils
from sourmash.picklist import SignaturePicklist

from collections import defaultdict, namedtuple

CompareResult = namedtuple('CompareResult',
                           'comparison_name, sigA_name, sigB_name, lowest_common_rank, lowest_common_lineage, alphabet, ksize, scaled, jaccard, max_containment, A_containment, B_containment, sigA_hashes, sigB_hashes, intersect_hashes')

def compare_sigs(sigA, sigB, comparison_name, lowest_common_rank, lowest_common_lineage, alpha, ksize, scaled):
    A_name = str(sigA).split(" ")[0]
    B_name = str(sigB).split(" ")[0]
    sigA_numhashes = len(sigA.minhash.hashes)
    sigB_numhashes = len(sigB.minhash.hashes)
    intersect_numhashes = sigA.minhash.count_common(sigB.minhash)
    jaccard = sigA.jaccard(sigB)
    containA = sigA.contained_by(sigB)
    containB = sigB.contained_by(sigA)
    max_contain = sigA.max_containment(sigB)
    #max_contain = max(containA,containB)
    return CompareResult(comparison_name, A_name, B_name, lowest_common_rank, lowest_common_lineage, alpha, ksize, scaled, jaccard, max_contain, containA, containB, sigA_numhashes, sigB_numhashes, intersect_numhashes)

def find_lca(linA,linB, reverse_taxlist=['species', 'genus', 'family', 'order', 'class', 'phylum', 'superkingdom']):
    lca_rank=None
    lca_lin=None
    for rank in reverse_taxlist:
            # do we have a lca lineage match at this rank?
        if lca_utils.is_lineage_match(linA, linB, rank):
            lca_rank = rank
            lca_lin = lca_utils.pop_to_rank(linA, lca_rank)
            break
    return lca_rank, lca_lin


def make_lineages(row):
    rank = row['rank']
    lin = lca_utils.make_lineage(rank)
    row['lineage'] = lin
    return row

def main(args):
    ksize=args.ksize
    scaled=args.scaled
    alphabet=args.alphabet
    if alphabet == "nucleotide":
        moltype = "DNA"
    else:
        moltype = alphabet

    # read in taxonomy csv and get lineages into sourmash format
    print('reading taxonomy')
    gtdb_taxonomy = pd.read_csv(args.gtdb_taxonomy)
    gtdb_taxonomy.set_index('ident', inplace=True)
    cols = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    gtdb_taxonomy['rank'] = gtdb_taxonomy[cols].agg(';'.join, axis=1)
    gtdb_taxonomy = gtdb_taxonomy.progress_apply(make_lineages, axis=1)

    print(f'loading database {args.database}')
    # for ident in taxonomy csv, compare to all remaining. Keep list of "done" and exclude via picklist EXCLUDE
    gtdb_zip = args.database
    #gtdb_loaded = sourmash.load_file_as_signatures(db)
    gtdb = sourmash.index.ZipFileLinearIndex.load(gtdb_zip)

    # compare to rest of idents
    comparison_idents = gtdb_taxonomy.index[1:]
    #loop through comparisons
    output_fieldnames = CompareResult._fields # may need to be list?
    for n, ident in enumerate(gtdb_taxonomy.index):
        #print(f'starting comparisons with ident {ident}')
        if n % 1000 == 0:
            print(f"... comparing {n}th ident, {ident}\n")

        # select and load anchor ident
        picklist = SignaturePicklist('ident')
        picklist.init([ident])
        anchor_sig = next(gtdb.select(picklist=picklist, ksize=ksize, moltype=moltype).signatures())
        #anchor_sig = next(anchor_sig_select)

        # grab anchor lineage
        anchor_lin = gtdb_taxonomy.at[ident, 'lineage']

        # make output file for this anchor
        if not os.path.exists(args.output_dir):
            os.makedirs(args.output_dir, exist_ok=True)
        this_outfile = os.path.join(args.output_dir, f"{ident}.gtdb-compare.csv")
        with open(this_outfile, 'w') as outF:
            w = csv.DictWriter(outF, fieldnames=output_fieldnames)

            # load this sig and make a comparison
            comparisons = []
            for compare_ident in tqdm(comparison_idents):
                picklist = SignaturePicklist('ident')
                picklist.init([compare_ident])
                compare_sig = next(gtdb.select(picklist=picklist, ksize=ksize, moltype=moltype).signatures())
                #compare_sig = next(compare_sig_gen)
                compare_lin = gtdb_taxonomy.at[compare_ident, 'lineage']

                lca_rank, lca_lineage = find_lca(anchor_lin, compare_lin)
                lca_lin = lca_utils.display_lineage(lca_lineage)

                comparison = compare_sigs(anchor_sig, compare_sig, f"{ident}_x_{compare_ident}", lca_rank, lca_lin, alphabet, ksize, scaled)
                # write comparison to file
                comparisons.append(comparison)
            for c in comparisons:
                w.writerow(c._asdict())

        # remove top ident, which will be next anchor
        comparison_idents = comparison_idents[1:]


def cmdline(sys_args):
    "Command line entry point w/argparse action."
    p = argparse.ArgumentParser()
    p.add_argument("--gtdb_taxonomy", default="/group/ctbrowngrp/gtdb/gtdb-rs202.taxonomy.v2.with-strain.csv")
    p.add_argument("--database", default="/group/ctbrowngrp/gtdb/databases/gtdb-rs202.protein.k10.zip")
    p.add_argument("--alphabet", default="protein")
    p.add_argument("--ksize", default=10, type=int)
    p.add_argument("--scaled", default=100, type=int)
    p.add_argument("--output-dir", default = "output.gtdb-rs202-compare")
    args = p.parse_args()
    return main(args)

if __name__ == '__main__':
    returncode = cmdline(sys.argv[1:])
    sys.exit(returncode)

import os
import sys
import csv
import argparse
import glob
import pprint

import pandas as pd

import screed
import sourmash
#import sourmash.lca.lca_utils as lca_utils
#import sourmash.tax.tax_utils as tax_utils
from sourmash.tax import tax_utils
from sourmash.lca import lca_utils
from collections import defaultdict



def main(args):

    ascending_taxlist = list(tax_utils.ascending_taxlist(include_strain=True))
    rankinfo_db = defaultdict(int)
    total_current_kmers=0
    current_db = ''
    with open(args.rankinfo, 'r') as rinfo:
        with open(args.output, 'w') as outF:
            header = ["database","ksize","alphabet","alpha_ksize","scaled","rank","unique_rank_kmers","f_rank_kmers","aggregated_kmers","f_aggregated_kmers"]
            outF.write(','.join(header) + "\n")
            for line in rinfo:
                # skip header
                if line.startswith("database"):
                    continue
                database,ksize,alphabet,scaled,rank,kmers,percent = line.strip().split(',')
                kmers=int(kmers)
                rankinfo_db[rank] = kmers
                total_current_kmers += kmers
                # assume input goes superkingdom --> strain; aggregate k-mers for this db
                # this will be inaccurate if this is not an accurate assumption
                if rank == 'strain':
                    print(f'aggregating info for database {database}')
                    # sum k-mers for this db
                    aggregated_kmers = 0
                    alpha_ksize = f"{alphabet}-k{ksize}"
                    for rank in ascending_taxlist:
                        unique_rank_kmers = rankinfo_db[rank]
                        f_rank_kmers = float(unique_rank_kmers)/total_current_kmers
                        aggregated_kmers += rankinfo_db[rank]
                        f_aggregated_kmers = float(aggregated_kmers)/total_current_kmers
                        outF.write(f"{database},{ksize},{alphabet},{alpha_ksize},{scaled},{rank},{unique_rank_kmers},{f_rank_kmers},{aggregated_kmers},{f_aggregated_kmers}\n")

                    # reset all
                    total_current_kmers = 0
                    rankinfo_db = defaultdict(int)


def cmdline(sys_args):
    "Command line entry point w/argparse action."
    p = argparse.ArgumentParser()
    p.add_argument("rankinfo")
    p.add_argument("--output")
    args = p.parse_args()
    return main(args)

if __name__ == '__main__':
    returncode = cmdline(sys.argv[1:])
    sys.exit(returncode)

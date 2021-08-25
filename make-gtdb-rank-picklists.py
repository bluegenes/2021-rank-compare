import os
import sys
import argparse

import pandas as pd
from sourmash.tax import tax_utils


def main(args):
    if not os.path.isdir(args.output_dir):
        os.mkdir(args.output_dir)
        for rank in tax_utils.ascending_taxlist(include_strain=False):
            os.mkdir(os.path.join(args.output_dir, f"{rank}"))

    taxDF = pd.read_csv(args.gtdb_taxonomy)
    for rank in tax_utils.ascending_taxlist(include_strain=False):
        rank_taxonomies = taxDF[rank].unique()
        for rt in rank_taxonomies:
            outfile = os.path.join(args.output_dir, f"{rank}", args.output_base + f".{rt}.csv")
            outfile = outfile.replace(' ', '_')
            subset_df = taxDF[taxDF[rank] == rt]
            subset_df.to_csv(outfile, index=False)


def cmdline(sys_args):
    "Command line entry point w/argparse action."
    p = argparse.ArgumentParser()
    p.add_argument("--gtdb-taxonomy",  default = "/group/ctbrowngrp/gtdb/gtdb-rs202.taxonomy.v2.csv")
    p.add_argument("--output-dir", default="gtdb-rs202-taxonomic-picklists")
    p.add_argument("--output-base", default="gtdb-rs202")
    args = p.parse_args()
    return main(args)

if __name__ == '__main__':
    returncode = cmdline(sys.argv[1:])
    sys.exit(returncode)

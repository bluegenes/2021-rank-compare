import os
import sys
import argparse

import numpy as np
import pandas as pd






def main(args):

    tax_info = pd.read_csv(args.taxonomy_csv, header=0)
    #For metadata info, pandas needs low_memory=False(DtypeWarning: Columns (61,65,74,82,83) have mixed types.Specify dtype option on import or set low_memory=False)
    metadata_info = pd.read_csv(args.metadata_csv, header=0, low_memory=False)
    representative_accessions = metadata_info[metadata_info["gtdb_representative"] == "t"]["accession"].str.replace("GB_", "").str.replace("RS_", "")
    tax_info['is_representative'] = np.where(tax_info["ident"].isin(representative_accessions), True, False)
    tax_info["signame"] = tax_info["ident"] + " " + tax_info["species"].str.replace("s__", "")
    tax_info.to_csv(args.output, index=False)


def cmdline(sys_args):
    "Command line entry point w/argparse action."
    p = argparse.ArgumentParser()
    p.add_argument("--taxonomy_csv",  default = "/group/ctbrowngrp/gtdb/gtdb-rs202.taxonomy.v2.csv")
    p.add_argument("--metadata_csv",  default = "/group/ctbrowngrp/gtdb/gtdb-rs202.metadata.csv.gz")
    p.add_argument("--output", default="gtdb-rs202.taxonomy.with-repinfo.csv")
    args = p.parse_args()
    return main(args)

if __name__ == '__main__':
    returncode = cmdline(sys.argv[1:])
    sys.exit(returncode)

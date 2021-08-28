import os
import sys
import screed
import argparse

import sourmash
from sourmash.nodegraph import Nodegraph, calc_expected_collisions#, extract_nodegraph_info

DEFAULT_N_TABLES = 4
DEFAULT_MAX_TABLESIZE = int(1e8)
DEFAULT_MAX_FALSE_POS = 0.2 # this is the default in https://github.com/sourmash-bio/sourmash/blob/latest/src/sourmash/nodegraph.py

def screed_open_fasta(fasta, strict_mode=False):
    records = []
    try:
        records = screed.open(fasta)
    except ValueError:
        print(f"File {fasta} doesn't seem to be a sequence file.\n")
        if strict_mode:
            print("Exiting.\n")
            sys.exit(-1)
        else:
            print("Skipping...\n")

    return records


def main(args):
    # handle alpha
    is_protein=True
    is_dayhoff=False
    is_hp=False

    tablesize = int(float(args.tablesize))

    alphabet = args.alphabet
    if alphabet in ['dna', 'rna', 'nucleotide']:
        is_protein=False
    elif alphabet == "dayhoff":
        is_dayhoff = True
    elif alphabet == "hp":
        is_hp = True

    # init bf
    bloom_filter = Nodegraph(args.ksize, tablesize, args.n_tables)
    mh = sourmash.MinHash(n=0, ksize=args.ksize, scaled=1, is_protein=is_protein, dayhoff=is_dayhoff, hp=is_hp)
    for fasta in args.input_files:
        records = screed_open_fasta(fasta, strict_mode=False)
        for record in records:
            if "*" in record.sequence:
                continue
            if is_protein:
                mh.add_protein(record.sequence)
            else:
                #handle invalid DNA chars that may exist in sequence
                if 'N' in record.sequence:
                    seqs = record.sequence.split('N')
                    for seq in seqs:
                        mh.add_sequence(seq)
                else:
                    mh.add_sequence(record.sequence)
            # upd after each seq? Or build mh of entire file, then update?
            bloom_filter.update(mh)
            mh.clear()

    # calculate expected collisions
    bf_info = calc_expected_collisions(bloom_filter, force=True, max_false_pos=args.max_false_pos)

    bloom_filter.save(args.output)




def cmdline(sys_args):
    "Command line entry point w/argparse action."
    p = argparse.ArgumentParser()
    p.add_argument("input_files", nargs='+')
    p.add_argument("-k", "--ksize",  type=int)
    p.add_argument("--alphabet")
    p.add_argument("--n_tables", type=int, default=DEFAULT_N_TABLES)
    p.add_argument("--tablesize", default=DEFAULT_MAX_TABLESIZE)
    p.add_argument("--max_false_pos", default=DEFAULT_MAX_FALSE_POS)
    p.add_argument("--output")
    args = p.parse_args()
    return main(args)

if __name__ == '__main__':
    returncode = cmdline(sys.argv[1:])
    sys.exit(returncode)


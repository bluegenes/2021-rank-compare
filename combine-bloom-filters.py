import os
import sys
import screed
import argparse

import sourmash
from sourmash.nodegraph import Nodegraph#, extract_nodegraph_info, calc_expected_collisions

DEFAULT_N_TABLES = 4
DEFAULT_MAX_TABLESIZE = int(1e8)
DEFAULT_MAX_FALSE_POS = 0.2 # this is the default in https://github.com/sourmash-bio/sourmash/blob/latest/src/sourmash/nodegraph.py


def main(args):
    ksize = args.ksize
    tablesize = int(float(args.tablesize))

   # alphabet = args.alphabet
   # if alphabet in ['dna', 'rna', 'nucleotide']:
   #     is_protein=False
   # elif alphabet == "dayhoff":
   #     is_dayhoff = True
   # elif alphabet == "hp":
   #     is_hp = True

    # make the full list of bf's to join
    all_bfs = []
    if args.input_files:
        all_bfs.extend(args.input_files)
    if args.from_file:
        ff = args.from_file
        file_bfs = [x.strip() for x in open(ff, 'r')]
        all_bfs.extend(file_bfs)

    # make sure no duplicates
    all_bfs = list(set(all_bfs))

    combined_bf = Nodegraph(args.ksize, tablesize, args.n_tables)
    sys.stderr.write(f"Building sourmash nodegraph, tablesize: {tablesize}; n_tables: {args.n_tables}\n")

    for bf_fname in all_bfs:
        #load bf and check that ksizes match
        bf = Nodegraph.load(bf_fname)
        this_ksize = int(bf.ksize())
        if this_ksize != ksize:
            print(f"Nodegraph {bf_fname} has ksize {combined_bf.ksize}, which does not match desired ksize {ksize}. Skipping ...")
            continue
        # update combined bf with this bf
        combined_bf.update(bf)
        expected_collisions = combined_bf.expected_collisions
        sys.stderr.write(f"Added ng {bf_fname}; current expected collisions: {expected_collisions}\n")
        if expected_collisions >= args.max_false_pos:
            sys.stderr.write('ERROR, FP rate {expected_collisions} is too high. Please increase tablesize.')
            print('ERROR, FP rate is too high. Please increase tablesize.')
            sys.exit(-1)

    sys.stderr.write(f"Nodegraph details:  tablesize: {tablesize}; n_tables: {args.n_tables}\n")
    sys.stderr.write(f"Expected collisions: {expected_collisions}\n")
    combined_bf.save(args.output)

    # once we've integrated all, save the combined bf
    combined_bf.save(args.output)


def cmdline(sys_args):
    "Command line entry point w/argparse action."
    p = argparse.ArgumentParser()
    p.add_argument("input_files", nargs='*')
    p.add_argument("--from-file")
    p.add_argument("-k", "--ksize",  type=int)
    #p.add_argument("--alphabet")
    p.add_argument("--n_tables", type=int, default=DEFAULT_N_TABLES)
    p.add_argument("--tablesize", default=DEFAULT_MAX_TABLESIZE)
    p.add_argument("--max_false_pos", default=DEFAULT_MAX_FALSE_POS)
    p.add_argument("--output")
    args = p.parse_args()
    return main(args)

if __name__ == '__main__':
    returncode = cmdline(sys.argv[1:])
    sys.exit(returncode)


import os
import sys
import screed
import argparse

import khmer
import sourmash
#from sourmash.nodegraph import Nodegraph, extract_nodegraph_info, calc_expected_collisions
#from sourmash.sbt import GraphFactory #SBT, GraphFactory, Leaf, Node
#from tqdm import tqdm

# sourmash GraphFacotry (~ khmer Nodegraph) features
DEFAULT_N_TABLES = 4
DEFAULT_MAX_TABLESIZE = int(1e8)


# GraphFactory is bloom filter generator
#factory = GraphFactory(ksize, starting_size(bytes), num_nodegraph_tables)
#factory = GraphFactory(ksize, tablesize, n_tables)
#factory = GraphFactory(31, 1e5, 4)
#peptide_bloom_filter = khmer.Nodegraph(peptide_ksize, tablesize, n_tables=n_tables)



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



#def khmer_bloom_filter(fasta_files,ksize,molecule,n_tables=DEFAULT_N_TABLES,tablesize=DEFAULT_MAX_TABLESIZE,is_protein=True,index_dir=None, strict=False):
#    """ Build GraphFactory bloom filter """
#    factory = GraphFactory(ksize, tablesize, n_tables)



def main(args):
    # handle alpha
    is_protein=True
    alphabet = args.alphabet
    if alphabet in ['dna', 'rna', 'nucleotide']:
        is_protein=False
    # ksize
    ksize = args.ksize

    # init bf
    peptide_bloom_filter = khmer.Nodegraph(ksize, args.tablesize, n_tables=args.n_tables)
    mh = sourmash.MinHash(n=0, ksize=ksize, scaled=1, is_protein=is_protein)
    for fasta in args.input_files:
        records = screed_open_fasta(fasta, strict_mode=False)
        #for record in tqdm(records):
        for record in records:
            if "*" in record["sequence"]:
                continue
            for hashval in mh.seq_to_hashes(record.sequence, is_protein=is_protein):
                peptide_bloom_filter.add(hashval)

                #for i in range(0, len(record) - ksize + 1):
                #    kmer = dnaseq[i:i+K]
                #hashes = seq_to_hashes(is_protein=is_protein)
    peptide_bloom_filter.save(args.output)


def cmdline(sys_args):
    "Command line entry point w/argparse action."
    p = argparse.ArgumentParser()
    p.add_argument("input_files", nargs='+')
    p.add_argument("-k", "--ksize",  type=int)
    p.add_argument("--alphabet")
    p.add_argument("--n_tables", default=DEFAULT_N_TABLES)
    p.add_argument("--tablesize", default=DEFAULT_MAX_TABLESIZE)
    p.add_argument("--output")
    args = p.parse_args()
    return main(args)

if __name__ == '__main__':
    returncode = cmdline(sys.argv[1:])
    sys.exit(returncode)


##
# once we have bf's, should be able to count directly from those, right?
# in the mean time...
# From a picklist csv of genomes/proteomes,
# build Counter of k-mers found in each genome/proteome.
# output as a csv.gz? Plus a report of # shared at each level.
##
import os
import sys
import screed
import argparse
from collections import Counter, defaultdict

import sourmash
#from sourmash.nodegraph import Nodegraph, calc_expected_collisions#, extract_nodegraph_info


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

    alphabet = args.alphabet
    if alphabet in ['dna', 'rna', 'nucleotide']:
        is_protein=False
    elif alphabet == "dayhoff":
        is_dayhoff = True
    elif alphabet == "hp":
        is_hp = True

    # init k-mer dictionary
    #kmers = defaultdict(lambda: defaultdict(int))
    kmers = defaultdict(int)

    mh = sourmash.MinHash(n=0, ksize=args.ksize, scaled=1, is_protein=is_protein, dayhoff=is_dayhoff, hp=is_hp, track_abundance=True)

    input_fasta = [x.strip() for x in open(args.fasta_list, 'r')]

    for n, fasta in enumerate(input_fasta):
        records = screed_open_fasta(fasta, strict_mode=False)
        for x, record in enumerate(records):
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

            # at this point, we should have added all k-mers in this genome/proteome to the minhash.
            # add these to our defaultdict, then continue with all other genomes/proteomes in this picklist
            hashes_with_abundance = mh.hashes
            import pdb;pdb.set_trace()
            for hashval, abund in hashes_with_abundance.items():
                kmers[hashval] += abund

            # now clear and move on to next genome/proteome
            mh.clear()
        if n % 1000 == 0:
            sys.stderr.write(f"Processed {x} records in {n}th fasta file\n")

    # summarize shared k-mers at this level:
    num_unique_kmers = len(kmers.keys())
    basename  = os.path.basename(args.basename)
    with open(args.output_num_unique, 'w') as counts:
        counts.write(f"{basename},{num_unique_kmers}")
    # now optionally build dataframe from this info and save for future use
    if args.output_kmer_counts:
        hashDF = pd.DataFrame.from_dict(kmers)
        hashDF.to_csv(args.output_kmer_counts, index_label=basename)




def cmdline(sys_args):
    "Command line entry point w/argparse action."
    p = argparse.ArgumentParser()
    p.add_argument("fasta_list")
    p.add_argument("-k", "--ksize",  type=int)
    p.add_argument("--alphabet", required=True)
    p.add_argument("--basename", required=True)
    p.add_argument("--output-num-unique", required=True)
    p.add_argument("--output-kmer-counts")
    args = p.parse_args()
    return main(args)

if __name__ == '__main__':
    returncode = cmdline(sys.argv[1:])
    sys.exit(returncode)


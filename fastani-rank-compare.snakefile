import os, sys
import pandas as pd
import glob

from snakemake.workflow import srcdir

import csv
import pandas as pd
from collections import defaultdict

### goal: compare representative --> all family members ANI

out_dir = config.get('output_dir', 'output.fastani-compare')
logs_dir = os.path.join(out_dir, "logs")
reps_only = config.get('representatives_only', True)

# requires taxonomy with representitive info
print('reading taxonomy')
gtdb_taxonomy=config.get('gtdb_taxonomy', 'conf/gtdb-rs207.taxonomy.with-repinfo.csv')
taxDF = pd.read_csv(gtdb_taxonomy)
# let's subset for testing
taxDF = taxDF.tail(100000)

# get fastapaths
fa_info = pd.read_csv("/home/ntpierce/2021-rank-compare/gtdb-rs207.fromfile.csv") # cols: ident, name, genome_filename, protein_filename

# set ident as index for each
fa_info.set_index("ident", inplace=True)

# get comparison groups
compare_rank = config.get('rank', "family")
alltax_at_rank = taxDF[compare_rank].unique().tolist()

ranktaxD = defaultdict(list)
for ranktax in alltax_at_rank:
    # get query_idents
    subDF = taxDF[(taxDF[compare_rank] == ranktax)]
    all_accs = subDF["ident"].to_list()
    num_acc = len(all_accs)
    # testing:let's start smaller: ignore >500 comparisons
    if 2 <= num_acc <= 50: # can't compare a single genome :)
        rt = ranktax.replace(' ', '_') # no spaces plsss
        ranktaxD[rt] = all_accs

num_rt_to_compare = len(ranktaxD.keys())
print(f"found {num_rt_to_compare} comparison groups")

onstart:
    print("------------------------------")
    print("Estimate ANI for comparisons at given GTDB rank")
    print("------------------------------")

onsuccess:
    print("\n--- Workflow executed successfully! ---\n")

onerror:
    print("Alas!\n")


rule all:
    input: 
        # fastani
        expand(os.path.join(out_dir, "fastani", f"{compare_rank}.fastani.csv.gz")),


### fastANI rules ###
localrules: build_filepaths_for_fastani
rule build_filepaths_for_fastani:
    output: os.path.join(out_dir, "fastani","{ranktax}.filepaths.txt")
    run:
        with open(str(output), "w") as out:
            acc_list = ranktaxD[wildcards.ranktax]
            for acc in acc_list:
                fn = fa_info.at[acc, 'genome_filename']
                out.write(f"{fn}\n")


rule compare_via_fastANI:
    input: os.path.join(out_dir, "fastani", "{ranktax}.filepaths.txt")
    output: os.path.join(out_dir, "fastani","{ranktax}.fastani.tsv"),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *3000,
        time=1200,
        partition="med2",
    log: os.path.join(logs_dir, "fastani", "{ranktax}.fastani.log")
    benchmark: os.path.join(logs_dir, "fastani", "{ranktax}.fastani.benchmark")
    conda: "conf/envs/fastani.yml"
    shell:
        """
        fastANI --ql {input} --rl {input} -o {output} 2> {log}
        """


localrules: parse_fastani_ranktax
rule parse_fastani_ranktax:
    input: os.path.join(out_dir, "fastani", "{ranktax}.fastani.tsv")
    output: os.path.join(out_dir, "fastani", "{ranktax}.fastani.parsed.csv")
    params:
        compare_rank= compare_rank,
    log: os.path.join(logs_dir, "fastani", "{ranktax}.parse_fastani.log")
    benchmark: os.path.join(logs_dir, "fastani", "{ranktax}.parse_fastani.benchmark")
    shell:
        """
        python aggregate-fastani.py --fastani-ranktax {input} --rank {params.compare_rank} \
                                    --ranktax-name {wildcards.ranktax} \
                                    --output-csv {output} 2> {log}
        """

localrules: aggregate_fastani_results
rule aggregate_fastani_results:
    input: expand(os.path.join(out_dir, "fastani", "{ranktax}.fastani.parsed.csv"), ranktax=ranktaxD.keys())
    output: os.path.join(out_dir, "fastani", f"{compare_rank}.fastani.csv.gz"),
    run:
        # aggreate all csv.gzs --> single csv
        aggDF = pd.concat([pd.read_csv(str(csv), sep=",") for csv in input])
        aggDF.to_csv(str(output), index=False)

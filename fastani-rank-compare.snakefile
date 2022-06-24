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
taxDF = taxDF.tail(60000)

# get fastapaths
fa_info = pd.read_csv("/home/ntpierce/2021-rank-compare/gtdb-rs207.fromfile.csv") # cols: ident, name, genome_filename, protein_filename

# set ident as index for each
fa_info.set_index("ident", inplace=True)

# get comparison groups
compare_rank = config.get('rank', "species")
alltax_at_rank = taxDF[compare_rank].unique().tolist()

#taxDF.set_index("ident", inplace=True)
all_comparisons = []
rankcompD = defaultdict(list)
for ranktax in alltax_at_rank:
    # get query_idents
    subDF = taxDF[(taxDF[compare_rank] == ranktax)]
    if reps_only:
        query_accs = subDF[subDF['is_representative'] == True]["ident"].to_list()
        compare_accs = subDF[subDF['is_representative'] == False]["ident"].to_list()
    else:
    # faster to just use itertools for this? probably can't run anyway, too many , so shrug.
        query_accs = subDF["ident"].apply(list)
        compare_accs = query_accs 
    these_comparisons = expand("{query}_x_{compare}", query=query_accs, compare=compare_accs)
    all_comparisons.extend(these_comparisons)
    rankcompD[ranktax].extend(these_comparisons)

#onstart:
#    print("------------------------------")
#    print("Estimate ANI for comparisons at given GTDB rank"
#    print("------------------------------")

onsuccess:
    print("\n--- Workflow executed successfully! ---\n")

onerror:
    print("Alas!\n")


rule all:
    input: 
        # pyani
        #os.path.join(out_dir, "pyani", f"{basename}.pyani-ANIm.csv.gz"),
        #os.path.join(out_dir, "pyani", f"{basename}.pyani-ANIb.csv.gz"),
        # fastani
        expand(os.path.join(out_dir, "fastani", f"{compare_rank}.fastani.csv.gz")),
        #expand(os.path.join(out_dir, "fastani-compare", "{basename}.fastani.tsv"), basename=basename),


### fastANI rules ###
#localrules: build_filepaths_for_fastani
#rule build_filepaths_for_fastani:
#    output: os.path.join(out_dir, "fastani", "{path}", "{path}.filepaths.txt")
#    run:
#        with open(str(output), "w") as out:
#            acc_list = path2acc[wildcards.path]
#            for acc in acc_list:
#                fn = fa_info.at[acc, 'genome_filename']
#                out.write(f"{fn}\n")


rule compare_via_fastANI:
    input:  
        rep_genome = lambda w: fa_info.at[w.rep, "genome_filename"],
        compare_genome = lambda w: fa_info.at[w.cmp, "genome_filename"],
    output: os.path.join(out_dir, "fastani", "comparisons", "{rep}_x_{cmp}.fastani.tsv"),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *3000,
        time=1200,
        partition="low2",
    log: os.path.join(logs_dir, "fastani", "comparisons/{rep}_x_{cmp}.fastani.log")
    benchmark: os.path.join(logs_dir, "fastani", "comparisons/{rep}_x_{cmp}.fastani.benchmark")
    conda: "conf/envs/fastani.yml"
    shell:
        """
        fastANI -q {input.rep_genome:q} -r {input.compare_genome:q} -o {output} > {log} 2>&1
        """

localrules: write_fastani_result_csv
rule write_fastani_result_csv:
    input: expand(os.path.join(out_dir, "fastani", "comparisons/{comparisons}.fastani.tsv"), comparisons=all_comparisons)
    output: os.path.join(out_dir, "fastani", f"{compare_rank}.fastani.filecsv"),
    run:
        with open(str(output), "w") as out:
            for inF in input:
                fn = os.path.abspath(str(inF))
                compInfo = str(inF).rsplit('.fastani.tsv')
                rep, comp= compInfo.split('_x_')
                out.write(f"{compInfo},{rep},{comp},{fn}\n")


localrules: aggregate_fastani_results
rule aggregate_fastani_results:
    input:
        fastani=os.path.join(out_dir, "fastani", f"{compare_rank}.fastani.filecsv")
    output: os.path.join(out_dir, "fastani", f"{compare_rank}.fastani.csv.gz"),
    log: os.path.join(logs_dir, "fastani", f"{compare_rank}.fastani.aggregate.log")
    benchmark: os.path.join(logs_dir, "fastani", f"{compare_rank}.fastani.aggregate.benchmark")
    shell:
        """
        python aggregate-fastani.py --fastani-filecsv {input.fastani} \
                                    --output-csv {output} > {log} 2>&1
        """

import os, sys
import pandas as pd
import glob

from snakemake.workflow import srcdir

import csv
import pandas as pd
from collections import defaultdict

### goal: compare representative --> all family members ANI

out_dir = config.get('output_dir', 'output.sourmash-compare')
logs_dir = os.path.join(out_dir, "logs")
database_dir = config.get('database_dir', "/group/ctbrowngrp/sourmash-db/gtdb-rs207")

reps_only = config.get('representatives_only', True)

# requires taxonomy with representitive info
print('reading taxonomy')
gtdb_taxonomy=config.get('gtdb_taxonomy', 'conf/gtdb-rs207.taxonomy.with-repinfo.csv')
taxDF = pd.read_csv(gtdb_taxonomy)

# get fastapaths
#fa_info = pd.read_csv("/home/ntpierce/2021-rank-compare/gtdb-rs207.fromfile.csv") # cols: ident, name, genome_filename, protein_filename

# set ident as index for each
#fa_info.set_index("ident", inplace=True)

# get comparison groups
compare_rank = config.get('rank', "family")
alltax_at_rank = taxDF[compare_rank].unique().tolist()

#taxDF.set_index("ident", inplace=True)
all_comparisons = []
ranktaxacc = []
grepD = {}
ranktaxD = defaultdict(list)
for ranktax in alltax_at_rank:
    # get query_idents
    subDF = taxDF[(taxDF[compare_rank] == ranktax)]
    all_accs = subDF["ident"].to_list()
    num_acc = len(all_accs)
    # testing:let's start smaller: ignore >50 comparisons
    if 2 <= num_acc <= 200: # can't compare a single genome :)
        rt = ranktax.replace(' ', '_') # no spaces plsss
        grepD[rt] = [ranktax]
        ranktaxD[rt] = all_accs
        these_rta = expand(f"{rt}/{{acc}}", acc = all_accs)
        ranktaxacc.extend(these_rta)

num_rt_to_compare = len(ranktaxD.keys())
print(f"found {num_rt_to_compare} comparison groups")

onstart:
    print("\n--- Starting ANI workflow ---\n")

onsuccess:
    print("\n--- Workflow executed successfully! ---\n")

onerror:
    print("Alas!\n")


rule all:
    input: 
        expand(f"{out_dir}/prefetch/{compare_rank}/gtdb-rs207.nucleotide-k{{ksize}}-scaled{{scaled}}.ani.csv", ksize = [21], scaled=[1000]) # 31,51
        #os.path.join(out_dir, "sourmash", f"{compare_rank}.sourmash-k{ksize}-scaled1000.ANI.csv.gz"),
        #expand(f"{out_dir}/prefetch/{compare_rank}/{{rta}}.nucleotide-k{{ksize}}-scaled1000.prefetch.csv", rta=ranktaxacc, ksize = [21,31,51]), #run: 21,31,51
        #expand(f"{out_dir}/prefetch/{compare_rank}/{{rank_tax}}.picklist.csv",rank_tax = ranktaxD.keys()),


rule write_picklist:
    input:
        taxonomy = gtdb_taxonomy,
    output:
        rt_picklist = f"{out_dir}/prefetch/{compare_rank}/{{rank_tax}}.picklist.csv",
    params:
        rt = lambda w: grepD[w.rank_tax],
    log: f"{logs_dir}/prefetch/{compare_rank}/{{rank_tax}}.picklist.log"
    benchmark: f"{logs_dir}/prefetch/{compare_rank}/{{rank_tax}}.picklist.benchmark",
    conda: "conf/envs/sourmash-4.4.yml"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 3000,
        time=30,
        runtime=60,
        partition="low2",
    shell:
        """
        echo "taxonomy is {input.taxonomy}"
        echo "taxonomy is {input.taxonomy}" > {log}

        # build picklist
        head -1 {input.taxonomy} > {output.rt_picklist} # grab header
        grep {params.rt} {input.taxonomy} >> {output.rt_picklist}
        cat {output.rt_picklist} | wc -l
        touch {output.rt_picklist}
        """

rule nucl_rank_prefetch:
    input:
        db = f"{database_dir}/gtdb-rs207.genomic.k{{ksize}}.zip", # scaled 1000
        rt_picklist = f"{out_dir}/prefetch/{compare_rank}/{{rank_tax}}.picklist.csv",
        taxonomy = gtdb_taxonomy,
    output: 
        prefetch=f"{out_dir}/prefetch/{compare_rank}/{{rank_tax}}/{{acc}}.nucleotide-k{{ksize}}-scaled{{scaled}}.prefetch.csv",
    params:
        alpha= "--dna",
        threshold_bp=10000,
        rt = lambda w: grepD[w.rank_tax],
        tmp_db=f"{{rank_tax}}.nucleotide-k{{ksize}}-scaled{{scaled}}.zip",
    log: f"{logs_dir}/prefetch/{compare_rank}/{{rank_tax}}/{{acc}}.nucleotide-k{{ksize}}-scaled{{scaled}}.prefetch.log"
    benchmark: f"{logs_dir}/prefetch/{compare_rank}/{{rank_tax}}/{{acc}}.nucleotide-k{{ksize}}-scaled{{scaled}}.prefetch.benchmark",
    conda: "conf/envs/sourmash-4.4.yml"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 3000,
        time=240,
        runtime=30,
        partition="low2",
    shell:
        """
        echo "DB is {input.db}"
        echo "DB is {input.db}" > {log}

        sourmash sig grep {wildcards.acc} {input.db} \
                 --ksize {wildcards.ksize} | sourmash prefetch - {input.db} \
                 -o {output.prefetch} -k {wildcards.ksize} {params.alpha} \
                 --picklist {input.rt_picklist}:ident:identprefix \
                 --threshold-bp={params.threshold_bp} --scaled {wildcards.scaled} 2>> {log}
        touch {output} # touch to prevent err with empty prefetch results
        """

rule write_prefetch_filelist:
    input: lambda w: expand(f"{out_dir}/prefetch/{compare_rank}/{{rta}}.nucleotide-k{{ksize}}-scaled{{scaled}}.prefetch.csv", rta=ranktaxacc, ksize = w.ksize, scaled=w.scaled)
    output: f"{out_dir}/prefetch/{compare_rank}/gtdb-rs207.nucleotide-k{{ksize}}-scaled{{scaled}}.prefetch.filelist",
    run:
        with open(str(output), 'w') as outF:
            for inF in input:
                fn = os.path.abspath(str(inF))
                outF.write(f"{fn}\n")

rule prefetch_to_ani_csv:
    input: 
        pf_filelist = f"{out_dir}/prefetch/{compare_rank}/gtdb-rs207.nucleotide-k{{ksize}}-scaled{{scaled}}.prefetch.filelist",
        taxonomy = gtdb_taxonomy,
    output: f"{out_dir}/prefetch/{compare_rank}/gtdb-rs207.nucleotide-k{{ksize}}-scaled{{scaled}}.ani.csv",
    log: f"{logs_dir}/prefetch-to-ani-csv/gtdb-rs207.nucleotide-k{{ksize}}-scaled{{scaled}}.log"
    benchmark: f"{logs_dir}/prefetch-to-ani-csv/gtdb-rs207.nucleotide-k{{ksize}}-scaled{{scaled}}.benchmark"
    conda: "conf/envs/sourmash-4.4.yml"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 3000,
        time=30,
        partition="low2",
    shell:
        """
        python prefetch-to-ani-csv.py --from-file {input.pf_filelist} --taxonomy {input.taxonomy} --recalculate-ani -o {output} 2> {log}
        """

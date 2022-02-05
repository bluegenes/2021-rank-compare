"""
Author: Tessa Pierce Ward, UC Davis Lab for Data Intensive Biology
Run: snakemake -j1 -n
"""

import pandas as pd
from collections import defaultdict
#from tqdm import tqdm
#tqdm.pandas()

gtdb_taxonomy=config.get('gtdb_taxonomy', '/group/ctbrowngrp/gtdb/gtdb-rs202.taxonomy.v2.with-strain.csv')

out_dir = config.get('output_dir', 'output.prefetch_compare')
logs_dir = os.path.join(out_dir, "logs")
basename = config.get('basename', 'gtdb-rs202')
database_dir = '/group/ctbrowngrp/gtdb/databases'

print('reading taxonomy')
taxDF = pd.read_csv(gtdb_taxonomy)

compare_rank = config.get('rank', "genus")
alltax_at_rank = taxDF[compare_rank].unique().tolist()
grouped_by_rank = taxDF.groupby(compare_rank)["ident"].apply(list)
picklist_dir = '/home/ntpierce/2021-rank-compare/gtdb-rs202-taxonomic-picklists'
# replace spaces with underscoare (e.g. for species names)
taxDF = taxDF.replace(r" ", "_", regex=True) # regex allows searching for partial str
taxDF.set_index('ident', inplace=True)

# check params are in the right format, build alpha-ksize combos
alpha_ksizes=[]
for alpha, info in config["alphabet_info"].items():
    scaled = info["scaled"]
    ksize = info["ksize"]
    if not isinstance(scaled, list):
        scaled = [scaled]
        config["alphabet_info"][alpha]["scaled"] = scaled
    if not isinstance(ksize, list):
        ksize=[ksize]
        config["alphabet_info"][alpha]["ksize"] = ksize
    # build a parameter for the right combinations
    #if alpha == 'protein':
    #    prot_alpha_ksizes += expand(f"{alpha}-k{{ksize}}", ksize = ksize)
    #if alpha == 'nucleotide':
    #    nucl_alpha_ksizes += expand(f"{alpha}-k{{ksize}}", ksize = ksize)
    alpha_ksizes += expand(f"{alpha}-k{{ksize}}", ksize = ksize)

rule all:
    input: 
        expand(f"{out_dir}/prefetch/{compare_rank}/{{rank_tax}}.{{ak}}.prefetch.csv", rank_tax = alltax_at_rank[:1], ak=alpha_ksizes)
         #expand(os.path.join(out_dir, 'prefetch', '{acc}.{alphak}.prefetch.csv'), acc=taxDF.index, alphak=alpha_ksizes),
         #expand(os.path.join(out_dir, 'prefetch', '{acc}.{alphak}.compare.csv'), acc=taxDF.index, alphak=alpha_ksizes),


# use prefetch to do comparison for each ident 
rule protein_rank_prefetch:
    input: 
        #db=f"{database_dir}/gtdb-rs202.{{alphabet}}.k{{ksize}}.zip", #
        db = f"{database_dir}/gtdb-rs202.protein.protein.k10.sbt.zip", # scaled 200
        picklist = f"{picklist_dir}/{compare_rank}/gtdb-rs202.{{rank_tax}}.csv"
    output: f"{out_dir}/prefetch/{compare_rank}/{{rank_tax}}/{{acc}}.protein-k{{ksize}}.prefetch.csv"
    params:
        alpha= "--protein",
        threshold_bp=30000,
        acc_identprefix = lambda w: w.acc.rsplit('.')[0]
    log: f"{logs_dir}/prefetch/{compare_rank}/{{rank_tax}}/{{acc}}.protein-k{{ksize}}.prefetch.log"
    benchmark: f"{logs_dir}/prefetch/{compare_rank}/{{rank_tax}}/{{acc}}.protein-k{{ksize}}.prefetch.benchmark",
    conda: "conf/envs/sourmash-dist-dbextract.yml"
    shell:
        """
        echo "DB is {input.db}"
        echo "DB is {input.db}" > {log}

        sourmash database extract {input.db} --identprefix {params.acc_identprefix} \
                 --ksize {wildcards.ksize} | sourmash prefetch - {input.db} \
                 --picklist {input.picklist}:ident:identprefix \
                 -o {output} -k {wildcards.ksize} {params.alpha} \
                 --threshold-bp={params.threshold_bp} 2>> {log}
        touch {output}
        """

rule nucl_rank_prefetch:
    input: 
        db = f"{database_dir}/ctb/gtdb-rs202.genomic.k{{ksize}}.sbt.zip", # scaled 1000
        picklist = f"{picklist_dir}/{compare_rank}/gtdb-rs202.{{rank_tax}}.csv"
    output: f"{out_dir}/prefetch/{compare_rank}/{{rank_tax}}/{{acc}}.nucleotide-k{{ksize}}.prefetch.csv"
    params:
        alpha= "--dna",
        threshold_bp=100000,
        acc_identprefix = lambda w: w.acc.rsplit('.')[0]
    log: f"{logs_dir}/prefetch/{compare_rank}/{{rank_tax}}/{{acc}}.nucleotide-k{{ksize}}.prefetch.log"
    benchmark: f"{logs_dir}/prefetch/{compare_rank}/{{rank_tax}}/{{acc}}.nucleotide-k{{ksize}}.prefetch.benchmark",
    conda: "conf/envs/sourmash-dist-dbextract.yml"
    shell:
        """
        echo "DB is {input.db}"
        echo "DB is {input.db}" > {log}

        sourmash database extract {input.db} --identprefix {params.acc_identprefix} \
                 --ksize {wildcards.ksize} | sourmash prefetch - {input.db} \
                 --picklist {input.picklist}:ident:identprefix \
                 -o {output} -k {wildcards.ksize} {params.alpha} \
                 --threshold-bp={params.threshold_bp} 2>> {log}
        touch {output}
        """

rule aggregate_nucl_prefetch:
    input: lambda w: expand(f"{out_dir}/prefetch/{compare_rank}/{{rt}}/{{acc}}.nucleotide-k{{ksize}}.prefetch.csv", acc = grouped_by_rank[w.rank_tax], rt = w.rank_tax, ksize = w.ksize)
    output: f"{out_dir}/prefetch/{compare_rank}/{{rank_tax}}.nucleotide-k{{ksize}}.prefetch.csv"
    log: f"{logs_dir}/prefetch/{compare_rank}/{{rank_tax}}.nucleotide-k{{ksize}}.prefetch.log" 
    benchmark: f"{logs_dir}/prefetch/{compare_rank}/{{rank_tax}}.nucleotide-k{{ksize}}.prefetch.benchmark" 
    shell:
        """
        cat {input} > {output} 2> {log}
        """

rule aggregate_prot_prefetch:
    input: lambda w: expand(f"{out_dir}/prefetch/{compare_rank}/{{rt}}/{{acc}}.protein-k{{ksize}}.prefetch.csv", acc = grouped_by_rank[w.rank_tax], rt = w.rank_tax, ksize = w.ksize)
    output: f"{out_dir}/prefetch/{compare_rank}/{{rank_tax}}.protein-k{{ksize}}.prefetch.csv"
    log: f"{logs_dir}/prefetch/{compare_rank}/{{rank_tax}}.protein-k{{ksize}}.prefetch.log" 
    benchmark: f"{logs_dir}/prefetch/{compare_rank}/{{rank_tax}}.protein-k{{ksize}}.prefetch.benchmark" 
    shell:
        """
        cat {input} > {output} 2> {log}
        """

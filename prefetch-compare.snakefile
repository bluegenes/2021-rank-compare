"""
Author: Tessa Pierce Ward, UC Davis Lab for Data Intensive Biology
Run: snakemake -j1 -n
"""

import pandas as pd
from collections import defaultdict
#from tqdm import tqdm
#tqdm.pandas()

# requires taxonomy with representitive info
gtdb_taxonomy=config.get('gtdb_taxonomy', 'gtdb-rs202.taxonomy.with-repinfo.csv')

out_dir = config.get('output_dir', 'output.prefetch_compare')
logs_dir = os.path.join(out_dir, "logs")
basename = config.get('basename', 'gtdb-rs202')
database_dir = '/group/ctbrowngrp/gtdb/databases'
reps_only = config.get('representatives_only', True)

print('reading taxonomy')
taxDF = pd.read_csv(gtdb_taxonomy)

compare_rank = config.get('rank', "genus")

alltax_at_rank, accs_to_prefetch = [],[]
if compare_rank.lower() == "all":
    if reps_only:
        accs_to_prefetch = taxDF[taxDF['is_representative'] == True]["ident"].tolist()
    else:
        accs_to_prefetch = taxDF["ident"].tolist()
else:
    # here, get representatives for the rank
    if reps_only:
        accs_to_prefetch = taxDF[taxDF['is_representative'] == True].groupby(compare_rank)["ident"].apply(list)
    else:
        accs_to_prefetch = taxDF.groupby(compare_rank)["ident"].apply(list)
    alltax_at_rank = taxDF[compare_rank].unique().tolist()
# group all accs at rank (useful for all x all comparisons)
#grouped_by_rank = taxDF.groupby(compare_rank)["ident"].apply(list)
    
picklist_dir = '/home/ntpierce/2021-rank-compare/gtdb-rs202-taxonomic-picklists'
# replace spaces with underscoare (e.g. for species names)
taxDF = taxDF.replace(r" ", "_", regex=True) # regex allows searching for partial str
taxDF.set_index('ident', inplace=True)
    
ranktaxacc=[]
if compare_rank != "all":
    for ranktax in alltax_at_rank:
        if reps_only:
            accs = accs_to_prefetch[ranktax][0],
        else:
            # does this work? Or stil need [0]?
            accs = accs_to_prefetch[ranktax],
        these_ranktaxacc = expand(f"{ranktax}/{{acc}}", acc=accs)
        ranktaxacc.extend(these_ranktaxacc)

# check params are in the right format, build alpha-ksize combos
alpha_ksize_scaled=[]
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
    alpha_ksize_scaled += expand(f"{alpha}-k{{ksize}}-scaled{{scaled}}", ksize = ksize, scaled=scaled)

wildcard_constraints:
    rank_tax="\w+",
    ksize="\w+"

outF = []
if compare_rank == "all":
    outF=expand(f"{out_dir}/prefetch/{basename}.{{aks}}.prefetch.csv.gz", aks=alpha_ksize_scaled)
    outF += expand(f"{out_dir}/prefetch/{basename}.{{aks}}.binned_containment_ani.csv", aks=alpha_ksize_scaled)
else:
    outF=expand(f"{out_dir}/prefetch/{basename}.{compare_rank}.{{aks}}.prefetch.csv", aks=alpha_ksize_scaled)


rule all:
    input: 
        outF
        #expand(f"{out_dir}/prefetch/{compare_rank}/{{rta}}.{{aks}}.prefetch.csv", rta = ranktaxacc, aks=alpha_ksize_scaled)
        #expand(f"{out_dir}/prefetch/{basename}.{compare_rank}.{{aks}}.prefetch.csv", aks=alpha_ksize_scaled)
         #expand(os.path.join(out_dir, 'prefetch', '{acc}.{alphak}.prefetch.csv'), acc=taxDF.index, alphak=alpha_ksizes),
         #expand(os.path.join(out_dir, 'prefetch', '{acc}.{alphak}.compare.csv'), acc=taxDF.index, alphak=alpha_ksizes),


# use prefetch to do comparison for each ident 
rule protein_rank_prefetch:
    input: 
        #db=f"{database_dir}/gtdb-rs202.{{alphabet}}.k{{ksize}}.zip", #
        db = f"{database_dir}/gtdb-rs202.protein.protein.k10.sbt.zip", # scaled 200, k 10
        picklist = f"{picklist_dir}/{compare_rank}/gtdb-rs202.{{rank_tax}}.csv"
    output: f"{out_dir}/prefetch/{compare_rank}/{{rank_tax}}/{{acc}}.protein-k{{ksize}}-scaled{{scaled}}.prefetch.csv"
    params:
        alpha= "--protein",
        threshold_bp=3000,
        acc_identprefix = lambda w: w.acc.rsplit('.')[0]
    log: f"{logs_dir}/prefetch/{compare_rank}/{{rank_tax}}/{{acc}}.protein-k{{ksize}}-scaled{{scaled}}.prefetch.log"
    benchmark: f"{logs_dir}/prefetch/{compare_rank}/{{rank_tax}}/{{acc}}.protein-k{{ksize}}-scaled{{scaled}}.prefetch.benchmark",
    conda: "conf/envs/sourmash-dist-dbextract.yml"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 6000,
        runtime=6000
    shell:
        """
        echo "DB is {input.db}"
        echo "DB is {input.db}" > {log}

        sourmash database extract {input.db} --identprefix {params.acc_identprefix} \
                 --ksize {wildcards.ksize} | sourmash prefetch - {input.db} \
                 --picklist {input.picklist}:ident:identprefix \
                 -o {output} -k {wildcards.ksize} {params.alpha} \
                 --threshold-bp={params.threshold_bp}  --scaled {wildcards.scaled} 2>> {log}
        touch {output}
        """

rule nucl_rank_prefetch:
    input: 
        db = f"{database_dir}/ctb/gtdb-rs202.genomic.k{{ksize}}.sbt.zip", # scaled 1000
        picklist = f"{picklist_dir}/{compare_rank}/gtdb-rs202.{{rank_tax}}.csv"
    output: f"{out_dir}/prefetch/{compare_rank}/{{rank_tax}}/{{acc}}.nucleotide-k{{ksize}}-scaled{{scaled}}.prefetch.csv"
    params:
        alpha= "--dna",
        threshold_bp=10000,
        acc_identprefix = lambda w: w.acc.rsplit('.')[0]
    log: f"{logs_dir}/prefetch/{compare_rank}/{{rank_tax}}/{{acc}}.nucleotide-k{{ksize}}-scaled{{scaled}}.prefetch.log"
    benchmark: f"{logs_dir}/prefetch/{compare_rank}/{{rank_tax}}/{{acc}}.nucleotide-k{{ksize}}-scaled{{scaled}}.prefetch.benchmark",
    conda: "conf/envs/sourmash-dist-dbextract.yml"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 6000,
        runtime=6000
    shell:
        """
        echo "DB is {input.db}"
        echo "DB is {input.db}" > {log}

        sourmash database extract {input.db} --identprefix {params.acc_identprefix} \
                 --ksize {wildcards.ksize} | sourmash prefetch - {input.db} \
                 --picklist {input.picklist}:ident:identprefix \
                 -o {output} -k {wildcards.ksize} {params.alpha} \
                 --threshold-bp={params.threshold_bp} --scaled {wildcards.scaled} 2>> {log}
        touch {output}
        """

rule aggregate_rank_nucl_prefetch:
    input: lambda w: expand(f"{out_dir}/prefetch/{compare_rank}/{{rta}}.nucleotide-k{{ksize}}-scaled{{scaled}}.prefetch.csv",  rta = ranktaxacc, ksize = w.ksize, scaled=w.scaled)
    output: f"{out_dir}/prefetch/{basename}.{compare_rank}.nucleotide-k{{ksize}}-scaled{{scaled}}.prefetch.csv"
    log: f"{logs_dir}/prefetch/{basename}.{compare_rank}.nucleotide-k{{ksize}}-scaled{{scaled}}.prefetch.log" 
    benchmark: f"{logs_dir}/prefetch/{basename}.{compare_rank}.nucleotide-k{{ksize}}-scaled{{scaled}}.prefetch.benchmark" 
    run:
        # aggregate all csvs --> single csv
        aggDF = pd.concat([pd.read_csv(str(csv), sep=",") for csv in input])
        aggDF["alpha-ksize"] =  "nucleotide-k" + str(wildcards.ksize)
        aggDF.to_csv(str(output), index=False)

rule aggregate_rank_prot_prefetch:
    input: lambda w: expand(f"{out_dir}/prefetch/{compare_rank}/{{rta}}.protein-k{{ksize}}-scaled{{scaled}}.prefetch.csv",  rta = ranktaxacc, ksize = w.ksize, scaled=w.scaled)
    output: f"{out_dir}/prefetch/{basename}.{compare_rank}.protein-k{{ksize}}-scaled{{scaled}}.prefetch.csv"
    log: f"{logs_dir}/prefetch/{basename}.{compare_rank}.protein-k{{ksize}}-scaled{{scaled}}.prefetch.log" 
    benchmark: f"{logs_dir}/prefetch/{basename}.{compare_rank}.protein-k{{ksize}}-scaled{{scaled}}.prefetch.benchmark" 
    run:
        # aggregate all csvs --> single csv
        aggDF = pd.concat([pd.read_csv(str(csv), sep=",") for csv in input])
        aggDF["alpha-ksize"] =  "protein-k" + str(wildcards.ksize)
        aggDF.to_csv(str(output), index=False)


# use prefetch to do comparison for each ident 
rule protein_all_prefetch:
    input: 
        #db=f"{database_dir}/gtdb-rs202.{{alphabet}}.k{{ksize}}.zip", #
        db = f"{database_dir}/gtdb-rs202.protein.protein.k10.sbt.zip", # scaled 200, k 10
        #picklist = gtdb_taxonomy
    output: f"{out_dir}/prefetch/gtdb-all/{{acc}}.protein-k{{ksize}}-scaled{{scaled}}.prefetch.csv"
    params:
        alpha= "--protein",
        threshold_bp=30000,
        acc_identprefix = lambda w: w.acc.rsplit('.')[0]
    log: f"{logs_dir}/prefetch/gtdb-all/{{acc}}.protein-k{{ksize}}-scaled{{scaled}}.prefetch.log"
    benchmark: f"{logs_dir}/prefetch/gtdb-all/{{acc}}.protein-k{{ksize}}-scaled{{scaled}}.prefetch.benchmark",
    conda: "conf/envs/sourmash-dist-dbextract.yml"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 15000,
        runtime=130
    shell:
        """
        echo "DB is {input.db}"
        echo "DB is {input.db}" > {log}

        sourmash database extract {input.db} --identprefix {params.acc_identprefix} \
                 --ksize {wildcards.ksize} | sourmash prefetch - {input.db} \
                 -o {output} -k {wildcards.ksize} {params.alpha} \
                 --threshold-bp={params.threshold_bp} --scaled {wildcards.scaled} > {log} 2>&1
        touch {output}
        """
                 #--picklist {input.picklist}:ident:identprefix \

# use prefetch to do comparison for each ident 
rule nucl_all_prefetch:
    input: 
        db = f"{database_dir}/ctb/gtdb-rs202.genomic.k{{ksize}}.sbt.zip", #scaled 1000
    output: f"{out_dir}/prefetch/gtdb-all/{{acc}}.nucleotide-k{{ksize}}-scaled{{scaled}}.prefetch.csv"
    params:
        alpha= "--dna",
        threshold_bp=10000,
        acc_identprefix = lambda w: w.acc.rsplit('.')[0]
    log: f"{logs_dir}/prefetch/gtdb-all/{{acc}}.nucleotide-k{{ksize}}-scaled{{scaled}}.prefetch.log"
    benchmark: f"{logs_dir}/prefetch/gtdb-all/{{acc}}.nucleotide-k{{ksize}}-scaled{{scaled}}.prefetch.benchmark",
    conda: "conf/envs/sourmash-dist-dbextract.yml"
    threads: 1
    resources:
        #mem_mb=lambda wildcards, attempt: attempt * 20000,
        mem_mb=lambda wildcards, attempt: attempt * 6000,
        disk_mb=lambda wildcards, attempt: attempt * 6000,
        runtime=120,
        time=120,
        partition="low2",
        #partition="med2",
    shell:
        """
        echo "DB is {input.db}"
        echo "DB is {input.db}" > {log}

        sourmash database extract {input.db} --identprefix {params.acc_identprefix} \
                 --ksize {wildcards.ksize} | sourmash prefetch - {input.db} \
                 -o {output} -k {wildcards.ksize} {params.alpha} \
                 --threshold-bp={params.threshold_bp} --scaled {wildcards.scaled} > {log} 2>&1
        touch {output}
        """
                 #--picklist {input.picklist}:ident:identprefix \

localrules: aggregate_allgtdb_prot_prefetch
rule aggregate_allgtdb_prot_prefetch:
    input: lambda w: expand(f"{out_dir}/prefetch/gtdb-all/{{acc}}.protein-k{{ksize}}-scaled{{scaled}}.prefetch.csv",  acc = accs_to_prefetch, ksize = w.ksize, scaled=w.scaled)
    output: f"{out_dir}/prefetch/{basename}.protein-k{{ksize}}-scaled{{scaled}}.prefetch.csv.gz"
    log: f"{logs_dir}/prefetch/{basename}.protein-k{{ksize}}-scaled{{scaled}}.prefetch.log" 
    benchmark: f"{logs_dir}/prefetch/{basename}.protein-k{{ksize}}-scaled{{scaled}}.prefetch.benchmark" 
    run:
        # aggregate all csvs --> single csv
        dfs_to_concat = []
        for inF in input:
            this_df = pd.read_csv(str(inF), sep=',')
            this_df = this_df[~(this_df["match_name"] == this_df["query_name"])]
            if not this_df.empty:
                this_df["alpha-ksize"] =  "protein-k" + str(wildcards.ksize)
                dfs_to_concat.append(this_df)

        aggDF = pd.concat(dfs_to_concat)
        # aggregate all csvs --> single csv
        #aggDF = pd.concat([pd.read_csv(str(csv), sep=",") for csv in input])
        #aggDF["alpha-ksize"] =  "protein-k" + str(wildcards.ksize)
        aggDF.to_csv(str(output), index=False)

localrules: aggregate_allgtdb_nucl_prefetch
rule aggregate_allgtdb_nucl_prefetch:
    input: lambda w: expand(f"{out_dir}/prefetch/gtdb-all/{{acc}}.nucleotide-k{{ksize}}-scaled{{scaled}}.prefetch.csv",  acc = accs_to_prefetch, ksize = w.ksize, scaled=w.scaled)
    output: 
        full_output=f"{out_dir}/prefetch/{basename}.nucleotide-k{{ksize}}-scaled{{scaled}}.prefetch.csv.gz",
        binned_containment_ani=f"{out_dir}/prefetch/{basename}.nucleotide-k{{ksize}}-scaled{{scaled}}.binned_containment_ani.csv"
    log: f"{logs_dir}/prefetch/{basename}.nucleotide-k{{ksize}}-scaled{{scaled}}.prefetch.log" 
    benchmark: f"{logs_dir}/prefetch/{basename}.nucleotide-k{{ksize}}-scaled{{scaled}}.prefetch.benchmark" 
    run:
        # aggregate all csvs --> single csv
        dfs_to_concat = []
        ksize = str(wildcards.ksize)
        for inF in input:
            this_df = pd.read_csv(str(inF), sep=',')
            this_df = this_df[~(this_df["match_name"] == this_df["query_name"])]
            if not this_df.empty:
                #this_df['check_dupes'] = this_df[["match_name", "query_name"]].apply(set, axis=1)
                this_df["alpha-ksize"] =  f"nucleotide-{ksize}"
                dfs_to_concat.append(this_df)

        # aggregate all
        aggDF = pd.concat(dfs_to_concat)
        print(aggDF.shape)
        # drop duplicated comparisons
        #aggDF.drop_duplicates(subset = "check_dupes", inplace=True)
        print(aggDF.shape)
        # get average containment/ani
        aggDF['avg_containment'] = aggDF[["f_query_match", "f_match_query"]].mean(axis=1)
        aggDF['avg_containment_ani'] = aggDF[["query_ani", "match_ani"]].mean(axis=1)
        # bin by containment; compare ANI
        containment_bins05 = np.arange(0, 1.01, 0.05)
        aggDF['binned_avg_containment05'] = pd.cut(x=aggDF['avg_containment'], bins=containment_bins05)
        contain05_info = aggDF.groupby('binned_avg_containment05').agg({'avg_containment_ani': ['count','min', 'max', 'mean']}).round(2)
        contain05_info.columns = contain05_info.columns.droplevel()
        contain05_info.reset_index(inplace=True)
        rename_cols = {"min": "minANI", "max": "maxANI", "mean": "meanANI"}
        contain05_info.rename(columns = rename_cols, inplace=True)
        # make separate col with str so we can modify it. keep categorical col for plotting.
        contain05_info['crange'] = contain05_info['binned_avg_containment05'].astype("str")
        contain05_info['highContainment'] = contain05_info['crange'].str.split(',', expand=True)[1].str.rstrip(']').str.strip()
        contain05_info.drop(columns=["crange"], inplace=True) # no need for this dup col. keep categorical col for plotting.
        contain05_info.to_csv(output.binned_containment_ani, index=False)
        # todo: maybe plot ANI?
        aggDF.to_csv(str(output.full_output), index=False)
        
        #aggDF = pd.concat([pd.read_csv(str(csv), sep=",") for csv in input])
        #aggDF["alpha-ksize"] =  "protein-k" + str(wildcards.ksize)
        #aggDF[~aggDF["match_name"] == aggDF["query_name"]]
        #aggDF.to_csv(str(output), index=False)


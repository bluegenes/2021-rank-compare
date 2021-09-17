import csv
import pandas as pd

basename = config.get('basename', 'gtdb-rs202')

database_dir = config.get("database_dir", "/group/ctbrowngrp/gtdb/databases")
out_dir = config.get('output_dir', 'output.rank-compare')
logs_dir = os.path.join(out_dir, 'logs')

alpha_ksizes,aks = [],[]
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
    #config["alphabet_info"][alpha]["select_params"] = expand("{alpha}-k{ksize}-scaled{scaled}", alpha=alpha, scaled=scaled, ksize=ksize)
    #these_alpha_ksize = expand("{alpha}-k{{ksize}}", ksize = ksize)
    alpha_ksizes += expand(f"{alpha}-k{{ksize}}", ksize = ksize)
    aks += expand(f"{alpha}-k{{ksize}}-scaled{{scaled}}", ksize = ksize, scaled=scaled)

# some easier vars
alphabet_info = config['alphabet_info']
#prot_ksizes = alphabet_info['protein']['ksize']
#nucl_ksizes = alphabet_info['nucleotide']['ksize']

lca_scaledD = {"protein": 1000, "nucleotide": 10000}

print(alpha_ksizes)

rule all:
    input: expand(f"{out_dir}/lca-unique-kmers/{basename}.{{alphak}}.lca-unique-kmers.txt",alphak=alpha_ksizes),


rule unique_kmers_lca:
    input: f"{database_dir}/lca/gtdb-rs202.{{alphabet}}.k{{ksize}}.lca.json.gz"
    #output: f"{out_dir}/lca-unique-kmers/{basename}.{{alphabet}}-k{{ksize}}-scaled{{scaled}}.lca-unique-kmers.txt",
    output: f"{out_dir}/lca-unique-kmers/{{basename}}.{{alphabet}}-k{{ksize}}.lca-unique-kmers.txt",
    params:
        lca_scaled = lambda w: lca_scaledD[w.alphabet]
    log: os.path.join(logs_dir, "lca-unique-kmers", "{basename}.{{alphabet}}-k{{ksize}}.log" )
    benchmark: os.path.join(logs_dir, "lca-unique-kmers", "{basename}.{{alphabet}}-k{{ksize}}.benchmark" )
    shell:
        """
        python unique-kmers-by-lineage.py --db {input} \
               --ksize {wildcards.ksize} --scaled {params.lca_scaled} \
               --output-csv {output} 2> {log}
        """

#rule intersect_shared_kmers_picklist:
#    input: "{database_dir}/gtdb-rs202.{alphabet}.k{ksize}.zip"
#    output: "{}"


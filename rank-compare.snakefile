"""
Author: Tessa Pierce Ward, UC Davis Lab for Data Intensive Biology
Run: snakemake -j1 -n
"""

import pandas as pd
from collections import defaultdict
from sourmash.tax import tax_utils

gtdb_taxonomy=config.get('gtdb_taxonomy', '/group/ctbrowngrp/gtdb/gtdb-rs202.taxonomy.v2.with-strain.csv')

#sample_info = pd.read_csv(config['sample_info'], sep = '\t')
#sample_column = config.get('sample_colname', 'Genome')
#lineage_column = config.get('sample_colname', 'Lineage')

#SAMPLES = sample_info[sample_column].to_list()
#data_dir = config['data_dir']
out_dir = config.get('output_dir', 'outputs')
logs_dir = os.path.join(out_dir, "logs")
basename = config.get('basename', 'gtdb-rs202')
database_dir = '/group/ctbrowngrp/gtdb/databases'
# how to incorporate db info??
search_databases = config.get('search_databases', ['gtdb-rs202'])
input_type = config.get('input_type', 'protein')

fasta_type= {'dna': "genomic", 'protein': 'protein', 'dayhoff': 'protein', 'hp': 'protein'}

taxDF = pd.read_csv(gtdb_taxonomy)
# replace spaces with underscoare
taxDF = taxDF.replace(r" ", "_", regex=True) # regex allows searching for partial str
# temp hack for doing a subset
#taxDF = taxDF[taxDF['phylum'] == 'p__Thermotogota']

taxonomyD = {}
rank_taxinfo = []

# make list of output bloom filter basenames
for rank in tax_utils.ascending_taxlist(include_strain=False):
    if rank == "phylum":
        rank_taxonomies = taxDF[rank].unique() # unique taxa at each rank
        #rank_taxonomies = [x.replace(' ', '_') for x in rank_taxonomies]
        taxonomyD[rank] = rank_taxonomies
        rank_tax = expand(os.path.join(rank, f'{basename}.{{taxon}}'), taxon = rank_taxonomies) 
        rank_taxinfo.extend(rank_tax)

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
    #if alpha in ['nucleotide', 'protein']:
    alpha_ksizes += expand(f"{alpha}-k{{ksize}}", ksize = ksize)

rule all:
    input: 
         expand(os.path.join(out_dir, 'unroll-compare', '{rt}.{alphak}.compare.csv'), rt=rank_taxinfo, alphak=alpha_ksizes),
#         expand(os.path.join(out_dir, 'unroll-compare', f'{basename}.{{alphak}}.compare.csv'), alphak=alpha_ksizes),
         #expand(os.path.join(out_dir, f'{basename}-taxonomic-picklists', '{rt}.csv'), rt=rank_taxinfo),
         #expand(os.path.join(out_dir, 'shared-kmers', '{ranktaxinf}.{alphak}.shared-kmers.csv'), ranktaxinf=rank_taxinfo, alphak=alpha_ksizes),


rule make_picklists:
    input: gtdb_taxonomy
    output: 
        expand(os.path.join(out_dir, '{{basename}}-taxonomic-picklists', '{rt}.csv'), rt=rank_taxinfo)
    params:
        output_dir = os.path.join(out_dir, '{basename}-taxonomic-picklists')
    log: os.path.join(logs_dir, "make_picklists", "{basename}.make_picklists.log")
    benchmark: os.path.join(logs_dir, "make_picklists", "{basename}.make_picklists.benchmark")
    shell:
        """
        python make-gtdb-rank-picklists.py --output-base {wildcards.basename} \
                                              --output-dir {params.output_dir} \
                                              --gtdb-taxonomy {input} 2> {log}
        """


# use sourmash sig intersect to get k-mers shared by taxonomy group
# this is k-mers shared by ALL members of a group = drops to zero v quickly. What we're looking for
# is pairwise jaccard/containment, I think
##rule shared_kmers_picklist_intersect:
#    input: 
#        db=f"{database_dir}/gtdb-rs202.{{alphabet}}.k{{ksize}}.zip",
#        picklist=f"{out_dir}/{{basename}}-taxonomic-picklists/{{rank}}/{{basename}}.{{taxon}}.csv",
#    output: f"{out_dir}/shared-kmers/{{rank}}/{{basename}}.{{taxon}}.{{alphabet}}-k{{ksize}}.shared-kmers.sig"
#    params:
#        alpha= lambda w: f"--{w.alphabet}",
#    log: f"{logs_dir}/shared-kmers/{{rank}}/{{basename}}.{{taxon}}.{{alphabet}}-k{{ksize}}.shared-kmers.log"
#    benchmark: f"{logs_dir}/shared-kmers/{{rank}}/{{basename}}.{{taxon}}.{{alphabet}}-k{{ksize}}.shared-kmers.benchmark",
#    shell:
#        """
#        sourmash sig intersect {input.db} \
#                 --ksize {wildcards.ksize} {params.alpha} \
#                 --picklist {input.picklist}:ident:identprefix -o {output} 2> {log}
#        """

# use sourmash sig intersect to get k-mers shared by taxonomy group
#rule describe_intersected_sig: 
#    input: rules.shared_kmers_picklist_intersect.output,
#    output: f"{out_dir}/shared-kmers/{{rank}}/{{basename}}.{{taxon}}.{{alphabet}}-k{{ksize}}.shared-kmers.csv"
#    params:
#        alpha= lambda w: f"--{w.alphabet}"
#    log: f"{logs_dir}/shared-kmers/{{rank}}/{{basename}}.{{taxon}}.{{alphabet}}-k{{ksize}}.describe.log"
#    benchmark: f"{logs_dir}/shared-kmers/{{rank}}/{{basename}}.{{taxon}}.{{alphabet}}-k{{ksize}}.describe.benchmark"
#    shell:
#        """
#        sourmash sig describe {input} --csv {output} \
#                 --ksize {wildcards.ksize} {params.alpha} 2> {log}
#        """


# compare using scaled = 200; unroll
rule jaccard_compare:
    input:
        picklist= os.path.join(out_dir, 'gtdb-rs202-taxonomic-picklists', '{rank}', '{basename}.{taxon}.csv'),
        zipfile=lambda w: expand(os.path.join(database_dir, 'gtdb-rs202.{fasta_type}.{{alphabet}}-k{{ksize}}.zip'), fasta_type = fasta_type[w.alphabet]),
    output:
        jaccard= os.path.join(out_dir, 'compare', '{rank}', '{basename}.{taxon}.{alphabet}-k{ksize}.jaccard.csv')
    log: os.path.join(logs_dir, 'compare', '{rank}', '{basename}.{taxon}.{alphabet}-k{ksize}.jaccard.log')
    benchmark: os.path.join(logs_dir, 'compare', '{rank}', '{basename}.{taxon}.{alphabet}-k{ksize}.jaccard.benchmark')
    shell:
        """
        sourmash compare --picklist {input.picklist}:ident:ident {input.zipfile} \
                             --csv {output.jaccard}  --{wildcards.alphabet} 2> {log}
        """


rule max_contain_compare:
    input:
        picklist= os.path.join(out_dir, 'gtdb-rs202-taxonomic-picklists', '{rank}', '{basename}.{taxon}.csv'),
        zipfile=lambda w: expand(os.path.join(database_dir, 'gtdb-rs202.{fasta_type}.{{alphabet}}-k{{ksize}}.zip'), fasta_type = fasta_type[w.alphabet]),
    output:
        max_contain= os.path.join(out_dir, 'compare', '{rank}', '{basename}.{taxon}.{alphabet}-k{ksize}.max_containment.csv')
    log: os.path.join(logs_dir, 'compare', '{rank}', '{basename}.{taxon}.{alphabet}-k{ksize}.max_contain.log')
    benchmark: os.path.join(logs_dir, 'compare', '{rank}', '{basename}.{taxon}.{alphabet}-k{ksize}.max_contain.benchmark')
    shell:
        """
        sourmash compare --picklist {input.picklist}:ident:ident {input.zipfile} \
                             --max-containment --csv {output.max_contain} 2>> {log}
        """


rule unroll_compare:
    input:
        max_contain= os.path.join(out_dir, 'compare', '{rank}', '{basename}.{taxon}.{alphabet}-k{ksize}.max_containment.csv'),
        jaccard= os.path.join(out_dir, 'compare', '{rank}', '{basename}.{taxon}.{alphabet}-k{ksize}.jaccard.csv')
    output:
        os.path.join(out_dir, 'unroll-compare', '{rank}', '{basename}.{taxon}.{alphabet}-k{ksize}.compare.csv')
    log: os.path.join(logs_dir, 'unroll', '{rank}', '{basename}.{taxon}.{alphabet}-k{ksize}.unroll-compare.log')
    benchmark: os.path.join(logs_dir, 'unroll', '{rank}', '{basename}.{taxon}.{alphabet}-k{ksize}.unroll-compare.benchmark')
    shell:
        """
        python  unroll-compare.py {input.max_contain} {input.jaccard} -o {output} 2> {log}
        """

#rule aggregate_unrolled:
#    input:
#        expand(os.path.join(out_dir, 'unroll-compare', '{rt}.{alphak}.compare.csv'), rt=rank_taxinfo, alphak=alpha_ksizes)
#        #os.path.join(out_dir, 'unroll-compare', '{rank}', '{basename}.{taxon}.{alphabet}-k{ksize}.compare.csv')
#    output:
#        os.path.join(out_dir, 'unroll-compare', '{basename}.{alphabet}-k{ksize}.compare.csv')
#    log: os.path.join(logs_dir, 'aggregate_unrolled', '{basename}.{alphabet}-k{ksize}.aggregate_unrolled.log')
#    benchmark: os.path.join(logs_dir, 'aggregate_unrolled', '{basename}.{alphabet}-k{ksize}.aggregate_unrolled.benchmark')
#    # could shell cat these, but not sure if we'd hit the cli limit..
#    #run:
#        with open(str(output), 'w') as outF):
#            for inF in input:
#                with open(str(inF), 'r') as f:
#                    for line in f:
#                        outF.write(line)
            
        

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
database_dir = 'databases'
# how to incorporate db info??
search_databases = config.get('search_databases', ['gtdb-rs202'])
input_type = config.get('input_type', 'protein')

fasta_type= {'dna': "genomic", 'protein': 'protein', 'dayhoff': 'protein', 'hp': 'protein'}

taxDF = pd.read_csv(gtdb_taxonomy)
taxonomyD = {}
rank_taxinfo = []

# make list of output bloom filter basenames
for rank in tax_utils.ascending_taxlist(include_strain=False):
    rank_taxonomies = taxDF[rank].unique() # unique taxa at each rank
    rank_taxonomies = [x.replace(' ', '_') for x in rank_taxonomies]
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
#        expand(os.path.join(out_dir, 'unroll-compare', '{rt}.{alphabet}-k{ksize}.csv'), rt=rank_taxinfo)
         expand(os.path.join(out_dir, f'{basename}-taxonomic-picklists', '{rt}.csv'), rt=rank_taxinfo),
         expand(os.path.join(out_dir, 'shared-kmers', '{ranktaxinf}.{alphak}.shared-kmers.csv'), ranktaxinf=rank_taxinfo, alphak=alpha_ksizes),
         #expand(os.path.join(out_dir, 'bloom-filters', '{ranktaxinf}.{alphak}.nodegraph'), ranktaxinf=rank_taxinfo, alphak=alpha_ksizes)


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
rule shared_kmers_picklist_intersect:
    input: 
        db=f"{database_dir}/gtdb-rs202.{{alphabet}}.k{{ksize}}.zip",
        picklist=f"{out_dir}/{{basename}}-taxonomic-picklists/{{rank}}/{{basename}}.{{taxon}}.csv",
    output: f"{out_dir}/shared-kmers/{{rank}}/{{basename}}.{{taxon}}.{{alphabet}}-k{{ksize}}.shared-kmers.sig"
    params:
        alpha= lambda w: f"--{w.alphabet}",
    log: f"{logs_dir}/shared-kmers/{{rank}}/{{basename}}.{{taxon}}.{{alphabet}}-k{{ksize}}.shared-kmers.log"
    benchmark: f"{logs_dir}/shared-kmers/{{rank}}/{{basename}}.{{taxon}}.{{alphabet}}-k{{ksize}}.shared-kmers.benchmark",
    shell:
        """
        sourmash sig intersect {input.db} \
                 --ksize {wildcards.ksize} {params.alpha} \
                 --picklist {input.picklist} -o {output} 2> {log}
        """

# use sourmash sig intersect to get k-mers shared by taxonomy group
rule describe_intersected_sig: 
    input: rules.shared_kmers_picklist_intersect.output,
    output: f"{out_dir}/shared-kmers/{{rank}}/{{basename}}.{{taxon}}.{{alphabet}}-k{{ksize}}.shared-kmers.csv"
    params:
        alpha= lambda w: f"--{w.alphabet}"
    log: f"{logs_dir}/shared-kmers/{{rank}}/{{basename}}.{{taxon}}.{{alphabet}}-k{{ksize}}.describe.log"
    benchmark: f"{logs_dir}/shared-kmers/{{rank}}/{{basename}}.{{taxon}}.{{alphabet}}-k{{ksize}}.describe.benchmark"
    shell:
        """
        sourmash sig describe {input} --csv {output} \
                 --ksize {wildcards.ksize} {params.alpha} 2> {log}
        """


localrules: write_config_species_bloom_filters
rule write_config_species_bloom_filters:
    input: 
        picklist= os.path.join(out_dir, 'gtdb-rs202-taxonomic-picklists', 'species', '{basename}.{taxon}.csv'),
    output: os.path.join(out_dir, 'conf', '{basename}.{taxon}.bf.yml'),
    threads: 1
    run:
       # location for all generated files
       inF = str(input)
       with open(str(output), 'w') as outF:
           outF.write(f"basename: {wildcards.taxon}\n")
           outF.write(f"output_dir: genbank")
           outF.write(f"taxonomy_csv: {inF}\n")
           outF.write("run_alphabets:\n")
           #outF.write("  - nucleotide\n")
           outF.write("  - protein\n")
           #outF.write("  - dayhoff\n")
           outF.write("alphabet_info:\n")
           outF.write("  nucleotide:\n")
           outF.write("    ksize: [21,31]\n")
           outF.write("  protein:\n")
           outF.write("    ksize: [6,7,8,9,10,11]\n")
           outF.write("  dayhoff:\n")
           outF.write("    ksize: [14,15,16,17,18,19]\n")


# download genomes, count kmers, make species-level bloom filters
rule make_species_bloom_filters:
    input: os.path.join(out_dir, 'conf', '{basename}.{taxon}.bf.yml'),
    output: os.path.join(out_dir, 'bloom-filters', 'species', '{basename}.{taxon}.nodegraph'),
    resources:
        mem_mb=lambda wildcards, attempt: attempt *50000,
    threads: 10
    log: os.path.join(logs_dir, 'make-bf', '{basename}.{taxon}.make-bf.log')
    benchmark: os.path.join(logs_dir, 'make-bf', '{basename}.{taxon}.make-bf.benchmark')
    shell:
        """
        snakemake -s picklist-bf.snakefile --configfile {input} --cores {threads}
        """

rule write_lists_higher_rank_bloom_filters:
    # for higher ranks, need mapping from ranktaxon: all species taxa within this group
    # f"{out_dir}/script-nodegraphs/{basename}.protein-k{{ksize}}.nodegraph"
    input: lambda w: expand(os.path.join(out_dir, 'bloom-filters', 'species', '{basename}.{taxon}.{{alphabet}}-k{{ksize}}.nodegraph'), taxon = taxDF[taxDF[w.rank] == w.ranktaxon]["species"])
    output: os.path.join(out_dir, 'bloom-filters', '{rank}', '{basename}.{ranktaxon}.{{alphabet}}-k{{ksize}}.bf-list.txt')
    run:
        with open(str(output), 'w') as outF:
            for inF in input:
                outF.write(str(inf) + "\n")


rule make_higher_rank_bloom_filters:
    input: os.path.join(out_dir, 'bloom-filters', '{rank}', '{basename}.{ranktaxon}.{alphabet}-k{ksize}bf-list.txt') 
    output: os.path.join(out_dir, 'bloom-filters', '{rank}', '{basename}.{ranktaxon}.{{alphabet}}-k{{ksize}}.nodegraph')
    log: os.path.join(logs_dir, "update-bloom-filters", "{rank}", "{basename}.{ranktaxon}.{alphabet}-k{ksize}.log")
    benchmark: os.path.join(logs_dir, "update-bloom-filters", "{rank}", "{basename}.{ranktaxon}.{alphabet}-k{ksize}.benchmark")
    shell:
        """
        python combine-bloom-filters.py --from-file {input} --ksize {wildcards.ksize} {wildcards.alphabet} 2> {log}
        """


# compare using scaled = 200; unroll
rule jaccard_compare:
    input:
        picklist= os.path.join(out_dir, 'gtdb-rs202-taxonomic-picklists', '{rank}', '{basename}.{taxon}.csv'),
        zipfile=lambda w: expand(os.path.join(db_dir, 'gtdb-rs202.{fasta_type}.{{alphabet}}-k{{ksize}}.zip'), fasta_type = fasta_type[w.alphabet]),
    output:
        jaccard= os.path.join(out_dir, 'compare', '{rank}', '{basename}.{taxon}.{alphabet}-k{ksize}.jaccard.csv')
    log: os.path.join(logs_dir, 'compare', '{rank}', '{basename}.{taxon}.{alphabet}-k{ksize}.jaccard.log')
    benchmark: os.path.join(logs_dir, 'compare', '{rank}', '{basename}.{taxon}.{alphabet}-k{ksize}.jaccard.benchmark')
    shell:
        """
        sourmash sig compare --picklist {input.picklist}:ident:ident {input.zipfile} 
                             --csv {output.jaccard}  --{wildcards.alphabet} 2> {log}
        """


rule max_contain_compare:
    input:
        picklist= os.path.join(out_dir, 'gtdb-rs202-taxonomic-picklists', '{rank}', '{basename}.{taxon}.csv'),
        zipfile=lambda w: expand(os.path.join(db_dir, 'gtdb-rs202.{fasta_type}.{{alphabet}}-k{{ksize}}.zip'), fasta_type = fasta_type[w.alphabet]),
    output:
        max_contain= os.path.join(out_dir, 'compare', '{rank}', '{basename}.{taxon}.{alphabet}-k{ksize}.max_containment.csv')
    log: os.path.join(logs_dir, 'compare', '{rank}', '{basename}.{taxon}.{alphabet}-k{ksize}.max_contain.log')
    benchmark: os.path.join(logs_dir, 'compare', '{rank}', '{basename}.{taxon}.{alphabet}-k{ksize}.max_contain.benchmark')
    shell:
        """
        sourmash sig compare --picklist {input.picklist}:ident:ident {input.zipfile} --max-containment --csv {output.max_contain} 2>> {log}
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
        python  unroll-compare.py {input.max_contain} {input.jaccard} -o {output} 
        """

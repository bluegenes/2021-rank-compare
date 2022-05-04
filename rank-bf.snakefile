"""
Author: Tessa Pierce Ward, UC Davis Lab for Data Intensive Biology
Run: snakemake -j1 -n
"""

import csv
import pandas as pd
from collections import defaultdict
from sourmash.tax import tax_utils

taxinfo = config['taxonomy_csv']
taxfile_with_fastainfo = 'genbank/' + os.path.basename(taxinfo).rsplit('.csv')[0] + '.with-fastainfo.csv'
tax_info = pd.read_csv(config['taxonomy_csv'], header=0)

# add intended fastapaths
tax_info['genome_fastafile'] = 'genbank/genomes/'+ tax_info['ident'] + "_genomic.fna.gz"
tax_info['protein_fastafile'] = 'genbank/proteomes/'+ tax_info['ident'] + "_protein.faa.gz"

basename = config['basename']

ACCS =  tax_info["ident"]
tax_info.set_index("ident", inplace=True)

out_dir = config.get('output_dir', 'outputs')
logs_dir = os.path.join(out_dir, "logs")
basename = config.get('basename', 'gtdb-rs202')
database_dir = '/group/ctbrowngrp/gtdb/databases'
# how to incorporate db info??
search_databases = config.get('search_databases', ['gtdb-rs202'])
input_type = config.get('input_type', 'protein')

fasta_type= {'dna': "genomic", 'protein': 'protein', 'dayhoff': 'protein', 'hp': 'protein'}

gtdb_taxonomy=config.get('gtdb_taxonomy', '/group/ctbrowngrp/gtdb/gtdb-rs202.taxonomy.v2.with-strain.csv')
taxDF = pd.read_csv(gtdb_taxonomy)
# replace spaces with underscoare
taxDF = taxDF.replace(r" ", "_", regex=True) # regex allows searching for partial str
# temp hack - just do these:
#taxDF = taxDF[taxDF['phylum'] == 'p__Thermotogota']
#taxDF = taxDF[taxDF['species'] == 's__Natrinema_sp013456555']

taxonomyD = {}
rank_taxinfo = []
higher_rank_taxinfo = []

#input_picklist = config['picklist']
#highest_rank = config.get('highest_rank', 'phylum')

# make list of output bloom filter basenames
for rank in tax_utils.ascending_taxlist(include_strain=False):
    if rank not in ['superkingdom', 'phylum', 'class', 'family', 'order', 'genus']: # make species and genus first :)
        rank_taxonomies = taxDF[rank].unique() # unique taxa at each rank
        taxonomyD[rank] = rank_taxonomies
        rank_tax = expand(os.path.join(rank, f'{basename}.{{taxon}}'), taxon = rank_taxonomies) 
        rank_taxinfo.extend(rank_tax)
        #if rank != 'species':
        #    higher_rank_taxinfo.extend(rank_tax)

# check params are in the right format, build alpha-ksize combos
alphabet_info = config['alphabet_info']
alpha_ksizes=[]
for alpha, info in alphabet_info.items():
    # hack for just protein
    if alpha == 'protein':
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
        #expand(os.path.join(out_dir, 'gtdb-rs202-taxonomic-picklists', '{ranktaxinf}.csv'), ranktaxinf=rank_taxinfo),
        #expand(os.path.join(out_dir, 'bf-conf', '{ranktaxinf}.{alphak}.bf-list.txt'), ranktaxinf=higher_rank_taxinfo, alphak=alpha_ksizes),
        expand(os.path.join(out_dir, 'sourmash-nodegraph', '{ranktaxinf}.{alphak}.nodegraph'), ranktaxinf=rank_taxinfo, alphak=alpha_ksizes),


tablesizeD = {'species': "1e7", 'genus': "1e8", 'family': "2e8", "order": "4e8", "class": "5e8", "phylum": "1e9", "superkingdom": "1.5e10"}

#localrules: write_config_species_bloom_filters
#rule write_config_species_bloom_filters:
rule write_config_bloom_filters:
    input: 
        picklist= os.path.join(out_dir, 'gtdb-rs202-taxonomic-picklists', '{rank}', '{basename}.{taxon}.csv'),
    output: temp(os.path.join(out_dir, 'bf-conf', '{rank}', '{basename}.{taxon}.bf.yml')),
    params:
        tablesize=lambda w: tablesizeD[w.rank],
        nucl_k = alphabet_info['nucleotide']['ksize'],
        prot_k = alphabet_info['protein']['ksize'],
        #dayhoff_k = alphabet_info['dayhoff']['ksize']
        nucl_sc = alphabet_info['nucleotide']['scaled'],
        prot_sc = alphabet_info['protein']['scaled'],
        #dayhoff_sc = alphabet_info['dayhoff']['scaled']
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *3000,
        runtime=30,
    run:
       # location for all generated files
       inF = str(input)
       with open(str(output), 'w') as outF:
           outF.write(f"basename: {wildcards.rank}/gtdb-rs202.{wildcards.taxon}\n")
           outF.write(f"output_dir: {out_dir}\n")
           outF.write(f"taxonomy_csv: {inF}\n")
           outF.write(f"nodegraph_tablesize: {params.tablesize}\n")
           outF.write("run_alphabets:\n")
           #outF.write("  - nucleotide\n")
           outF.write("  - protein\n")
           #outF.write("  - dayhoff\n")
           outF.write("alphabet_info:\n")
           outF.write("  nucleotide:\n")
           outF.write(f"    ksize: {params.nucl_k}\n")
           outF.write(f"    scaled: {params.nucl_sc}\n")
           outF.write("  protein:\n")
           outF.write(f"    ksize: {params.prot_k}\n")
           outF.write(f"    scaled: {params.prot_sc}\n")
           #outF.write("  dayhoff:\n")
           #outF.write(f"    ksize: {params.dayhoff_k}\n")


# download genomes, count kmers, make species-level bloom filters
#rule make_species_bloom_filters:
rule make_bloom_filters:
    input: os.path.join(out_dir, 'bf-conf', '{rank}', '{basename}.{taxon}.bf.yml'),
    output: os.path.join(out_dir, 'sourmash-nodegraph', '{rank}', '{basename}.{taxon}.{alphabet}-k{ksize}.nodegraph'),
    resources:
        mem_mb=lambda wildcards, attempt: attempt *3000,
        runtime=3000,
    threads: 1
    log: os.path.join(logs_dir, 'sourmash-nodegraph', '{rank}', '{basename}.{taxon}.{alphabet}-k{ksize}.sourmash-nodegraph.log')
    benchmark: os.path.join(logs_dir, 'sourmash-nodegraph', '{rank}', '{basename}.{taxon}.{alphabet}-k{ksize}.sourmash-nodegraph.benchmark')
    shell:
        """
        snakemake -s picklist-bf.snakefile --configfile {input} --cores {threads} --until make_sourmash_nodegraph_protein --nolock 2> {log}
        """

#rule write_lists_higher_rank_bloom_filters:
#    # for higher ranks, need mapping from ranktaxon: all species taxa within this group
#    # f"{out_dir}/script-nodegraphs/{basename}.protein-k{{ksize}}.nodegraph"
#    input: 
#        lambda w: expand(os.path.join(out_dir, 'sourmash-nodegraph', 'species', '{{basename}}.{taxon}.{{alphabet}}-k{{ksize}}.nodegraph'), taxon = taxDF[taxDF[w.rank] == w.ranktaxon]["species"])
#    output: os.path.join(out_dir, 'bf-conf', '{rank}', '{basename}.{ranktaxon}.{alphabet}-k{ksize}.bf-list.txt')
#    run:
#        with open(str(output), 'w') as outF:
#            for inF in input:
#                outF.write(str(inF) + "\n")
#


#rule make_higher_rank_bloom_filters:
#    input: 
#        os.path.join(out_dir, 'bf-conf', '{rank}', '{basename}.{ranktaxon}.{alphabet}-k{ksize}.bf-list.txt') 
#    output: 
#        os.path.join(out_dir, 'sourmash-nodegraph', '{rank}', '{basename}.{ranktaxon}.{alphabet}-k{ksize}.nodegraph')
#    log: 
#        os.path.join(logs_dir, "combine-nodegraph", "{rank}", "{basename}.{ranktaxon}.{alphabet}-k{ksize}.log")
#    benchmark: 
#        os.path.join(logs_dir, "combine-nodegraph", "{rank}", "{basename}.{ranktaxon}.{alphabet}-k{ksize}.benchmark")
#    params:
#        tablesize = lambda w: tablesizeD[w.rank],
#    wildcard_constraints:
#        #rank='superkingdom|phylum|class|order|family|genus'
#        rank='phylum|class|order|family|genus'
#    shell:
#        """
#        python combine-bloom-filters.py --from-file {input} \
#               --ksize {wildcards.ksize} --tablesize {params.tablesize} \
#               --output {output} 2> {log}
#        """

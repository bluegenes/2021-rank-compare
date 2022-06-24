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

#taxDF.set_index("ident", inplace=True)
all_comparisons = []
ranktaxD = defaultdict(list)
for ranktax in alltax_at_rank[:50]:
    # get query_idents
    subDF = taxDF[(taxDF[compare_rank] == ranktax)]
    all_accs = subDF["ident"].to_list()
    num_acc = len(all_accs)
    # testing:let's start smaller: ignore >100 comparisons
    if 2 <= num_acc <= 500: # can't compare a single genome :)
        rt = ranktax.replace(' ', '_') # no spaces plsss
        ranktaxD[rt] = all_accs

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
        # pyani
        os.path.join(out_dir, "pyani", f"{compare_rank}.pyani-ANIm.csv.gz"),
        os.path.join(out_dir, "pyani", f"{compare_rank}.pyani-ANIb.csv.gz"),


def get_genomes_for_pyani(w):
    these_accs = ranktaxD[w.ranktax]
    genomes = []
    for acc in these_accs:
        genomes.append(fa_info.at[acc, 'genome_filename'])
    return genomes 


### split into genome folders + unzip fna files and generate classes/labels, then run pyANI index and compare 
# make a folder with just the genomes in a single ranktax
localrules: make_pyani_ranktax_folder
rule make_pyani_ranktax_folder:
    input:
        get_genomes_for_pyani
    output: 
        classes=os.path.join(out_dir, "pyani", "{ranktax}", "py.classes.txt"),
        labels=os.path.join(out_dir, "pyani", "{ranktax}", "py.labels.txt")
    params:
        acc_list = lambda w: ranktaxD[w.ranktax],
        taxdir = lambda w: os.path.join(out_dir, 'pyani', w.ranktax),
    run:
        import hashlib
        os.makedirs(params.taxdir, exist_ok=True)
        with open(output.classes, 'w') as out_classes:
            with open(output.labels, 'w') as out_labels:
                for acc in params.acc_list:
                    fn = fa_info.at[acc, 'genome_filename']
                    dest = os.path.join(params.taxdir, f"{acc}_genomic.fna.gz")
                    dest_unz = os.path.join(params.taxdir, f"{acc}_genomic.fna")
                    md5_file = os.path.join(params.taxdir, f"{acc}_genomic.md5")
                    shell("cp {fn} {dest}")
                    shell("gunzip -c {fn} > {dest_unz}")
                    # get md5 of unzipped fna
                    with open(dest_unz, "rb") as f:
                        bytes = f.read()
                        md5 = hashlib.md5(bytes).hexdigest()
                    # write to md5 file
                    with open(md5_file, 'w') as mfile:
                        mfile.write(f"{md5}\t{dest_unz}\n")
                    fna_base = os.path.basename(dest_unz).rsplit('.fna')[0]
                    out_classes.write(f"{md5}\t{fna_base}\t{acc}\n")
                    out_labels.write(f"{md5}\t{fna_base}\t{acc}\n")


rule pyani_index_and_createdb:
    input: 
        classes=os.path.join(out_dir, "pyani", "{ranktax}", "py.classes.txt"),
        labels=os.path.join(out_dir, "pyani", "{ranktax}", "py.labels.txt")
    output:
        classes=os.path.join(out_dir, "pyani", "{ranktax}/classes.txt"),
        labels=os.path.join(out_dir, "pyani", "{ranktax}/labels.txt"),
        db=os.path.join(out_dir, "pyani", "{ranktax}/.pyani-{ranktax}/pyanidb")
    params:
        taxdir = lambda w: os.path.join(out_dir, 'pyani', w.ranktax),
        pyanidb = lambda w: os.path.join(out_dir, 'pyani', w.ranktax, f".pyani-{w.ranktax}/pyanidb"),
        classes_basename = "classes.txt",
        labels_basename = "labels.txt"
    log: os.path.join(logs_dir, "pyani", "{ranktax}.index-and-createdb.log")
    benchmark: os.path.join(logs_dir, "pyani", "{ranktax}.index-and-createdb.benchmark")
    conda: "conf/envs/pyani0.3.yml"
    shell:
        """
        pyani index -i {params.taxdir} --classes {params.classes_basename} --labels {params.labels_basename} 
        pyani createdb --dbpath {params.pyanidb} -v -l {log}
        """

    
rule pyani_ANIm:
    input:  
        classes=os.path.join(out_dir, "pyani","{ranktax}/py.classes.txt"),
        labels=os.path.join(out_dir, "pyani", "{ranktax}/py.labels.txt"),
        idx_classes=os.path.join(out_dir, "pyani", "{ranktax}/classes.txt"),
        idx_labels=os.path.join(out_dir, "pyani", "{ranktax}/labels.txt")
    output: 
        directory(os.path.join(out_dir, "pyani", "{ranktax}/ANIm_results/nucmer_output")),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *5000,
        runtime=1200,
    params:
        pyanidb = lambda w: os.path.join(out_dir, 'pyani', w.ranktax, f".pyani-{w.ranktax}/pyanidb"),
        genome_dir = lambda w: os.path.join(out_dir, 'pyani', w.ranktax),
        output_dir = lambda w: os.path.join(out_dir, 'pyani', w.ranktax, "ANIm_results"),
    log: os.path.join(logs_dir, "pyani_anim", "{ranktax}.pyANI-aniM.log")
    benchmark: os.path.join(logs_dir, "pyani_anim", "{ranktax}.pyANI-aniM.benchmark")
    conda: "conf/envs/pyani0.3.yml"
    shell:
        """
        pyani anim --dbpath {params.pyanidb} --name {wildcards.ranktax} \
             --classes {input.classes} --labels {input.labels} \
            -i {params.genome_dir} -o {params.output_dir} \
            -l {log} -v
        """

rule pyANI_ANIb:
    input:  
        classes=os.path.join(out_dir, "pyani", "{ranktax}/py.classes.txt"),
        labels=os.path.join(out_dir, "pyani", "{ranktax}/py.labels.txt"),
        idx_classes=os.path.join(out_dir, "pyani", "{ranktax}/classes.txt"),
        idx_labels=os.path.join(out_dir, "pyani", "{ranktax}/labels.txt")
    output: 
        covF= os.path.join(out_dir, "pyani/{ranktax}/ANIb_results","ANIb_alignment_coverage.tab"),
        lenF= os.path.join(out_dir, "pyani/{ranktax}/ANIb_results","ANIb_alignment_lengths.tab"),
        hadF= os.path.join(out_dir, "pyani/{ranktax}/ANIb_results","ANIb_hadamard.tab"),
        idF=  os.path.join(out_dir, "pyani/{ranktax}/ANIb_results","ANIb_percentage_identity.tab"),
        seF=  os.path.join(out_dir, "pyani/{ranktax}/ANIb_results","ANIb_similarity_errors.tab"),
        bn =  os.path.join(out_dir, "pyani/{ranktax}/ANIb_results","blastn_output.tar.gz"),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *5000,
        runtime=1200,
    params:
        #pyanidb = lambda w: os.path.join(out_dir, 'pyani/paths', w.path, f".pyani-{w.path}/pyanidb"),
        genome_dir = lambda w: os.path.join(out_dir, 'pyani', w.ranktax),
        output_dir = lambda w: os.path.join(out_dir, 'pyani', w.ranktax, "ANIb_results"),
    log: os.path.join(logs_dir, "pyani_anib", "{ranktax}.pyANI-anib.log")
    benchmark: os.path.join(logs_dir, "pyani_anib", "{ranktax}.pyANI-anib.benchmark")
    conda: "conf/envs/pyani0.2.yml"
    shell:
        """
        average_nucleotide_identity.py -i {params.genome_dir} \
             -o {params.output_dir} -m ANIb -v \
             --labels {input.labels} --classes {input.classes} \
             --force > {log}
        """

rule pyANI_report_ANIm:
    input:  
        os.path.join(out_dir, "pyani", "{ranktax}/ANIm_results/nucmer_output"),
    output: 
        idF = os.path.join(out_dir, "pyani/{ranktax}/ANIm_results/matrix_identity_1.tab"),
        lenF = os.path.join(out_dir, "pyani/{ranktax}/ANIm_results/matrix_aln_lengths_1.tab"),
        covF = os.path.join(out_dir, "pyani/{ranktax}/ANIm_results/matrix_coverage_1.tab"),
        seF = os.path.join(out_dir, "pyani/{ranktax}/ANIm_results/matrix_sim_errors_1.tab"),
        hadF = os.path.join(out_dir, "pyani/{ranktax}/ANIm_results/matrix_hadamard_1.tab"),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *5000,
        runtime=1200,
    params:
        ANIm_dir = lambda w: os.path.join(out_dir, 'pyani', w.ranktax, "ANIm_results"),
        pyanidb = lambda w: os.path.join(out_dir, 'pyani', w.ranktax, f".pyani-{w.ranktax}/pyanidb"),
    log: os.path.join(logs_dir, "pyani", "{ranktax}/{ranktax}.pyANI-aniM.log")
    benchmark: os.path.join(logs_dir, "pyani", "{ranktax}/{ranktax}.pyANI-aniM.benchmark")
    conda: "conf/envs/pyani0.3.yml"
    shell:
        """
        pyani report -o {params.ANIm_dir} --dbpath {params.pyanidb} --runs --run_results 1 --formats stdout -l {log}
        pyani report -v -o {params.ANIm_dir} --formats stdout --run_matrices 1  \
        --dbpath {params.pyanidb} -l {log}
        """

localrules: aggregate_ranktax_anim
rule aggregate_ranktax_anim:
    input:
        idF = os.path.join(out_dir, "pyani/{ranktax}/ANIm_results/matrix_identity_1.tab"),
        lenF = os.path.join(out_dir, "pyani/{ranktax}/ANIm_results/matrix_aln_lengths_1.tab"),
        covF = os.path.join(out_dir, "pyani/{ranktax}/ANIm_results/matrix_coverage_1.tab"),
        seF = os.path.join(out_dir, "pyani/{ranktax}/ANIm_results/matrix_sim_errors_1.tab"),
        hadF = os.path.join(out_dir, "pyani/{ranktax}/ANIm_results/matrix_hadamard_1.tab"),
    output:
        os.path.join(out_dir, "pyani", "{ranktax}/ANIm_results/{ranktax}.pyani.csv"),
    params:
        results_dir =  os.path.join(out_dir, "pyani/{ranktax}/ANIm_results"),
        compare_rank = compare_rank,
    shell:
        """
        python aggregate-pyani-results.py {params.results_dir} --rank {params.compare_rank} --ranktax-name {wildcards.ranktax} --output-csv {output}
        """

localrules: aggregate_ranktax_anib
rule aggregate_ranktax_anib:
    input:
        covF= os.path.join(out_dir, "pyani/{ranktax}/ANIb_results", "ANIb_alignment_coverage.tab"),
        lenF= os.path.join(out_dir, "pyani/{ranktax}/ANIb_results","ANIb_alignment_lengths.tab"),
        hadF= os.path.join(out_dir, "pyani/{ranktax}/ANIb_results","ANIb_hadamard.tab"),
        idF=  os.path.join(out_dir, "pyani/{ranktax}/ANIb_results","ANIb_percentage_identity.tab"),
        seF=  os.path.join(out_dir, "pyani/{ranktax}/ANIb_results","ANIb_similarity_errors.tab"),
    output:
        os.path.join(out_dir, "pyani", "{ranktax}/ANIb_results/{ranktax}.pyani.csv"),
    params:
        results_dir =  os.path.join(out_dir, "pyani/{ranktax}/ANIb_results"),
        compare_rank = compare_rank,
    shell:
        """
        python aggregate-pyani-results.py {params.results_dir} --rank {params.compare_rank} --ranktax-name {wildcards.ranktax} --output-csv {output} --pyani-version v0.2
        """

localrules: aggregate_all_anim
rule aggregate_all_anim:
    input:
        expand(os.path.join(out_dir, "pyani", "{ranktax}/ANIm_results/{ranktax}.pyani.csv"), ranktax=ranktaxD.keys())
    output:
        os.path.join(out_dir, "pyani", f"{compare_rank}.pyani-ANIm.csv.gz")
    run:
        # aggreate all csv.gzs --> single csv
        aggDF = pd.concat([pd.read_csv(str(csv), sep=",") for csv in input])
        aggDF.to_csv(str(output), index=False)

localrules: aggregate_all_anib
rule aggregate_all_anib:
    input:
        expand(os.path.join(out_dir, "pyani", "{ranktax}/ANIb_results/{ranktax}.pyani.csv"), ranktax=ranktaxD.keys())
    output:
        os.path.join(out_dir, "pyani", f"{compare_rank}.pyani-ANIb.csv.gz")
    run:
        # aggreate all csv.gzs --> single csv
        aggDF = pd.concat([pd.read_csv(str(csv), sep=",") for csv in input])
        aggDF.to_csv(str(output), index=False)


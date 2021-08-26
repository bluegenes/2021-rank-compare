import csv
import pandas as pd

#configfile: "conf/gtdb-r202.yml"
configfile: "conf/R_gnavus.yml"
taxinfo = config['taxonomy_csv']
taxfile_with_fastainfo = 'genbank/' + os.path.basename(taxinfo).rsplit('.csv')[0] + '.with-fastainfo.csv'
tax_info = pd.read_csv(config['taxonomy_csv'], header=0)

# add intended fastapaths
tax_info['genome_fastafile'] = 'genbank/genomes/'+ tax_info['ident'] + "_genomic.fna.gz"
tax_info['protein_fastafile'] = 'genbank/proteomes/'+ tax_info['ident'] + "_protein.faa.gz"

basename = config['basename']

ACCS =  tax_info["ident"]
tax_info.set_index("ident", inplace=True)

logs_dir = 'genbank/logs'
prodigal_dir = 'genbank/prodigal'
out_dir = config['output_dir']

class Checkpoint_MakePattern:
    def __init__(self, pattern):
        self.pattern = pattern

    def update_fapaths(self, basename=basename, acc=None, ksize=None):
        with open(f"genbank/{basename}.prodigal-list.txt", "rt") as fp:
            accessions = [ os.path.basename(x).rstrip().rsplit("_protein.faa.gz")[0] for x in fp ]
            for accession in accessions:
                # update protein fastapath
                tax_info.at[accession, "protein_fastafile"] = f"{prodigal_dir}/{accession}_protein.faa.gz"
            # write fastapaths to file
            tax_info.to_csv(taxfile_with_fastainfo)
            
            return tax_info["protein_fastafile"]
    
    def find_protein_fasta(self, acc=None, ksize=None, alphabet=None):
        if acc:
            return tax_info.at[acc, "protein_fastafile"]

    def __call__(self, w):
        global checkpoints

        # wait for the results of 'check_proteins'; this will trigger an
        # exception until that rule has been run.
        checkpoints.check_proteins.get(**w)
        prot_fastafiles = self.update_fapaths()
        single_fasta = self.find_protein_fasta(**w)
        if single_fasta:
            prot_fastafiles = single_fasta
        pattern = expand(self.pattern, fastafile=prot_fastafiles, **w)
        return pattern


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
    #config["alphabet_info"][alpha]["select_params"] = expand("{alpha}-k{ksize}-scaled{scaled}", alpha=alpha, scaled=scaled, ksize=ksize)
    #these_alpha_ksize = expand("{alpha}-k{{ksize}}", ksize = ksize)
    if alpha in ['nucleotide', 'protein']:
        alpha_ksizes += expand(f"{alpha}-k{{ksize}}", ksize = ksize)

# some easier vars
alphabet_info = config['alphabet_info']
prot_ksizes = alphabet_info['protein']['ksize']
nucl_ksizes = alphabet_info['nucleotide']['ksize']

print(alpha_ksizes)

rule all:
    input:
        Checkpoint_MakePattern("{fastafile}"),
        ancient(expand("genbank/genomes/{acc}_genomic.fna.gz", acc=ACCS)),
        expand(f"{out_dir}/count-kmers/{{acc}}.nucleotide-k{{ksize}}.unique-kmers.txt", acc=ACCS, ksize = nucl_ksizes),
        expand(f"{out_dir}/count-kmers/{{acc}}.protein-k{{ksize}}.unique-kmers.txt", acc=ACCS, ksize = prot_ksizes),
        f"genbank/kmer-counts.csv",
        #expand(f"{out_dir}/{basename}.{{alphak}}.nodegraph", alphak=alpha_ksizes),
        expand(f"{out_dir}/script-nodegraphs/{basename}.{{alphak}}.nodegraph",alphak=alpha_ksizes),
        #expand("genbank/genomes/{acc}_genomic.fna.gz", acc=ACCS),
        #expand("genbank/proteomes/{acc}_protein.faa.gz", acc=ACCS)

# download genbank genome details; make an info.csv file for entry.
rule make_genome_info_csv:
    output:
        csvfile = 'genbank/info/{acc}.info.csv'
    shell: """
        python -Werror genbank_genomes.py {wildcards.acc} \
            --output {output.csvfile}
    """

# download actual genomes!
rule download_matching_genome_wc:
     input:
         csvfile = ancient('genbank/info/{acc}.info.csv')
     output:
         genome = "genbank/genomes/{acc}_genomic.fna.gz"
     run:
         with open(input.csvfile, 'rt') as infp:
             r = csv.DictReader(infp)
             rows = list(r)
             assert len(rows) == 1
             row = rows[0]
             acc = row['acc']
             assert wildcards.acc.startswith(acc)
             url = row['genome_url']
             name = row['ncbi_tax_name']

             print(f"downloading genome for acc {acc}/{name} from NCBI...",
                   file=sys.stderr)
             with open(output.genome, 'wb') as outfp:
                 with urllib.request.urlopen(url) as response:
                     content = response.read()
                     outfp.write(content)
                     print(f"...wrote {len(content)} bytes to {output.genome}",
                           file=sys.stderr)

# download proteome if possible
rule download_matching_proteome_wc:
    input:
        csvfile = ancient('genbank/info/{acc}.info.csv')
    output:
        proteome = "genbank/proteomes/{acc}_protein.faa.gz"
    run:
        with open(input.csvfile, 'rt') as infp:
            r = csv.DictReader(infp)
            rows = list(r)
            assert len(rows) == 1
            row = rows[0]
            acc = row['acc']
            assert wildcards.acc.startswith(acc)
            url = row['protein_url']
            name = row['ncbi_tax_name']

            print(f"downloading proteome for acc {acc}/{name} from NCBI...",
                file=sys.stderr)
            with open(output.proteome, 'wb') as outfp:
                try:
                    with urllib.request.urlopen(url) as response:
                        content = response.read()
                        outfp.write(content)
                        print(f"...wrote {len(content)} bytes to {output.proteome}",
                              file=sys.stderr)
                except:
                    shell('touch {output}')


rule check_protein_fastas:
    input:
        fasta=ancient(expand("genbank/proteomes/{acc}_protein.faa.gz", acc=ACCS))
    output: f"genbank/{basename}.prodigal-list.txt"
    run:
        with open(str(output), 'w') as out:
            for fasta in input.fasta:
                fasta_str = str(fasta)
                if os.path.getsize(fasta_str) == 0:
                    out.write(f"{fasta_str}\n")


# make the checkpoint work
localrules: check_proteins
checkpoint check_proteins:
    input:
        f"genbank/{basename}.prodigal-list.txt"
    output: temp(touch(f"genbank/.{basename}.check_proteins"))


rule unzip_genome_for_prodigal:
    input:
        fasta="genbank/genomes/{acc}_genomic.fna.gz",
        checkpt=f"genbank/.{basename}.check_proteins",
    output:
        unzipped=temp("genbank/genomes/{acc}_genomic.fna"),
    log: os.path.join(logs_dir, "gzip", "{acc}.genome_unzip.log")
    benchmark: os.path.join(logs_dir, "gzip", "{acc}.genome_unzip.benchmark")
    shell:
        """
        gunzip -c {input.fasta} > {output.unzipped} 2> {log}
        """

rule prodigal_genomes_for_failed_protein_downloads:
    input:
        fasta="genbank/genomes/{acc}_genomic.fna",
        checkpt=f"genbank/.{basename}.check_proteins",
    output:
        proteins=f"{prodigal_dir}/{{acc}}_protein.faa",
        genes=f"{prodigal_dir}/{{acc}}_genes.fna"
    conda: "conf/envs/prodigal-env.yml"
    log: os.path.join(logs_dir, "prodigal", "{acc}.prodigal.log")
    benchmark: os.path.join(logs_dir, "prodigal", "{acc}.prodigal.benchmark")
    shell:
        """
         prodigal -i {input.fasta} -a {output.proteins} -o {output.genes} 2> {log}
        """

rule gzip_prodigal_proteins:
    input:
        proteins=f"{prodigal_dir}/{{acc}}_protein.faa",
    output:
        f"{prodigal_dir}/{{acc}}_protein.faa.gz"
    log: os.path.join(logs_dir, 'gzip_prodigal', '{acc}.gzip.log')
    shell:
        """
        gzip -c {input} > {output} 2> {log}
        """

rule unique_kmers_protein:
    # need to use checkpoint here to make sure we get the updated proteome location
    input: fasta=Checkpoint_MakePattern("{fastafile}")
    output:
        f"{out_dir}/count-kmers/{{acc}}.protein-k{{ksize}}.unique-kmers.txt"
    log: f"{logs_dir}/unique_kmers/{{acc}}.protein-k{{ksize}}.log"
    conda: "conf/envs/khmer-env.yml"
    shell:
        """
        unique-kmers.py -k {wildcards.ksize} -e 0.005 {input.fasta} --report {output}
        """

rule unique_kmers_nucleotide:
    input: 
        fasta=lambda w: tax_info.at[w.acc, "genome_fastafile"]
    output:
        f"{out_dir}/count-kmers/{{acc}}.nucleotide-k{{ksize}}.unique-kmers.txt"
    log: f"{logs_dir}/unique_kmers/{{acc}}.nucleotide-k{{ksize}}.log"
    conda: "conf/envs/khmer-env.yml"
    shell:
        """
        unique-kmers.py -k {wildcards.ksize} -e 0.005 {input.fasta} --report {output}
        """

rule aggregate_unique_kmer_info:
    input:
        nucl=expand(f"{out_dir}/count-kmers/{{acc}}.nucleotide-k{{ksize}}.unique-kmers.txt", acc= ACCS, ksize=nucl_ksizes),
        prot=expand(f"{out_dir}/count-kmers/{{acc}}.protein-k{{ksize}}.unique-kmers.txt", acc=ACCS, ksize=prot_ksizes),
    output:
        f"genbank/kmer-counts.csv",
    run:
        from collections import defaultdict
        kcount = defaultdict(dict)
        for inf in input.nucl:
            inf_str = str(inf)
            with open(inf_str, 'r') as fp:
                alphabet='nucleotide'
                accession = os.path.basename(inf_str).rsplit('.nucleotide')[0]
                info = fp.readline()
                n_kmers,ksize= info.split(' ')[:2]
                kcount[accession][f'nucleotide-k{ksize}'] = n_kmers
        for inf in input.prot:
            inf_str = str(inf)
            with open(inf_str, 'r') as fp:
                alphabet='protein'
                accession = os.path.basename(inf_str).rsplit('.protein')[0]
                info = fp.readline()
                n_kmers,ksize= info.split(' ')[:2]
                kcount[accession][f'protein-k{ksize}'] = n_kmers
        
        # make pandas dataframe form dict
        df = pd.DataFrame.from_dict(kcount, orient="index")
        df.to_csv(str(output), index_label='accession')

rule make_bloom_filter_cmdline:
    #how do I make a bigger tablesize!??
    input: fasta=Checkpoint_MakePattern("{fastafile}")
    output: f"{out_dir}/{basename}.{{alphabet}}-k{{ksize}}.nodegraph"
    log: f"{logs_dir}/nodegraph/{basename}.{{alphabet}}-k{{ksize}}.log"
    threads: 10
    params:
        cmd_mem=1e9, #bc conversion...!?
    resources:
        mem=100000,
    conda: "conf/envs/khmer-env.yml"
    shell:
        """
        load-graph.py -k {wildcards.ksize} -T {threads} -M {params.cmd_mem} {output} {input} 2> {log}
        """

rule make_protein_bloom_filter_script:
    input: fasta=Checkpoint_MakePattern("{fastafile}")
    output: f"{out_dir}/script-nodegraphs/{basename}.protein-k{{ksize}}.nodegraph"
    log: f"{logs_dir}/script-nodegraphs/{basename}.protein-k{{ksize}}.log"
    threads: 10
    resources:
        mem=100000,
    #conda: "conf/envs/bf.yml"
    shell:
        """
        python protein-bf.py {input} --output {output} -k {wildcards.ksize} --alphabet protein 2> {log}
        """

rule make_nucl_bloom_filter_script:
    input: expand("genbank/genomes/{acc}_genomic.fna.gz", acc=ACCS),
    output: f"{out_dir}/script-nodegraphs/{basename}.nucleotide-k{{ksize}}.nodegraph"
    log: f"{logs_dir}/script-nodegraphs/{basename}.nucleotide-k{{ksize}}.log"
    threads: 10
    resources:
        mem=100000,
    #conda: "conf/envs/bf.yml"
    shell:
        """
        python protein-bf.py {input} --output {output} -k {wildcards.ksize} --alphabet nucleotide 2> {log}
        """

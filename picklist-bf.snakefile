import csv
import pandas as pd

taxinfo = config['taxonomy_csv']
taxfile_with_fastainfo = 'genbank/' + os.path.basename(taxinfo).rsplit('.csv')[0] + '.with-fastainfo.csv'
tax_info = pd.read_csv(config['taxonomy_csv'], header=0)

# add intended fastapaths
tax_info['genome_fastafile'] = 'genbank/genomes/'+ tax_info['ident'] + "_genomic.fna.gz"
tax_info['protein_fastafile'] = 'genbank/proteomes/'+ tax_info['ident'] + "_protein.faa.gz"

basename = config['basename']

ACCS =  tax_info["ident"]
tax_info.set_index("ident", inplace=True)

#logs_dir = 'genbank/logs'
out_dir = config['output_dir']
logs_dir = os.path.join(out_dir, 'logs')
#prodigal_dir = os.path.join(out_dir, 'prodigal')
#logs_dir = 'genbank/logs'
prodigal_dir = 'genbank/prodigal'

class Checkpoint_MakePattern:
    def __init__(self, pattern):
        self.pattern = pattern
        self.prot_fastafiles = None

    def update_fapaths(self, basename=basename, acc=None, ksize=None):
        with open(f"genbank/{basename}.prodigal-list.txt", "rt") as fp:
            accessions = [ os.path.basename(x).rstrip().rsplit("_protein.faa.gz")[0] for x in fp ]
            for accession in accessions:
                # update protein fastapath
                tax_info.at[accession, "protein_fastafile"] = f"{prodigal_dir}/{accession}_protein.faa.gz"
            # write fastapaths to file
            tax_info.to_csv(taxfile_with_fastainfo)

            return tax_info["protein_fastafile"], tax_info

    def find_protein_fasta(self, acc=None, ksize=None, alphabet=None):
        if acc:
            return tax_info.at[acc, "protein_fastafile"]

    def __call__(self, w):
        global checkpoints
        global updated_taxD

        # wait for the results of 'check_proteins'; this will trigger an
        # exception until that rule has been run.
        checkpoints.check_proteins.get(**w)
        if self.prot_fastafiles is None:
            self.prot_fastafiles, updated_taxD = self.update_fapaths()

        # single fasta instead
        single_fasta = self.find_protein_fasta(**w)
        if single_fasta:
            pattern = expand(self.pattern, fastafile=single_fasta, **w)
        else:
            pattern = expand(self.pattern, fastafile=self.prot_fastafiles, **w)
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
    if alpha in ['protein']: #['nucleotide']: #['protein']:#, 'dayhoff']: #['nucleotide', 'protein']:
        alpha_ksizes += expand(f"{alpha}-k{{ksize}}", ksize = ksize)

# some easier vars
alphabet_info = config['alphabet_info']
prot_ksizes = alphabet_info['protein']['ksize']
nucl_ksizes = alphabet_info['nucleotide']['ksize']

print(alpha_ksizes)

rule all:
    input:
     #   Checkpoint_MakePattern("{fastafile}"),
     #   ancient(expand("genbank/genomes/{acc}_genomic.fna.gz", acc=ACCS)),
     #   expand(f"{out_dir}/count-kmers/{{acc}}.nucleotide-k{{ksize}}.unique-kmers.txt", acc=ACCS, ksize = nucl_ksizes),
        #expand(f"{out_dir}/count-kmers/{{acc}}.protein-k{{ksize}}.unique-kmers.txt", acc=ACCS, ksize = prot_ksizes),
        #f"genbank/kmer-counts.csv",
        #expand(f"{out_dir}/sourmash-nodegraph/{basename}.{{alphak}}.nodegraph",alphak=alpha_ksizes),
        #expand("genbank/genomes/{acc}_genomic.fna.gz", acc=ACCS),
        #expand("genbank/proteomes/{acc}_protein.faa.gz", acc=ACCS)
        expand(f"{out_dir}/sourmash-nodegraph/{basename}.{{alphak}}.nodegraph",alphak=alpha_ksizes),
        #expand(f"{out_dir}/unique-kmers/{basename}.{{alphak}}.unique-kmers.txt",alphak=alpha_ksizes),
        #expand(f"{out_dir}/count-unique-kmers/{basename}.{{alphak}}.unique-kmers.txt", alphak=alpha_ksizes)

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
    output: touch(f"genbank/.{basename}.check_proteins")


rule unzip_genome_for_prodigal:
    input:
        fasta="genbank/genomes/{acc}_genomic.fna.gz",
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


rule write_protein_fastalist: 
    # need to use checkpoint here to make sure we get the updated proteome location
    input: fasta=ancient(Checkpoint_MakePattern("{fastafile}"))
    output:
        f"{out_dir}/fastalists/{basename}.protein.fastalist.txt"
    log: f"{logs_dir}/write-protein-fastalist/{basename}.protein.log"
    run:
        with open(str(output), "w") as out:
            for inF in input:
                fasta_file = str(inF)
                out.write(f"{fasta_file}\n")


#rule picklist_count_kmers:
#    # need to use checkpoint here to make sure we get the updated proteome location
#    input:  f"{out_dir}/fastalists/{basename}.protein.fastalist.txt" 
#    output: 
#        num_unique=f"{out_dir}/count-unique-kmers/{basename}.{{alphabet}}-k{{ksize}}.unique-kmers.txt",
#        fasta_nunique=f"{out_dir}/kmer-counts/{basename}.{{alphabet}}-k{{ksize}}.kmer-counts.txt",
#    params:
#        basename = basename,
#    log: f"{logs_dir}/count-kmers/{basename}.{{alphabet}}-k{{ksize}}.log"
#    shell:
#        #python picklist-shared-kmers.py {input} \
#        """
#        python hll-count-kmers.py {input} \
#          --basename {params.basename} \
#          --alphabet {wildcards.alphabet} \
#          --ksize {wildcards.ksize} \
#          --output-num-unique {output.num_unique} \
#          --output-fasta-nunique {output.fasta_nunique} 2> {log}
#        """


#rule write_nucleotide_fastalist: 
#    # need to use checkpoint here to make sure we get the updated proteome location
#    input: fasta=ancient(lambda w: tax_info.at[w.acc, "genome_fastafile"]) 
#    output:
#        f"{out_dir}/fastalists/{{acc}}.nucl.fastalist.txt"
#    log: f"{logs_dir}/write-nucl-fastalist/{{acc}}.protein.log"
#    run:
#        with open(str(output), "w") as out:
#            for inF in input:
#                fasta_file = str(inF)
#                out.write(f"{fasta_file}\n")


# this is for each individual genome/proteome
rule unique_kmers_protein:
    # need to use checkpoint here to make sure we get the updated proteome location
    input: fasta=ancient(Checkpoint_MakePattern("{fastafile}"))
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
        fasta=ancient(lambda w: tax_info.at[w.acc, "genome_fastafile"])
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


### this will probably not work with _all_the 250k genomes -- will run out of room on the command line
rule picklist_unique_kmers_protein_khmer:
    # need to use checkpoint here to make sure we get the updated proteome location
    input: fasta=ancient(Checkpoint_MakePattern("{fastafile}"))
    output: f"{out_dir}/unique-kmers/{basename}.protein-k{{ksize}}.unique-kmers.txt"
    log: f"{logs_dir}/unique_kmers/{basename}.protein-k{{ksize}}.log"
    conda: "conf/envs/khmer-env.yml"
    shell:
        """
        unique-kmers.py -k {wildcards.ksize} -e 0.005 {input.fasta} --report {output}
        """


#rule write_fasta_csv:
#    input: fasta=ancient(Checkpoint_MakePattern("{fastafile}"))
#    output: taxfile_with_fastainfo

# read new fastas from file instead of getting each time?
rule make_sourmash_nodegraph_protein:
#    input: fasta=ancient(Checkpoint_MakePattern("{fastafile}"))
    input: 
         fastalist=f"{out_dir}/fastalists/{basename}.protein.fastalist.txt",
    output: f"{out_dir}/sourmash-nodegraph/{basename}.protein-k{{ksize}}.nodegraph"
    log: f"{logs_dir}/sourmash-nodegraph/{basename}.protein-k{{ksize}}.log"
    benchmark: f"{logs_dir}/sourmash-nodegraph/{basename}.protein-k{{ksize}}.benchmark"
    params:
        tablesize=config.get('nodegraph_tablesize', 1e10),
        # full gtdb: 2e11 was plenty, could even go smaller
    threads: 1 
    resources:
        mem=600000,
    shell:
        """
        python sourmash-nodegraph.py --input-filelist {input.fastalist} \
               --output {output} -k {wildcards.ksize} \
               --alphabet protein --tablesize {params.tablesize} 2> {log}
        """

rule make_sourmash_nodegraph_dayhoff:
    input: fasta=ancient(Checkpoint_MakePattern("{fastafile}"))
    output: f"{out_dir}/sourmash-nodegraph/{basename}.dayhoff-k{{ksize}}.nodegraph"
    log: f"{logs_dir}/sourmash-nodegraph/{basename}.dayhoff-k{{ksize}}.log"
    benchmark: f"{logs_dir}/sourmash-nodegraph/{basename}.dayhoff-k{{ksize}}.benchmark"
    threads: 1
    resources:
        mem=150000,
    shell:
        """
        python sourmash-nodegraph.py {input} --output {output} -k {wildcards.ksize} --alphabet dayhoff --tablesize 1e10 2> {log}
        """

rule make_sourmash_nodegraph_nucl:
    input: ancient(expand("genbank/genomes/{acc}_genomic.fna.gz", acc=ACCS)),
    output: f"{out_dir}/sourmash-nodegraph/{basename}.nucleotide-k{{ksize}}.nodegraph"
    log: f"{logs_dir}/sourmash-nodegraph/{basename}.nucleotide-k{{ksize}}.log"
    benchmark: f"{logs_dir}/sourmash-nodegraph/{basename}.nucleotide-k{{ksize}}.benchmark"
    threads: 1
    resources:
        mem=150000,
    shell:
        """
        python sourmash-nodegraph.py {input} --output {output} -k {wildcards.ksize} --alphabet nucleotide --tablesize 1e12 2> {log}
        """

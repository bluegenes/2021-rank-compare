import csv
import pandas as pd

configfile: "conf/gtdb-rs207.yml"
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
        self.prot_fastafiles=None #self.update_fapaths()
        if os.path.exists("genbank/{basename}.prodigal-list.txt"):
            self.prot_fastafiles=self.update_fapaths()

    def update_fapaths(self, basename=basename, acc=None, ksize=None):
        if os.path.exists("genbank/{basename}.prodigal-list.txt"):
            with open(f"genbank/{basename}.prodigal-list.txt", "rt") as fp:
                accessions = [ os.path.basename(x).rstrip().rsplit("_protein.faa.gz")[0] for x in fp ]
                for accession in accessions:
                    # update protein fastapath
                    tax_info.at[accession, "protein_fastafile"] = f"{prodigal_dir}/{accession}_protein.faa.gz"
                # write fastapaths to file
                tax_info.to_csv(taxfile_with_fastainfo)
                # finally, return list of protein files
        return tax_info["protein_fastafile"]
    
    def find_protein_fasta(self, acc=None, ksize=None, alphabet=None):
        if acc:
            return tax_info.at[acc, "protein_fastafile"]

    def __call__(self, w):
        global checkpoints

        # wait for the results of 'check_proteins'; this will trigger an
        # exception until that rule has been run.
        checkpoints.check_proteins.get(**w)
        if self.prot_fastafiles is None:
            self.prot_fastafiles = self.update_fapaths()

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
    if alpha in ['nucleotide', 'protein']:
        alpha_ksizes += expand(f"{alpha}-k{{ksize}}", ksize = ksize)

# some easier vars
alphabet_info = config['alphabet_info']
prot_ksizes = alphabet_info['protein']['ksize']
nucl_ksizes = alphabet_info['nucleotide']['ksize']

print(alpha_ksizes)

rule all:
    input:
        #f"genbank/{basename}.prodigal-list.txt",
        os.path.join(out_dir, f"{basename}.fromfile.csv"),
        expand(f"{out_dir}/databases/{basename}.protein-k{{ksize}}.zip", ksize=prot_ksizes),
        expand(f"{out_dir}/databases/species/{basename}.protein-k{{ksize}}.species.zip", ksize=prot_ksizes),
        #Checkpoint_MakePattern("{fastafile}"),
        #ancient(expand("genbank/genomes/{acc}_genomic.fna.gz", acc=ACCS)),
        #ancient(expand("genbank/proteomes/{acc}_protein.faa.gz", acc=ACCS)),
        #expand("genbank/genomes/{acc}_genomic.fna.gz", acc=ACCS),
        #expand("genbank/proteomes/{acc}_protein.faa.gz", acc=ACCS)
        #expand(f"{basename}.{ak}.zip", ak=alpha_ksizes)
        #expand(f"{basename}.protein-k{k}.zip", k=prot_ksizes)

# download genbank genome details; make an info.csv file for entry.
rule make_genome_info_csv:
    output:
        csvfile = 'genbank/info/{acc}.info.csv'
    threads: 1
    resources:
        mem_mb=3000,
        disk_mb=5000,
        runtime=90,
        time=60,
        partition="bml", #low2
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
     threads: 1
     resources:
         mem_mb=3000,
         disk_mb=5000,
         runtime=60,
         time=90,
         partition="bml",#low2
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
    threads: 1
    resources:
         mem_mb=3000,
         runtime=60,
         time=90,
         partition="bml",#"low2", # bml
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


localrules: check_protein_fastas
rule check_protein_fastas:
    input:
        fasta=ancient(expand("genbank/proteomes/{acc}_protein.faa.gz", acc=ACCS))
    output: f"genbank/{basename}.prodigal-list.txt"
    threads: 1
    resources:
         mem_mb=3000,
         runtime=60,
         time=90,
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
    threads: 1
    output: temp(touch(f"genbank/.{basename}.check_proteins"))


rule unzip_genome_for_prodigal:
    input:
        fasta="genbank/genomes/{acc}_genomic.fna.gz",
        checkpt=f"genbank/.{basename}.check_proteins",
    output:
        unzipped=temp("genbank/genomes/{acc}_genomic.fna"),
    resources:
         mem_mb=3000,
         runtime=60,
         time=90,
         partition="bml",#"low2", # bml
    log: os.path.join(logs_dir, "gzip", "{acc}.genome_unzip.log")
    benchmark: os.path.join(logs_dir, "gzip", "{acc}.genome_unzip.benchmark")
    threads: 1
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
    threads: 1
    resources:
         mem_mb=3000,
         runtime=60,
         time=90,
         partition="bml",#"low2", # bml
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
    resources:
         mem_mb=3000,
         runtime=60,
         time=90,
         partition="bml",#"low2", # bml
    threads: 1
    log: os.path.join(logs_dir, 'gzip_prodigal', '{acc}.gzip.log')
    shell:
        """
        gzip -c {input} > {output} 2> {log}
        """

# need updated tax_info, so need to make sure checkpoint makepattern runs
rule build_fromfile_csv:
    input:  
        #proteomes=ancient(expand("genbank/proteomes/{acc}_protein.faa.gz", acc=ACCS)),
        plist = f"genbank/{basename}.prodigal-list.txt",
        proteomes = ancient(Checkpoint_MakePattern("{fastafile}")),
    resources:
         mem_mb=3000,
         runtime=60,
         time=90,
         partition="bml",#"low2", # bml
    output:
        fromfile_csv = os.path.join(out_dir, f"{basename}.fromfile.csv")
    run:
        fromfile_csvDF = tax_info[["signame", "genome_fastafile", "protein_fastafile"]]
        with open(str(input.plist), "rt") as fp:
            accessions = [ os.path.basename(x).rstrip().rsplit("_protein.faa.gz")[0] for x in fp ]
            for accession in accessions:
                # update protein fastapath
                fromfile_csvDF.at[accession, "protein_fastafile"] = f"{prodigal_dir}/{accession}_protein.faa.gz"

        # select columns from tax_info for sketch fromfile (cols needed: name,genome_filename,protein_filename)
        fromfile_csvDF.rename(columns={"signame":"name", "genome_fastafile": "genome_filename", "protein_fastafile": "protein_filename"}, inplace=True)
        fromfile_csvDF.to_csv(str(output.fromfile_csv))


rule sketch_protein_fromfile:
## just do proteins, bc titus already built genomic db
    input: 
        fromfile_csv = os.path.join(out_dir, f"{basename}.fromfile.csv")
    output: f"{out_dir}/databases/{basename}.protein-k{{ksize}}.zip"
    log: f"{logs_dir}/databases/{basename}.protein-k{{ksize}}.log"
    benchmark: f"{logs_dir}/sketch-fromfile/{basename}.protein-k{{ksize}}.benchmark"
    threads: 1
    resources:
        mem_mb=10000,
        time=20000,
    conda: "conf/envs/sourmash-4.4.yml" 
    shell:
        """
        sourmash sketch fromfile {input.fromfile_csv} -p protein,scaled=200,k={wildcards.ksize} -o {output} 2> {log}
        """
        #sourmash sketch fromfile output.rank-compare/gtdb-rs207.fromfile.csv -o gtdb-rs207.no-prodigal.zip -p k=10,scaled=200,protein

rule merge_to_species_db: 
    input: 
        database= f"{out_dir}/databases/{basename}.protein-k{{ksize}}.zip",
        taxonomy= config['taxonomy_csv'],
    output: f"{out_dir}/databases/species/{basename}.protein-k{{ksize}}.species.zip"
    log: f"{logs_dir}/databases/species/{basename}.protein-k{{ksize}}.merge-to-species.log"
    benchmark: f"{logs_dir}/merge-to-species/{basename}.protein-k{{ksize}}.merge-to-species.benchmark"
    threads: 1
    resources:
        mem_mb=10000,
        time=20000,
    conda: "conf/envs/sourmash-4.4.yml"
    shell:
        """
        python ../database-examples/mass-merge.py {input.database} --from-spreadsheet {input.taxonomy} \
                                                  --merge-col "species" -k {wildcards.ksize} --protein \
                                                  -o {output} 2> {log}
        """
        #sourmash sketch fromfile output.rank-compare/gtdb-rs207.fromfile.csv -o gtdb-rs207.no-prodigal.zip -p k=10,scaled=200,protein

01-simplesnakemake
================
2026-07-13

# Intro

In the previous section we showed how to perform and document a simple
bioinformatics pipeline. However, this previous section had one major
weakness: we manually performed all analysis steps for each of the
files. This is a waste of time and source of error. So it may be a good
idea to turn this into a Snakemake pipeline. In this section we show how
to generate SnakeMake pipelines. It’s more of a tutorial, starting very
simple and gradually adding complexity. Nevertheless it also
re-emphasizes some of the bioinformatics principles introduced before
(in section 01).

# the data

The analysis will be performed in a novel folder
‘/home/robert-kofler/analysis/2026-tutorial/02-snakemake’. We will use
the same files than used in the previous analysis
‘/home/robert-kofler/analysis/2026-tutorial/01-compFewStrains’ A major
bioinformatics principle introduced previously is to have all files
necessary for an analysis bundled in the folder where the analysis will
be performed. To avoid duplicating data we will use hard links:

``` bash
mkdir refg
mkdir refg-wg
mkdir rawreads
for i in ../../01-compFewStrains/rawdata/*.fastq; do ln $i . ; done # rawreads
for i in ../../01-compFewStrains/refg/*fasta*; do ln $i . ; done    # refg
for i in ../../01-compFewStrains/refg-wg/*fasta*; do ln $i . ; done # refg-wg 
```

Lets install snakemake

``` bash
# install snakemake; it needs an older Python version... (and i updated the sra-tools as well)
conda activate tutorial # see previous tutorial 01
conda install python=3.12 snakemake sra-tools
# this installed: snakemake                                9.23.1           hdfd78af_2            bioconda
```

# analysis

We follow the major principle: start as simple as possible and gradually
add complexity

## mapping single file with snakemake

First lets map just a single file, that is explicitely mentioned in the
Snakefile

``` python
# create the file
micro Snakefile
```

the content of the Snakefile

``` python
# Snakefiles are Python, so intendation matters
rule bwa_map:
    input:
        ref = "../refg/dmel-tes-scg.fasta",
        fq = "../rawdata/2004-I38_1.fastq" 
    output:
        "2004-I38.sort.bam"
    shell:
        "bwa mem {input.ref} {input.fq} | samtools sort -@ 4 -o {output} -"
# input means: these files must exist before I, the rule, can proceed
# output means: I will generate these files
```

lets run it

``` python
snakemake --cores 8
```

## mapping some specified files

lets map some files that we specify manually in the command line; so
this is already much more versatile than the previous version

``` python
mkdir map-some
micro Snakefile
```

Here is the snakemake file

``` python
rule bwa_map:
    input:
        ref = "../refg/dmel-tes-scg.fasta",
        fq = "../rawdata/{sample}_1.fastq" 
    output:
        "{sample}.sort.bam"
    shell:
        "bwa mem {input.ref} {input.fq} | samtools sort -@ 4 -o {output} -"
```

Now run with the following command. Note how Snakemake **pulls** the
requested files

``` python
snakemake --cores 8 2004-I38.sort.bam 2004-CO1.sort.bam
```

## mapping all files in folder

To map all files in a folder we need an additional rule; Note rule all
is just some random name, what matters is the order of the rules.
Snakemake always targets the very first rule

``` python
# find all samples in the folder rawdata; I only use the first read per this rule
SAMPLES = glob_wildcards("../rawdata/{sample}_1.fastq").sample

# glob_wildcards explanation:
# glob_wildcards uses regular expression pattern matching to check for files in a disc; it uses a capture group
# lets assume we have a folder with the files A_1.fastq  A_2.fastq  B_1.fastq  B_2.fastq  C_1.fastq 
# then glob_wildcards("rawdata/{sample}_1.fastq") returns
# Wildcards(sample=['C', 'B', 'A'])
# we could also specify several wildcards
# r = glob_wildcards("rawdata/{sample}_{read}.fastq")
# r.sample  # ['B', 'C', 'B', 'A', 'A']
# r.read    # ['2', '1', '1', '2', '1']


# request an output file for each sample
rule all:
    input:
        expand("{sample}.sort.bam", sample=SAMPLES)
# note that input means these files must exist before I can proceed;
# this rule all does not generate anything; hence not output
# this is a consumer rule stating: I want these files -> generate them
# so bascially rule all: should list all files that we want in the end

# rule bwa_map is a producer rule; it generates something
rule bwa_map:
    input:
        ref = "../refg/dmel-tes-scg.fasta",
        fq = "../rawdata/{sample}_1.fastq" 
    output:
        "{sample}.sort.bam"
    shell:
        "bwa mem {input.ref} {input.fq} | samtools sort -@ 4 -o {output} -"
```

Now this mapping would actually take a lot of time. So lets make a dry
run first and see if everything works (-n), also lets display all
commands (-p).

``` python
snakemake -n -p
```

if this works without problems lets un

``` python
snakemake --cores 8
```

## mapping all files in folder - building block

We will improve the snakemake pipeline for mapping all fastq files in a
folder by adding two important parts. First, we make the pipeline more
configurable with a yaml file. Second, we index the refernce
automatically when it is not yet indexed. Also we add some sanity
checks, like the assertion that files exist; and if they exist print
them

We will upload the pipeline to github in the folder bb-map. So it can
serve as a building block for quickly mapping several fastq files to a
reference genome.

When using it first configure in config.yaml than run snakemake

**Config.yaml**

``` yaml
ref: /home/robert-kofler/analysis/2026-tutorial/02-snakemake/refg/dmel-tes-scg.fasta
readdir: /home/robert-kofler/analysis/2026-tutorial/02-snakemake/rawdata
outdir: /home/robert-kofler/analysis/2026-tutorial/02-snakemake/bb-map
```

**Snakefile**

``` python
configfile: "config.yaml"

# define some frequently used variables
# this is merely for convenience
REF     = config["ref"]
READDIR = config["readdir"]
OUTDIR  = config["outdir"]
BWAIDX = multiext(REF, ".amb", ".ann", ".bwt", ".pac", ".sa")

# find all samples in the folder 
SAMPLES = glob_wildcards(READDIR+ "/{sample}.fastq").sample

# check if samples exists
# and print the names of samples
assert SAMPLES, f"no *.fastq found in {READDIR}"
print("samples:")
[print(i) for i in SAMPLES]
# also print other parameters
print()
print("reference: "+ REF)
print("input directory: "+ READDIR)
print("output directory: "+ OUTDIR)

# request an output file for each sample in the folder
rule all:
    input:
        expand(OUTDIR+"/{sample}.sort.bam", sample=SAMPLES)


rule bwa_index:
    input:
        REF
    output:
        BWAIDX
    shell:
        "bwa index {input}"


rule bwa_map:
    input:
        ref = REF,
        fq = READDIR+"/{sample}.fastq",
        idx = BWAIDX
    output:
        OUTDIR+"/{sample}.sort.bam"
    shell:
        "bwa mem {input.ref} {input.fq} | samtools sort -@ 4 -o {output} -"
```

## mapping + visualization

Now we can finally turn the pipeline from the first section (01) into a
full snakemake pipeline. We need to add the seqvista steps

**yaml**

``` yaml
ref: /home/robert-kofler/analysis/2026-tutorial/02-snakemake/refg/dmel-tes-scg.fasta
readdir: /home/robert-kofler/analysis/2026-tutorial/02-snakemake/rawdata
outdir: /home/robert-kofler/analysis/2026-tutorial/02-snakemake/map-seqvista
```

**Snakefile**

``` python
configfile: "config.yaml"

# define some frequently used variables
# this is merely for convenience
REF     = config["ref"]
READDIR = config["readdir"]
OUTDIR  = config["outdir"]
BWAIDX = multiext(REF, ".amb", ".ann", ".bwt", ".pac", ".sa")

# find all samples in the folder 
SAMPLES = glob_wildcards(READDIR+ "/{sample}_1.fastq").sample

# check if samples exist and print them
assert SAMPLES, f"no *.fastq found in {READDIR}"
print("samples:")
[print(i) for i in SAMPLES]
print()
print("reference: "+ REF)
print("input directory: "+ READDIR)
print("output directory: "+ OUTDIR)

# request an output file for each sample in the folder
rule all:
    input:
        expand(OUTDIR+"/{sample}.sort.bam.bai", sample=SAMPLES),
        expand(OUTDIR+"/{sample}.norm.so", sample=SAMPLES)



rule normalize_so:
    input:
        OUTDIR+"/{sample}.so"
    output:
        OUTDIR+"/{sample}.norm.so"
    shell:
        "normalize-so.py --so {input} --scg-begin 'Dmel' --output-file {output}"    


rule bam2so:
    input:
        ref = REF,
        bam = OUTDIR + "/{sample}.sort.bam"
    output:
        bam = OUTDIR + "/{sample}.so"
    shell:
        "bam2so.py --infile {input.bam} --fasta {input.ref} --output-file {output}"             
        # # bam2so.py --infile map-se/1958-Hikone.sort.bam --fasta refg/dmel-tes-scg.fasta --output-file sose2/tmp/1958-Hikone.so 


rule bwa_index:
    input:
        REF
    output:
        BWAIDX
    shell:
        "bwa index {input}"


rule samtools_index:
    input:
        OUTDIR + "/{sample}.sort.bam"
    output:
        OUTDIR + "/{sample}.sort.bam.bai"
    shell:
        "samtools index {input}"


rule bwa_map:
    input:
        ref = REF,
        fq = READDIR+"/{sample}_1.fastq",
        idx = BWAIDX
    output:
        OUTDIR+"/{sample}.sort.bam"
    shell:
        "bwa mem {input.ref} {input.fq} | samtools sort -@ 4 -o {output} -"
```

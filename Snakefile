###  Snakefile for allelspefic transcription
###
### input parent1, parent2 and F1

# Configuration
###############
configfile: "config.yaml" # where to find parameters
WORKING_DIR = config["working_dir"]
RESULT_DIR = config["result_dir"]

# fetch URL to transcriptome multi fasta from configfile
genome_P1_url = config["refs"]["genome_P1"]
transcriptome_P1_url_gtf= config["refs"]["transcriptome_P1_gtf"]
transcriptome_P1_fasta_url= config["refs"]["transcriptome_P1_fasta"]
genome_P2_url = config["refs"]["genome_P2"]
transcriptome_P2_gtf_url= config["refs"]["transcriptome_P2_gtf"]
transcriptome_P2_gtf_url= config["refs"]["transcriptome_P2_fasta"]

########################
# Samples and conditions
########################

# read the tabulated separated table containing the sample, condition and fastq file information∂DE
units = pd.read_table(config["units"], dtype=str).set_index(["sample"], drop=False)

# create lists containing the sample names and conditions
SAMPLES = units.index.get_level_values('sample').unique().tolist()
samples = pd.read_csv(config["units"], dtype=str,index_col=0,sep="\t")
#CONDITIONS = list(pd.read_table(config["units"])["condition"])
samplefile = config["units"]


###########################
# Input functions for rules
###########################

def sample_is_single_end(sample):
    """This function detect missing value in the column 2 of the units.tsv"""
    if "fq2" not in samples.columns:
        return True
    else:
        return pd.isnull(samples.loc[(sample), "fq2"])

def get_fastq(wildcards):
    """ This function checks if the sample has paired end or single end reads
    and returns 1 or 2 names of the fastq files """
    if sample_is_single_end(wildcards.sample):
        return samples.loc[(wildcards.sample), ["fq1"]].dropna()
    else:
        return samples.loc[(wildcards.sample), ["fq1", "fq2"]].dropna()

def get_trimmed(wildcards):
    """ This function checks if sample is paired end or single end
    and returns 1 or 2 names of the trimmed fastq files """
    if sample_is_single_end(wildcards.sample):
        return WORKING_DIR + "trimmed/" + wildcards.sample + "_R1_trimmed.fq.gz"
    else:
        return [WORKING_DIR + "trimmed/" + wildcards.sample + "_R1_trimmed.fq.gz", WORKING_DIR + "trimmed/" + wildcards.sample + "_R2_trimmed.fq.gz"]


#################
# Desired outputs
#################
rule all:
    input:


    message:
        "Job done! Removing temporary directory"

#######
# Rules
#######


#####################
# Download references
#####################

rule get_genome_P1_fasta:
    output:
        WORKING_DIR + "genome/genomeP1.fasta"
    message:
        "downloading the required genomic fasta file"
    conda:
        "envs/wget.yaml"
    shell:
        "wget -O {output} {genome_P1_url}"

rule get_transcriptome_P2_gtf:
    output:
        WORKING_DIR + "genome/transcriptomeP1.gff"
    message:
        "downloading required transcriptome gtf file"
    conda:
        "envs/wget.yaml"
    shell:
        "wget -O {output} {transcriptome_P1_url}"

rule get_genome_P1_fasta:
    output:
        WORKING_DIR + "genome/genomeP2.fasta"
    message:
        "downloading the required genomic fasta file"
    conda:
        "envs/wget.yaml"
    shell:
        "wget -O {output} {genome_P2_url}"

rule get_transcriptome_P2_gtf:
    output:
        WORKING_DIR + "genome/transcriptome_P2.gff"
    message:
        "downloading required transcriptome gtf file"
    conda:
        "envs/wget.yaml"
    shell:
        "wget -O {output} {transcriptome_P2_url}"

rule merge_genomes:
    input:
        P1  = WORKING_DIR + "genome/genomeP1.fasta"
        P2  = WORKING_DIR + "genome/genomeP2.fasta"
    output:
        WORKING_DIR + "genome/total.fasta"
    shell:
        "cat {input.P1} {input.P2} > {output}"


##################################################################################
#  Blast transcriptomes to get homologous gene couples between both parental lines
##################################################################################

# create transcriptome index, for blasting
rule get_transcriptome_P1_index:
    input:
        WORKING_DIR + "genome/transcriptome_P1.fasta"
    output:
        [WORKING_DIR + "genome/transcriptome_P1.fasta." + i for i in ("psq", "phr", "pin")]
    conda:
        "envs/blast.yaml"
    shell:
        "makeblastdb -in {input} -dbtype prot"

rule get_transcriptome_P2_index:
    input:
        WORKING_DIR + "genome/transcriptome_P2.fasta"
    output:
        [WORKING_DIR + "genome/transcriptome_P2.fasta." + i for i in ("psq", "phr", "pin")]
    conda:
        "envs/blast.yaml"
    shell:
        "makeblastdb -in {input} -dbtype prot"

# Do the blasts
rule blast_P1_to_P2:
    input:
        newTct     = WORKING_DIR + "genome/transcriptome_P1.fasta",
        refTct     = WORKING_DIR + "genome/transcriptome_P2.fasta",
        indexFiles = [WORKING_DIR + "genome/transcriptome_P2.fasta." + i for i in ("psq", "phr", "pin")]
    output:
        WORKING_DIR + "results/P1_to_P2_blast.txt"
    params:
        evalue     = str(config['blast']['evalue']),     # 1e-10
        outFmt     = str(config['blast']['outFmt']),     # 6 qseqid qlen slen evalue salltitles
        maxTargets = str(config['blast']['maxTargets'])  # 1bin/bash: indent: command not found
    threads:
        5
    conda:
        "envs/blast.yaml"
    shell:
        "blastx "
        "-query {input.newTct} "
        "-db {input.refTct} "
        "-outfmt \"{params.outFmt}\" "
        "-evalue {params.evalue} "
        "-out {output} "
        "-num_threads {threads} "
        "-max_target_seqs {params.maxTargets}"

rule blast_P2_to_P1:
    input:
        newTct     = WORKING_DIR + "genome/transcriptome_P2.fasta",
        refTct     = WORKING_DIR + "genome/transcriptome_P1.fasta",
        indexFiles = [WORKING_DIR + "genome/transcriptome_P1.fasta." + i for i in ("psq", "phr", "pin")]
    output:
        WORKING_DIR + "results/P2_to_P1_blast.txt"
    params:
        evalue     = str(config['blast']['evalue']),     # 1e-10
        outFmt     = str(config['blast']['outFmt']),     # 6 qseqid qlen slen evalue salltitles
        maxTargets = str(config['blast']['maxTargets'])  # 1bin/bash: indent: command not found
    threads:
        5
    conda:
        "envs/blast.yaml"
    shell:
        "blastx "
        "-query {input.newTct} "
        "-db {input.refTct} "
        "-outfmt \"{params.outFmt}\" "
        "-evalue {params.evalue} "
        "-out {output} "
        "-num_threads {threads} "
        "-max_target_seqs {params.maxTargets}"



##################################
# Fastp
##################################

rule fastp:
    input:
        get_fastq
    output:
        fq1  = WORKING_DIR + "trimmed/" + "{sample}_R1_trimmed.fq.gz",
        fq2  = WORKING_DIR + "trimmed/" + "{sample}_R2_trimmed.fq.gz",
        html = RESULT_DIR + "fastp/{sample}.html"
    message:"trimming {wildcards.sample} reads"
    threads: 10
    log:
        RESULT_DIR + "fastp/{sample}.log.txt"
    params:
        sampleName = "{sample}",
        qualified_quality_phred = config["fastp"]["qualified_quality_phred"]
    run:
        if sample_is_single_end(params.sampleName):
            shell("fastp --thread {threads}  --html {output.html} \
            --qualified_quality_phred {params.qualified_quality_phred} \
            --in1 {input} --out1 {output} \
            2> {log}; \
            touch {output.fq2}")
        else:
            shell("fastp --thread {threads}  --html {output.html} \
            --qualified_quality_phred {params.qualified_quality_phred} \
            --detect_adapter_for_pe \
            --in1 {input[0]} --in2 {input[1]} --out1 {output.fq1} --out2 {output.fq2}; \
            2> {log}")


###################################################################################################
# RNA-Seq read alignement to parental genomes, accepting multiple missmatches and with softclipping
###################################################################################################



# If you pass one file, STAR will consider these as single-end reads: --readFilesIn single_reads.fastq.

# If you pass two files, STAR will consider these as paired reads: --readFilesIn pair_1.fastq pair_2.fastq.



rule star_index_p1:
    input:
        WORKING_DIR + "genome/genome_P1.fasta"
    output:
        [WORKING_DIR + "genome/genome_P1." + str(i) + ".ht2" for i in range(1,9)]
    message:
        "indexing genome"
    params:
        WORKING_DIR + "genome/genome"
    threads: 10
    shell:
        "hisat2-build -p {threads} {input} {params} --quiet"

rule star_mapping_P1:
    input:
        get_trimmed,
        indexFiles = [WORKING_DIR + "genome/genome." + str(i) + ".ht2" for i in range(1,9)]
    output:
        bams  = WORKING_DIR + "mapped/{sample}.bam",
        sum   = RESULT_DIR + "logs/{sample}_sum.txt",
        met   = RESULT_DIR + "logs/{sample}_met.txt"
    params:
        indexName = WORKING_DIR + "genome/genome",
        sampleName = "{sample}"
    # conda:
    #     "envs/hisat_mapping.yaml"
    message:
        "mapping reads to genome to bam files."
    threads: 10
    run:
        if sample_is_single_end(params.sampleName):
            shell("hisat2 -p {threads} --summary-file {output.sum} --met-file {output.met} -x {params.indexName} \
            -U {input} | samtools view -Sb -F 4 -o {output.bams}")
        else:
            shell("hisat2 -p {threads} --summary-file {output.sum} --met-file {output.met} -x {params.indexName} \
            -1 {input[0]} -2 {input[1]} | samtools view -Sb -F 4 -o {output.bams}")




rule star_index_p2:
    input:
        WORKING_DIR + "genome/genome_P2.fasta"
    output:
        [WORKING_DIR + "genome/genome_P2." + str(i) + ".ht2" for i in range(1,9)]
    message:
        "indexing genome"
    params:
        WORKING_DIR + "genome/genome"
    threads: 10
    shell:
        "hisat2-build -p {threads} {input} {params} --quiet"

rule star_mapping_P2:
    input:
        get_trimmed,
        indexFiles = [WORKING_DIR + "genome/genome." + str(i) + ".ht2" for i in range(1,9)]
    output:
        bams  = WORKING_DIR + "mapped/{sample}.bam",
        sum   = RESULT_DIR + "logs/{sample}_sum.txt",
        met   = RESULT_DIR + "logs/{sample}_met.txt"
    params:
        indexName = WORKING_DIR + "genome/genome",
        sampleName = "{sample}"
    # conda:
    #     "envs/hisat_mapping.yaml"
    message:
        "mapping reads to genome to bam files."
    threads: 10
    run:
        if sample_is_single_end(params.sampleName):
            shell("hisat2 -p {threads} --summary-file {output.sum} --met-file {output.met} -x {params.indexName} \
            -U {input} | samtools view -Sb -F 4 -o {output.bams}")
        else:
            shell("hisat2 -p {threads} --summary-file {output.sum} --met-file {output.met} -x {params.indexName} \
            -1 {input[0]} -2 {input[1]} | samtools view -Sb -F 4 -o {output.bams}")



rule star_index_merged:
    input:
        WORKING_DIR + "genome/total.fasta"
    output:
        [WORKING_DIR + "genome/total." + str(i) + ".ht2" for i in range(1,9)]
    message:
        "indexing genome"
    params:
        WORKING_DIR + "genome/genome"
    threads: 10
    shell:
        "hisat2-build -p {threads} {input} {params} --quiet"






rule star_mapping_merged:
    input:
        get_trimmed,
        indexFiles = [WORKING_DIR + "genome/genome." + str(i) + ".ht2" for i in range(1,9)]
    output:
        bams  = WORKING_DIR + "mapped/{sample}.bam",
        sum   = RESULT_DIR + "logs/{sample}_sum.txt",
        met   = RESULT_DIR + "logs/{sample}_met.txt"
    params:
        indexName = WORKING_DIR + "genome/genome",
        sampleName = "{sample}"
    # conda:
    #     "envs/hisat_mapping.yaml"
    message:
        "mapping reads to genome to bam files."
    threads: 10
    run:
        if sample_is_single_end(params.sampleName):
            shell("hisat2 -p {threads} --summary-file {output.sum} --met-file {output.met} -x {params.indexName} \
            -U {input} | samtools view -Sb -F 4 -o {output.bams}")
        else:
            shell("hisat2 -p {threads} --summary-file {output.sum} --met-file {output.met} -x {params.indexName} \
            -1 {input[0]} -2 {input[1]} | samtools view -Sb -F 4 -o {output.bams}")




rule create_counts_table_P1:
    input:
        bams = expand(WORKING_DIR + "mapped/{sample}.bam", sample = SAMPLES),
        gff  = WORKING_DIR + "genome/transcriptome_P1.gff"
    output:
        RESULT_DIR + "counts_P1.txt"
    conda:
        "envs/subread.yaml"
    shell:
        "featureCounts -O -t mRNA -g ID -F 'gtf' -a {input.gff} -o {output} {input.bams}"


rule create_counts_table_P2:
    input:
        bams = expand(WORKING_DIR + "mapped/{sample}.bam", sample = SAMPLES),
        gff  = WORKING_DIR + "genome/transcriptome_P2.gff"
    output:
        RESULT_DIR + "counts_P2.txt"
    conda:
        "envs/subread.yaml"
    shell:
        "featureCounts -O -t mRNA -g ID -F 'gtf' -a {input.gff} -o {output} {input.bams}"


rule create_counts_table_merged_P1:
    input:
        bams = expand(WORKING_DIR + "mapped/{sample}.bam", sample = SAMPLES),
        gff  = WORKING_DIR + "genome/transcriptome_P1.gff"
    output:
        RESULT_DIR + "counts_merged_P1.txt"
    conda:
        "envs/subread.yaml"
    shell:
        "featureCounts -O -t mRNA -g ID -F 'gtf' -a {input.gff} -o {output} {input.bams}"



rule create_counts_table_merged_P2:
    input:
        bams = expand(WORKING_DIR + "mapped/{sample}.bam", sample = SAMPLES),
        gff  = WORKING_DIR + "genome/transcriptome_P2.gff"
    output:
        RESULT_DIR + "counts_merged_P2.txt"
    conda:
        "envs/subread.yaml"
    shell:
        "featureCounts -O -t mRNA -g ID -F 'gtf' -a {input.gff} -o {output} {input.bams}"

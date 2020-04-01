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
proteome_P1_fasta_url= config["refs"]["proteins_P1"]
genome_P2_url = config["refs"]["genome_P2"]
transcriptome_P2_gtf_url= config["refs"]["transcriptome_P2_gtf"]
proteome_P2_gtf_url= config["refs"]["proteins_P2"]

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
        WORKING_DIR + "genome/merged.fasta"'
        #[WORKING_DIR + "genome/proteomeP1.fasta." + i for i in ("psq", "phr", "pin")]
        #[WORKING_DIR + "genome/proteomeP2.fasta." + i for i in ("psq", "phr", "pin")]
        #WORKING_DIR + "results/P2_to_P1_blast.txt"
        #WORKING_DIR + "results/P1_to_P2_blast.txt"
        #RESULT_DIR + "counts_P1.txt",
        #RESULT_DIR + "counts_P2.txt",
        #RESULT_DIR + "counts_merged_P1.txt",
        #RESULT_DIR + "counts_merged_P2.txt"
    message:
        "Job done! Removing temporary directory"


#######
# Rules
#######


#####################
# Download references
#####################

rule get_genome_fastas:
    output:
        P1 = WORKING_DIR + "genome/genomeP1.fasta",
        P2 = WORKING_DIR + "genome/genomeP2.fasta"
    message:
        "downloading the required genomic fasta file"
    conda:
        "envs/wget.yaml"
    shell:
        "wget -O {output.P1} {genome_P1_url}; wget -O {output.P2} {genome_P2_url}"


rule get_transcriptome_gtfs:
    output:
        P1 = WORKING_DIR + "genome/transcriptomeP1.gff",
        P2 = WORKING_DIR + "genome/transcriptomeP2.gff"
    message:
        "downloading required transcriptome gtf file"
    conda:
        "envs/wget.yaml"
    shell:
        "wget -O {output.P1} {transcriptome_P1_url}; wget -O {output.P2} {transcriptome_P2_url}"


rule get_protein_fastas:
    output:
        P1 = WORKING_DIR + "genome/proteomeP1.fasta",
        P2 = WORKING_DIR + "genome/proteomeP2.fasta"
    message:
        "downloading required proteome fasta file"
    conda:
        "envs/wget.yaml"
    shell:
        "wget -O {output.P1} {proteome_P1_gtf_url}; wget -O {output.P2} {proteome_P2_gtf_url}"


rule merge_genomes:
    input:
        P1  = WORKING_DIR + "genome/genomeP1.fasta"
        P2  = WORKING_DIR + "genome/genomeP2.fasta"
    output:
        WORKING_DIR + "genome/merged.fasta"
    shell:
        "cat {input.P1} {input.P2} > {output}"


##################################################################################
#  Blast transcriptomes to get homologous gene couples between both parental lines
##################################################################################

""" 
To get the homologous genes from both perantal lines, blast will be done in to ways (P1 -> P2 and P2 -> P1)
If the in a blast two genes hit each other in both directions it is presumed a "homologous gene couple"
"""

# create transcriptome index, for blasting
rule get_transcriptome_P1_index:
    input:
        WORKING_DIR + "genome/proteomeP1.fasta"
    output:
        [WORKING_DIR + "genome/proteomeP1.fasta." + i for i in ("psq", "phr", "pin")]
    conda:
        "envs/blast.yaml"
    shell:
        "makeblastdb -in {input} -dbtype prot"


rule get_transcriptome_P2_index:
    input:
        WORKING_DIR + "genome/proteomeP2.fasta"
    output:
        [WORKING_DIR + "genome/proteomeP2.fasta." + i for i in ("psq", "phr", "pin")]
    conda:
        "envs/blast.yaml"
    shell:
        "makeblastdb -in {input} -dbtype prot"


# Do the blasts
rule blast_P1_to_P2:
    input:
        P1 = WORKING_DIR + "genome/proteomeP1.fasta",
        P2 = WORKING_DIR + "genome/proteomeP2.fasta",
        indexFiles = [WORKING_DIR + "genome/proteomeP2.fasta." + i for i in ("psq", "phr", "pin")]
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
        "-query {input.P1} "
        "-db {input.P2} "
        "-outfmt \"{params.outFmt}\" "
        "-evalue {params.evalue} "
        "-out {output} "
        "-num_threads {threads} "
        "-max_target_seqs {params.maxTargets}"


rule blast_P2_to_P1:
    input:
        P1 = WORKING_DIR + "genome/proteomeP1.fasta",
        P2 = WORKING_DIR + "genome/proteomeP2.fasta",
        indexFiles = [WORKING_DIR + "genome/proteomeP1.fasta." + i for i in ("psq", "phr", "pin")]
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
        "-query {input.P2} "
        "-db {input.P1} "
        "-outfmt \"{params.outFmt}\" "
        "-evalue {params.evalue} "
        "-out {output} "
        "-num_threads {threads} "
        "-max_target_seqs {params.maxTargets}"


#######################################
# Fastp (quality controle and trimmimg)
#######################################

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

rule star_index_P1:
    input:
        WORKING_DIR + "genome/genomeP1.fasta"
    output:
        WORKING_DIR + "genome/genomeP1/SAindex"
    message:
        "indexing genome"
    params:
        dir    = WORKING_DIR + "genome/genomeP1",
        prefix = "genomeP1"
    conda:
        "envs/star.yaml"
    threads: 10
    shell:
        "STAR --runMode genomeGenerate "
        "--genomeDir {params.dir} "
        "--outFileNamePrefix {params.prefix} " 
        "--genomeFastaFiles {input} "
        "--runThreadN {threads} "
        "--limitGenomeGenerateRAM 100000000000"


rule star_mapping_P1:
    input:
        fq         = get_trimmed,
        WORKING_DIR + "genome/genomeP1/SAindex"
    output:
        WORKING_DIR + "mapped/{sample}.bam"
    params:
        indexName = WORKING_DIR + "genome/genomeP1",
        sampleName = "{sample}"
    conda:
        "envs/star.yaml"
    message:
        "mapping reads to genome to bam files."
    threads: 10
    shell:
    	"STAR --genomeDir {params.indexName} "
        "--runThreadN {threads} "
        "--readFilesIn {fq} "
        "--readFilesCommand zcat "
        "--outFileNamePrefix {params.sampleName} "
        "--outSAMtype BAM SortedByCoordinate "
        "--outSAMunmapped None "
        "--outSAMattributes Standard"


rule star_index_p2:
    input:
        WORKING_DIR + "genome/genomeP2.fasta"
    output:
        WORKING_DIR + "genome/genomeP2/SAindex"
    message:
        "indexing genome"
    params:
        dir    = WORKING_DIR + "genome/genomeP2",
        prefix = "genomeP2"
    conda:
        "envs/star.yaml"
    threads: 10
    shell:
        "STAR --runMode genomeGenerate "
        "--genomeDir {params.dir} "
        "--outFileNamePrefix {params.prefix} " 
        "--genomeFastaFiles {input} "
        "--runThreadN {threads} "
        "--limitGenomeGenerateRAM 100000000000"


rule star_mapping_P2:
    input:
        fq    = get_trimmed,
        index = WORKING_DIR + "genome/genomeP2/SAindex"
    output:
        WORKING_DIR + "mapped/{sample}.bam"
    params:
        indexName = WORKING_DIR + "genome/genome",
        sampleName = "{sample}"
    conda:
        "envs/star.yaml"
    message:
        "mapping reads to genome to bam files."
    threads: 10
    shell:
    	"STAR --genomeDir {params.indexName} "
        "--runThreadN {threads} "
        "--readFilesIn {fq} "
        "--readFilesCommand zcat "
        "--outFileNamePrefix {params.sampleName} "
        "--outSAMtype BAM SortedByCoordinate "
        "--outSAMunmapped None "
        "--outSAMattributes Standard"


#######################################################################################
# RNA-Seq read alignement to merged genome, accepting no mismatches and no softclipping
#######################################################################################

rule star_index_merged:
    input:
        WORKING_DIR + "genome/merged.fasta"
    output:
         WORKING_DIR + "genome/merged/SAindex"
    message:
        "indexing genome"
    params:
        dir    = WORKING_DIR + "genome/merged",
        prefix = "merged"
    conda:
        "envs/star.yaml"
    threads: 10
    shell:
        "STAR --runMode genomeGenerate --genomeDir {params.dir} "
        "--outFileNamePrefix {params.prefix} " 
        "--genomeFastaFiles {input} "
        "--runThreadN {threads} "
        "--limitGenomeGenerateRAM 200000000000"


rule star_mapping_merged:
    input:
        get_trimmed,
        indexFiles = WORKING_DIR + "genome/merged/SAindex"
    output:
        bams  = WORKING_DIR + "mapped/{sample}.bam"
    params:
        indexName = WORKING_DIR + "genome/merged",
        sampleName = "{sample}"
    conda:
        "envs/star.yaml"
    message:
        "mapping reads to genome to bam files."
    threads: 10
    shell:
    	"STAR --genomeDir {params.indexName} "
        "--runThreadN {threads} "
        "--readFilesIn {fq} "
        "--readFilesCommand zcat "
        "--outFileNamePrefix {params.sampleName} "
        "--outSAMtype BAM SortedByCoordinate "
        "--outSAMunmapped None "
        "--outSAMattributes Standard "
        "--outFilterMismatchNmax 0 "
        "--alignEndsType EndToEnd "
        "--outFilterMultimapNmax 1"


############################
## create raw counts tables.
############################

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

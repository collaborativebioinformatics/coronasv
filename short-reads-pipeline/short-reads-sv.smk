# ---------------------------------------------------------------
# Snakemake file for detection of SVs in PE short-read sequencing
# Baylor Hackathon | Oct 2020 | CoronaSV Team
# ---------------------------------------------------------------

# CONFIG

REFERENCE = config['reference']

PAIRED = True
RUNS, = glob_wildcards("reads/{run}_1.fastq.gz")

# TARGETS

rule all:
  input:
     expand("results/variants/{run}.bwa.rg.dedup.delly.vcf", run = RUNS)

# REFERENCE PREP

rule reference_faidx:
  input:
    REFERENCE
  output:
    REFERENCE+".fai"
  shell:
    """
    samtools faidx {input}
    """

rule bwa_index:
  input:
    REFERENCE
  output:
    REFERENCE+".bwt"
  shell:
    """
    bwa index {input}
    """

# QC TRIMMING

rule trimming:
  input:
    forward = "reads/{run}_1.fastq.gz",
    reverse = "reads/{run}_2.fastq.gz"
  output:
    forward = "results/trimmed/{run}_1_val_1.fq.gz",
    reverse = "results/trimmed/{run}_2_val_2.fq.gz"
  params:
    "results/trimmed"
  threads: 5
  shell:
    """
    trim_galore -j {threads} -q 20 --phred33 --length 50 --illumina --paired -o {params} {input.forward} {input.reverse}
    """

# MAPPING

rule bwa:
  input:
    reverse = "results/trimmed/{run}_1_val_1.fq.gz",
    forward = "results/trimmed/{run}_2_val_2.fq.gz",
    reference = REFERENCE+".bwt"
  output:
    sam = temp("results/alignments/{run}.bwa.sam")
  params:
    reference = REFERENCE
  threads: 5
  shell:
    """
    bwa mem -t {threads} {params.reference} {input.reverse} {input.forward} > {output.sam}
    """

rule sam2bam:
  input:
    sam = "results/alignments/{run}.bwa.sam"
  output:
    bam = temp("results/alignments/{run}.bwa.bam")
  shell:
    """
    samtools view -Sb {input.sam} | samtools sort > {output.bam}
    """

# MAPPING POSTPROCESSING

rule read_groups:
  input:
    bam = "results/alignments/{run}.bwa.primary.bam"
  output:
    bam = temp("results/alignments/{run}.bwa.primary.rg.bam")
  shell:
    """
    picard AddOrReplaceReadGroups I={input.bam} O={output.bam} RGID=1 RGLB=lib RGPL=ILLUMINA RGPU=unit1 RGSM=Sample1
    """

rule deduplicate:
  input:
    bam = "results/alignments/{run}.bwa.rg.bam"
  output:
    bam = "results/alignments/{run}.bwa.rg.dedup.bam"
  params:
    txt = "results/alignments/{run}.bwa.rg.dedup.txt"
  shell:
    """
    picard MarkDuplicates REMOVE_DUPLICATES=true ASSUME_SORTED=true I={input.bam} O={output.bam} M={params.txt}
    """

rule bam_index:
  input:
    bam = "results/alignments/{run}.bwa.rg.dedup.bam"
  output:
    bam_index = "results/alignments/{run}.bwa.rg.dedup.bam.bai"
  shell:
    """
    samtools index {input.bam}
    """

# SV CALLING

rule delly:
  input:
    bam = "results/alignments/{run}.bwa.rg.dedup.bam",
    bai = "results/alignments/{run}.bwa.rg.dedup.bam.bai",
    ref = REFERENCE
  output:
    vcf = "results/variants/{run}.bwa.rg.dedup.delly.vcf"
  params:
    bcf = "results/variants/{run}.bwa.rg.dedup.delly.bcf"
  shell:
    """
    delly call -g {REFERENCE} {input.bam} -o {params.bcf}
    bcftools view {params.bcf} > {output.vcf}
    """

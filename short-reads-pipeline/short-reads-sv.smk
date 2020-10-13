# ---------------------------------------------------------------
# Snakemake file for detection of SVs in PE short-read sequencing
# Baylor Hackathon | Oct 2020 | CoronaSV Team
# ---------------------------------------------------------------

# CONFIG

REFERENCE = config['reference']
TYPE = config['type']

if TYPE == "paired":
  RUNS, = glob_wildcards("reads/{run}_1.fastq")
else:
  RUNS, = glob_wildcards("reads/{run}.fastq")

# TARGETS

rule all:
  input:
     expand("results/variants_delly/{run}.bwa.rg.dedup.delly.vcf", run=RUNS),
     expand("results/variants_manta/{run}/results/variants/tumorSV.vcf", run=RUNS),
     expand("results/variants_lumpy/{run}.bwa.rg.dedup.lumpy.vcf", run=RUNS)

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
    forward = "reads/{run}_1.fastq",
    reverse = "reads/{run}_2.fastq"
  output:
    forward = "results/trimmed/{run}_1_val_1.fq",
    reverse = "results/trimmed/{run}_2_val_2.fq"
  params:
    "results/trimmed"
  shell:
    """
    trim_galore -q 20 --phred33 --length 50 --illumina --paired -o {params} {input.forward} {input.reverse}
    """

# MAPPING

rule bwa:
  input:
    reverse = "results/trimmed/{run}_1_val_1.fq",
    forward = "results/trimmed/{run}_2_val_2.fq",
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
    bam = "results/alignments/{run}.bwa.bam"
  output:
    bam = temp("results/alignments/{run}.bwa.rg.bam")
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

# SV calling with Delly

rule delly:
  input:
    bam = "results/alignments/{run}.bwa.rg.dedup.bam",
    bai = "results/alignments/{run}.bwa.rg.dedup.bam.bai",
    ref = REFERENCE
  output:
    vcf = "results/variants_delly/{run}.bwa.rg.dedup.delly.vcf"
  params:
    bcf = "results/variants_delly/{run}.bwa.rg.dedup.delly.bcf"
  shell:
    """
    delly call -g {REFERENCE} {input.bam} -o {params.bcf}
    bcftools view {params.bcf} > {output.vcf}
    """

# SV calling with Manta

rule manta:
  input:
    bam = "results/alignments/{run}.bwa.rg.dedup.bam",
    bai = "results/alignments/{run}.bwa.rg.dedup.bam.bai",
    ref = REFERENCE
  output:
    "results/variants_manta/{run}/results/variants/tumorSV.vcf"
  params:
    outdir = "results/variants_manta/{run}",
    output = "results/variants_manta/{run}/results/variants/tumorSV.vcf.gz"
  shell:
    """
    configManta.py --tumorBam {input.bam} --referenceFasta {REFERENCE} --runDir {params.outdir}
    {params.outdir}/runWorkflow.py
    gunzip {params.output}
    """

# SV calling with Lumpy

rule lumpy:
  input:
    bam = "results/alignments/{run}.bwa.rg.dedup.bam"
  output:
    "results/variants_lumpy/{run}.bwa.rg.dedup.lumpy.vcf"
  shell:
    """
    lumpyexpress \
    -B {input.bam} \
    -o {output}
    """ 


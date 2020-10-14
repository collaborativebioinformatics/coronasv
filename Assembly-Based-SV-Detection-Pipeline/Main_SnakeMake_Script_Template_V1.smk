# This is a snakemake file/script for SV Detection within the SARS-COV2 genome using an assembly based approach

### Input Data Types:

### Maximillian Marin (mgmarin@g.harvard.edu)


#### Import Statements ####
import pandas as pd


# Define PATH to the reference genome to be used: 
refGenome_FA_PATH = config["RefGenome_FA_PATH"]
refGenome_GFF_PATH = config["RefGenome_GFF_PATH"]


# Define PATH of OUTPUT directory
output_Dir = config["output_dir"]


# Read in meta data regarding input 
df = pd.read_csv( config["inputSampleData_TSV"], sep='\t')



CoronaSV_Metadata_All_DF = pd.read_csv(  config["inputSampleData_TSV"]  , sep = "\t")

CoronaSV_Metadata_ONT_DNA_DF = CoronaSV_Metadata_All_DF[  (CoronaSV_Metadata_All_DF["Platform"] == "OXFORD_NANOPORE") & (CoronaSV_Metadata_All_DF["Assay_Type"] != "RNA-Seq") ]           

CoronaSV_Metadata_Illumina_PE_DF = CoronaSV_Metadata_All_DF[ (CoronaSV_Metadata_All_DF["LibraryLayout"] == "PAIRED") & (CoronaSV_Metadata_All_DF["Platform"] == "ILLUMINA") & (CoronaSV_Metadata_All_DF["Assay_Type"] != "RNA-Seq") ]           



# Save a list of SRA/ENA "Run" Accessions
input_SampleIDs_WiIllumina = list( CoronaSV_Metadata_Illumina_PE_DF["Run"].values )
input_SampleIDs_WiNanopore = list( CoronaSV_Metadata_ONT_DNA_DF["Run"].values )

#print(len(input_SampleIDs_WiIllumina))
#print(len(input_SampleIDs_WiNanopore))


rule all:
    input:
        expand(output_Dir + "/Illumina_PE_FQs/{sampleID_WiIllumina}_1.fastq", sampleID_WiIllumina=input_SampleIDs_WiIllumina),
        expand(output_Dir + "/Illumina_PE_FQs/{sampleID_WiIllumina}_2.fastq", sampleID_WiIllumina=input_SampleIDs_WiIllumina),
        expand(output_Dir + "/{sampleID_WiIllumina}/IlluminaWGS/Unicycler_SPAdes_Assembly/{sampleID_WiIllumina}.SPAdes.Assembly.fasta.fai", sampleID_WiIllumina=input_SampleIDs_WiIllumina),
        expand(output_Dir + "/{sampleID_WiIllumina}/VariantCalling/Minimap2_Alignment/{sampleID_WiIllumina}.minimap2.paftools.vcf", sampleID_WiIllumina=input_SampleIDs_WiIllumina),
        expand(output_Dir + "/{sampleID_WiIllumina}/VariantCalling/NucDiff_Analysis_{sampleID_WiIllumina}_V2_WiVCFout/results/{sampleID_WiIllumina}.NucDiff_ref_snps.vcf", sampleID_WiIllumina=input_SampleIDs_WiIllumina),
        expand(output_Dir + "/{sampleID_WiIllumina}/VariantCalling/SVanalyzer_SVrefine_Ill_SPAdes_Assembly_SV_Calling/{sampleID_WiIllumina}.MUMmer.SVrefine.vcf", sampleID_WiIllumina=input_SampleIDs_WiIllumina),


        expand(output_Dir + "/Nanopore_FQs/{sampleID_WiNanopore}.fastq", sampleID_WiNanopore=input_SampleIDs_WiNanopore),
        expand(output_Dir + "/Nanopore_MM2_sniffles/{sampleID_WiNanopore}.sniffles.vcf", sampleID_WiNanopore=input_SampleIDs_WiNanopore),
        expand(output_Dir + "/Nanopore_MM2_svim/{sampleID_WiNanopore}/variants.vcf", sampleID_WiNanopore=input_SampleIDs_WiNanopore),
        expand(output_Dir + "/Nanopore_MM2_cuteSV/{sampleID_WiNanopore}.cuteSV.vcf", sampleID_WiNanopore=input_SampleIDs_WiNanopore),

        expand(output_Dir +"results/variants_delly/{sampleID_WiIllumina}.bwa.rg.dedup.delly.vcf", sampleID_WiIllumina=input_SampleIDs_WiIllumina),
        expand(output_Dir +"results/variants_manta/{sampleID_WiIllumina}/results/variants/tumorSV.vcf", sampleID_WiIllumina=input_SampleIDs_WiIllumina),
        expand(output_Dir +"results/variants_lumpy/{sampleID_WiIllumina}.bwa.rg.dedup.lumpy.vcf", sampleID_WiIllumina=input_SampleIDs_WiIllumina)




rule download_Illumina_FQ_FromSRA_RunID:
    output: 
        FQ_1 = output_Dir + "/Illumina_PE_FQs/{sampleID_WiIllumina}_1.fastq",
        FQ_2 = output_Dir + "/Illumina_PE_FQs/{sampleID_WiIllumina}_2.fastq",

    params:
        target_Download_Dir = output_Dir + "/Illumina_PE_FQs/"
    #conda: "Envs/sratools_2_10_7_Conda.yml"
    
    shell:
        "fastq-dump --split-files {wildcards.sampleID_WiIllumina} --outdir {params.target_Download_Dir} \n"




rule download_Nanopore_FQ_FromSRA_RunID:
    output: 
        ONT_FQ = output_Dir + "/Nanopore_FQs/{sampleID_WiNanopore}.fastq",

    params:
        target_Download_Dir = output_Dir + "/Nanopore_FQs/"
    #conda: "Envs/sratools_2_10_7_Conda.yml"
    
    shell:
        "fastq-dump {wildcards.sampleID_WiNanopore} --outdir {params.target_Download_Dir} \n"




##################################################################
############ SV Calling w/ long read (ONT) alignment #############
##################################################################


########## A) Long Reads Quality Control ##########




rule nanoplot_QC:
    input:
        ONT_FQ = output_Dir + "/Nanopore_FQs/{sampleID_WiNanopore}.fastq",
    output:
         output_Dir + "/Nanopore/Nanopore_QC/{sampleID_WiNanopore}.NanoPlot/NanoPlot-report.html"
    threads: 2
    shell:
        "NanoPlot -t {threads} --fastq {input.ONT_reads_fq} -o {output_Dir}/{wildcards.sampleID_WiNanopore}/Nanopore/Nanopore_QC/{wildcards.sampleID_WiNanopore}.NanoPlot/"



rule nanofilt_QC:
    input:
        ONT_FQ = output_Dir + "/Nanopore_FQs/{sampleID_WiNanopore}.fastq",
    output:
        ONT_Filtered_FQ = output_Dir + "/Nanopore_FQs_Filtered/{sampleID_WiNanopore}.Filtered.fastq",

    threads: 2
    shell:
        "NanoFilt -q 10 --headcrop 50 {input.ONT_FQ} > {output.ONT_Filtered_FQ}"






########## B) Mapping long reads to reference genome ##########

rule ONT_Minimap2_Alignment_And_Processing_ONT_Reads:
    input:
        ONT_Filtered_FQ = output_Dir + "/Nanopore_FQs_Filtered/{sampleID_WiNanopore}.Filtered.fastq",
        Ref_FA = refGenome_FA_PATH,
    output:
        ONT_MM2_SAM = output_Dir + "/Nanopore_MM2_Alignments/{sampleID_WiNanopore}.minimap2.sam",

    threads: 1

    shell:
        "minimap2 -ax map-ont {input.Ref_FA} {input.ONT_Filtered_FQ} > {output.ONT_MM2_SAM}  \n"



rule samtools_Processing_ONT_Reads:
    input:
        ONT_MM2_SAM = output_Dir + "/Nanopore_MM2_Alignments/{sampleID_WiNanopore}.minimap2.sam",
        Ref_FA = refGenome_FA_PATH,
    output:
        ONT_MM2_Sorted_BAM = output_Dir + "/Nanopore_MM2_Alignments/{sampleID_WiNanopore}.minimap2.sorted.bam",
        ONT_MM2_Final_BAM = output_Dir + "/Nanopore_MM2_Alignments/{sampleID_WiNanopore}.minimap2.final.bam",
        ONT_MM2_Final_BAM_BAI = output_Dir + "/Nanopore_MM2_Alignments/{sampleID_WiNanopore}.minimap2.final.bam.bai",
    shell:
        "samtools view -bS {input.ONT_MM2_SAM} | samtools sort - > {output.ONT_MM2_Sorted_BAM}  \n"
        "samtools calmd {output.ONT_MM2_Sorted_BAM} {input.Ref_FA} | samtools view -bS > {output.ONT_MM2_Final_BAM}  \n"
        "samtools index {output.ONT_MM2_Final_BAM}  \n"




########## C) SV Calling using long reads aligned to reference genome ##########

####### Sniffles SV Calling ####### 

rule ONT_MM2_Sniffles:
    input:
        ONT_MM2_Final_BAM = output_Dir + "/Nanopore_MM2_Alignments/{sampleID_WiNanopore}.minimap2.final.bam",
    output:
        ONT_Sniffles_VCF = output_Dir + "/Nanopore_MM2_sniffles/{sampleID_WiNanopore}.sniffles.vcf",

    threads: 1

    shell:
        "sniffles -m {input.ONT_MM2_Final_BAM} -v {output.ONT_Sniffles_VCF}" 


####### Sniffles SV Calling ####### 

rule ONT_MM2_SVIM:
    input:
        ONT_MM2_Final_BAM = output_Dir + "/Nanopore_MM2_Alignments/{sampleID_WiNanopore}.minimap2.final.bam",
        Ref_FA = refGenome_FA_PATH,
    output:
        ONT_SVIM_VCF = output_Dir + "/Nanopore_MM2_svim/{sampleID_WiNanopore}/variants.vcf",
    shell:
        "svim alignment {output_Dir}/Nanopore_MM2_svim/{wildcards.sampleID_WiNanopore} {input.ONT_MM2_Final_BAM} {input.Ref_FA}" 



####### cuteSV SV Calling ####### 

rule ONT_MM2_cuteSV:
    input:
        ONT_MM2_Final_BAM = output_Dir + "/Nanopore_MM2_Alignments/{sampleID_WiNanopore}.minimap2.final.bam",
        Ref_FA = refGenome_FA_PATH,
    output:
        ONT_cuteSV_VCF = output_Dir + "/Nanopore_MM2_cuteSV/{sampleID_WiNanopore}.cuteSV.vcf",

    shell:
        "mkdir -p tmp/{wildcards.sampleID_WiNanopore}  \n"
        "cuteSV {input.ONT_MM2_Final_BAM} {input.Ref_FA} {output.ONT_cuteSV_VCF} tmp/{wildcards.sampleID_WiNanopore}" 





# adapter list from: https://github.com/stephenturner/adapters/blob/master/adapters_combined_256_unique.fasta
# can move adapter list to config file later
rule trimmomatic_Illumina_PE_Trimming:
    input:
        r1 = output_Dir + "/Illumina_PE_FQs/{sampleID_WiIllumina}_1.fastq",
        r2 = output_Dir + "/Illumina_PE_FQs/{sampleID_WiIllumina}_2.fastq",
    output:
        r1 = output_Dir + "/{sampleID_WiIllumina}/IlluminaWGS/FASTQs/{sampleID_WiIllumina}_1_trimmed.fastq",
        r2 = output_Dir + "/{sampleID_WiIllumina}/IlluminaWGS/FASTQs/{sampleID_WiIllumina}_2_trimmed.fastq",
        # reads where trimming entirely removed the mate
        r1_unpaired = output_Dir + "/{sampleID_WiIllumina}/IlluminaWGS/FASTQs/{sampleID_WiIllumina}_1_trimmed.unpaired.fastq",
        r2_unpaired = output_Dir + "/{sampleID_WiIllumina}/IlluminaWGS/FASTQs/{sampleID_WiIllumina}_2_trimmed.unpaired.fastq",
    log:
        "logs/trimmomatic/{sampleID_WiIllumina}.log"
    params:
        # list of trimmers (see manual)
        trimmer=["ILLUMINACLIP:'References/CustomTrimmoatic_IlluminaWGS_AdapterList.fasta':2:30:10:2:true SLIDINGWINDOW:4:20 MINLEN:75"],
        # optional parameters
        # extra=" "
    threads: 1
    wrapper:
        "0.38.0/bio/trimmomatic/pe"



################################################################################
######### UniCycler (Using SPAdes Internally) for short read assembly  #########
################################################################################

# Adding SPADes assembly with UNIcycler to analysis

### A) Assembly with SPAdes through Unicycler
rule unicycler_SPAdes_Assemble_IlluminaWGS:
    input:
        fq1_trimmed = output_Dir + "/{sampleID_WiIllumina}/IlluminaWGS/FASTQs/{sampleID_WiIllumina}_1_trimmed.fastq",
        fq2_trimmed = output_Dir + "/{sampleID_WiIllumina}/IlluminaWGS/FASTQs/{sampleID_WiIllumina}_2_trimmed.fastq",
    output:
        assembly_GFA = output_Dir + "/{sampleID_WiIllumina}/IlluminaWGS/Unicycler_SPAdes_Assembly/assembly.gfa",
        Ill_SPAdes_assembly_fa = output_Dir + "/{sampleID_WiIllumina}/IlluminaWGS/Unicycler_SPAdes_Assembly/assembly.fasta",
    conda:
        "Envs/unicycler_4_8.yml"
    threads: 4
    shell:
        "unicycler -t {threads}   "  # --vcf
        " -1 {input.fq1_trimmed} -2 {input.fq2_trimmed} "
        " -o {output_Dir}/{wildcards.sampleID_WiIllumina}/IlluminaWGS/Unicycler_SPAdes_Assembly/ "
        
        # " ‑‑mode conservative " # This version of unicycler doesn't support the --mode arguement


rule rename_SPAdesAssembly_WithBioAWK_And_IndexFA:
    input:
        Ill_SPAdes_assembly_fa = output_Dir + "/{sampleID_WiIllumina}/IlluminaWGS/Unicycler_SPAdes_Assembly/assembly.fasta",
    output:
        Ill_SPAdes_Assembly_Renamed_fa = output_Dir + "/{sampleID_WiIllumina}/IlluminaWGS/Unicycler_SPAdes_Assembly/{sampleID_WiIllumina}.SPAdes.Assembly.fasta",
        Ill_SPAdes_Assembly_Renamed_fa_fai = output_Dir + "/{sampleID_WiIllumina}/IlluminaWGS/Unicycler_SPAdes_Assembly/{sampleID_WiIllumina}.SPAdes.Assembly.fasta.fai",

    shell:
        "bioawk -c fastx '{{ print \">{wildcards.sampleID_WiIllumina}_\"$name \"\\n\" $seq }}' {input} > {output.Ill_SPAdes_Assembly_Renamed_fa} \n"
        "samtools faidx {output.Ill_SPAdes_Assembly_Renamed_fa}"






##############################################################################
############ MM2: Assembly To H37rv Alignment & Variant Calling ##############
##############################################################################


MinAlnLen_ForCoverage = 5000  # The minimum length of an alignment to call a 
MinAlnLen_ForVariantCalling = 5000


rule Minimap2_Ill_SPAdes_Assembly_AlignTo_Ref:
    input:
        Ill_SPAdes_Assembly_Renamed_fa = output_Dir + "/{sampleID_WiIllumina}/IlluminaWGS/Unicycler_SPAdes_Assembly/{sampleID_WiIllumina}.SPAdes.Assembly.fasta",
        Ref_FA = refGenome_FA_PATH,
    output:
        MM2_Ill_SPAdes_Assembly_To_Ref_SAM = output_Dir + "/{sampleID_WiIllumina}/VariantCalling/Minimap2_Alignment/{sampleID_WiIllumina}.minimap2.paftools.sam",
        MM2_Ill_SPAdes_Assembly_To_Ref_BAM = output_Dir + "/{sampleID_WiIllumina}/VariantCalling/Minimap2_Alignment/{sampleID_WiIllumina}.minimap2.paftools.bam",
        MM2_Ill_SPAdes_Assembly_To_Ref_bai = output_Dir + "/{sampleID_WiIllumina}/VariantCalling/Minimap2_Alignment/{sampleID_WiIllumina}.minimap2.paftools.bam.bai",
        MM2_Ill_SPAdes_Assembly_To_Ref_VCF = output_Dir + "/{sampleID_WiIllumina}/VariantCalling/Minimap2_Alignment/{sampleID_WiIllumina}.minimap2.paftools.vcf",
        MM2_Ill_SPAdes_Assembly_To_Ref_PAF = output_Dir + "/{sampleID_WiIllumina}/VariantCalling/Minimap2_Alignment/{sampleID_WiIllumina}.minimap2.paftools.paf",

    #conda:"envs/PacBio_Software_py27_Conda.yml"

    threads: 1
    shell: 
        "minimap2 -ax asm10 --cs {input.Ref_FA} {input.Ill_SPAdes_Assembly_Renamed_fa} | awk '$5 != 0' > {output.MM2_Ill_SPAdes_Assembly_To_Ref_SAM} \n"
        "samtools view -bS {output.MM2_Ill_SPAdes_Assembly_To_Ref_SAM} | samtools sort - > {output.MM2_Ill_SPAdes_Assembly_To_Ref_BAM} \n"
        "samtools index {output.MM2_Ill_SPAdes_Assembly_To_Ref_BAM} \n"
        "minimap2 -cx asm10 --cs {input.Ref_FA} {input.Ill_SPAdes_Assembly_Renamed_fa} | awk '$5 != 0' | sort -k6,6 -k8,8n | paftools.js call -s {wildcards.sampleID_WiIllumina} -L {MinAlnLen_ForVariantCalling} -l {MinAlnLen_ForCoverage} -f {input.Ref_FA} - > {output.MM2_Ill_SPAdes_Assembly_To_Ref_VCF} \n"
        "minimap2 -cx asm10 --cs {input.Ref_FA} {input.Ill_SPAdes_Assembly_Renamed_fa} | awk '$5 != 0' | sort -k6,6 -k8,8n | paftools.js call -s {wildcards.sampleID_WiIllumina} -L {MinAlnLen_ForVariantCalling} -l {MinAlnLen_ForCoverage} - > {output.MM2_Ill_SPAdes_Assembly_To_Ref_PAF} \n"






rule NucDiff_Ill_SPAdes_Assembly_AlignTo_Ref_WithVCFoutput:
    input:
        Ill_SPAdes_Assembly_Renamed_fa = output_Dir + "/{sampleID_WiIllumina}/IlluminaWGS/Unicycler_SPAdes_Assembly/{sampleID_WiIllumina}.SPAdes.Assembly.fasta",
        Ref_FA = refGenome_FA_PATH,
    output:
        NucDiff_SmallVariants_GFF = output_Dir + "/{sampleID_WiIllumina}/VariantCalling/NucDiff_Analysis_{sampleID_WiIllumina}_V2_WiVCFout/results/{sampleID_WiIllumina}.NucDiff_ref_snps.gff",
        NucDiff_SmallVariants_VCF = output_Dir + "/{sampleID_WiIllumina}/VariantCalling/NucDiff_Analysis_{sampleID_WiIllumina}_V2_WiVCFout/results/{sampleID_WiIllumina}.NucDiff_ref_snps.vcf",
        NucDiff_Struct_GFF = output_Dir + "/{sampleID_WiIllumina}/VariantCalling/NucDiff_Analysis_{sampleID_WiIllumina}_V2_WiVCFout/results/{sampleID_WiIllumina}.NucDiff_ref_struct.gff",
        MUMmer_Delta_Aln_File = output_Dir + "/{sampleID_WiIllumina}/VariantCalling/NucDiff_Analysis_{sampleID_WiIllumina}_V2_WiVCFout/{sampleID_WiIllumina}.NucDiff.delta",
    conda:
        "Envs/nucdiff_2_0_3_Conda.yml"
    threads: 1
    shell:
        "nucdiff {input.Ref_FA} {input.Ill_SPAdes_Assembly_Renamed_fa} {output_Dir}/{wildcards.sampleID_WiIllumina}/VariantCalling/NucDiff_Analysis_{wildcards.sampleID_WiIllumina}_V2_WiVCFout {wildcards.sampleID_WiIllumina}.NucDiff --vcf yes \n"



rule SVanalyzer_SVrefine_Ill_SPAdes_Assembly_AlignTo_Ref:
    input:
        MUMmer_Delta_Aln_File = output_Dir + "/{sampleID_WiIllumina}/VariantCalling/NucDiff_Analysis_{sampleID_WiIllumina}_V2_WiVCFout/{sampleID_WiIllumina}.NucDiff.delta",
        Ill_SPAdes_Assembly_Renamed_fa = output_Dir + "/{sampleID_WiIllumina}/IlluminaWGS/Unicycler_SPAdes_Assembly/{sampleID_WiIllumina}.SPAdes.Assembly.fasta",
        Ref_FA = refGenome_FA_PATH,
    output:
        SVrefine_Ill_SPAdes_Assembly_VCF = output_Dir + "/{sampleID_WiIllumina}/VariantCalling/SVanalyzer_SVrefine_Ill_SPAdes_Assembly_SV_Calling/{sampleID_WiIllumina}.MUMmer.SVrefine.vcf"
    threads: 1
    shell:
        "SVrefine --delta {input.MUMmer_Delta_Aln_File} "
        "--ref_fasta {input.Ref_FA} "
        "--query_fasta {input.Ill_SPAdes_Assembly_Renamed_fa} "
        "--outvcf {output.SVrefine_Ill_SPAdes_Assembly_VCF} "
        " --samplename {wildcards.sampleID_WiIllumina} "




rule reference_faidx:
  input:
    REFERENCE = refGenome_FA_PATH
  output:
    refGenome_FA_PATH +".fai"
  shell:
    """
    samtools faidx {input}
    """

rule bwa_index:
  input:
    REFERENCE = refGenome_FA_PATH
  output:
    refGenome_FA_PATH +".bwt"
  shell:
    """
    bwa index {input}
    """


# MAPPING

rule bwa:
  input:
    fq1_trimmed = output_Dir + "/{sampleID_WiIllumina}/IlluminaWGS/FASTQs/{sampleID_WiIllumina}_1_trimmed.fastq",
    fq2_trimmed = output_Dir + "/{sampleID_WiIllumina}/IlluminaWGS/FASTQs/{sampleID_WiIllumina}_2_trimmed.fastq",
    reference = refGenome_FA_PATH +".bwt",
  output:
    sam = temp(output_Dir + "results/alignments/{sampleID_WiIllumina}.bwa.sam")
  params:
    reference = refGenome_FA_PATH
  threads: 5
  shell:
    """
    bwa mem -t {threads} {params.reference} {input.fq1_trimmed} {input.fq2_trimmed} > {output.sam}
    """

rule sam2bam:
  input:
    sam = output_Dir + "results/alignments/{sampleID_WiIllumina}.bwa.sam"
  output:
    bam = temp(output_Dir +"results/alignments/{sampleID_WiIllumina}.bwa.bam")
  shell:
    """
    samtools view -Sb {input.sam} | samtools sort > {output.bam}
    """

# MAPPING POSTPROCESSING

rule read_groups:
  input:
    bam = output_Dir + "results/alignments/{sampleID_WiIllumina}.bwa.bam"
  output:
    bam = temp(output_Dir +"results/alignments/{sampleID_WiIllumina}.bwa.rg.bam")
  shell:
    """
    picard AddOrReplaceReadGroups I={input.bam} O={output.bam} RGID=1 RGLB=lib RGPL=ILLUMINA RGPU=unit1 RGSM=Sample1
    """

rule deduplicate:
  input:
    bam = output_Dir + "results/alignments/{sampleID_WiIllumina}.bwa.rg.bam"
  output:
    bam = output_Dir + "results/alignments/{sampleID_WiIllumina}.bwa.rg.dedup.bam"
  params:
    txt = output_Dir + "results/alignments/{sampleID_WiIllumina}.bwa.rg.dedup.txt"

  shell:
    """
    picard MarkDuplicates REMOVE_DUPLICATES=true ASSUME_SORTED=true I={input.bam} O={output.bam} M={params.txt}
    """

rule bam_index:
  input:
    bam = output_Dir + "results/alignments/{sampleID_WiIllumina}.bwa.rg.dedup.bam"
  output:
    bam_index = output_Dir + "results/alignments/{sampleID_WiIllumina}.bwa.rg.dedup.bam.bai"
  shell:
    """
    samtools index {input.bam}
    """

# SV calling with Delly

rule delly:
  input:
    bam = output_Dir + "results/alignments/{sampleID_WiIllumina}.bwa.rg.dedup.bam",
    bai = output_Dir + "results/alignments/{sampleID_WiIllumina}.bwa.rg.dedup.bam.bai",
    ref_FA = refGenome_FA_PATH
  output:
    vcf = output_Dir + "results/variants_delly/{sampleID_WiIllumina}.bwa.rg.dedup.delly.vcf"
  params:
    bcf = output_Dir + "results/variants_delly/{sampleID_WiIllumina}.bwa.rg.dedup.delly.bcf"
  conda:
    "Envs/SV_Hack_V1_ShortRead_SV_Calling.yml"
  shell:
    """
    delly call -g {input.ref_FA} {input.bam} -o {params.bcf}
    bcftools view {params.bcf} > {output.vcf}
    """

# SV calling with Manta

rule manta:
  input:
    bam = output_Dir + "results/alignments/{sampleID_WiIllumina}.bwa.rg.dedup.bam",
    bai = output_Dir + "results/alignments/{sampleID_WiIllumina}.bwa.rg.dedup.bam.bai",
    ref_FA = refGenome_FA_PATH
  output:
    output_Dir + "results/variants_manta/{sampleID_WiIllumina}/results/variants/tumorSV.vcf"
  params:
    outdir = output_Dir + "results/variants_manta/{sampleID_WiIllumina}",
    output = output_Dir + "results/variants_manta/{sampleID_WiIllumina}/results/variants/tumorSV.vcf.gz"
  conda:
    "Envs/SV_Hack_V1_ShortRead_SV_Calling.yml"
  shell:
    """
    configManta.py --tumorBam {input.bam} --referenceFasta {input.ref_FA} --runDir {params.outdir}
    {params.outdir}/runWorkflow.py
    gunzip {params.output}
    """

# SV calling with Lumpy

rule lumpy:
  input:
    bam = output_Dir + "results/alignments/{sampleID_WiIllumina}.bwa.rg.dedup.bam"
  output:
    output_Dir + "results/variants_lumpy/{sampleID_WiIllumina}.bwa.rg.dedup.lumpy.vcf"
  conda:
    "Envs/SV_Hack_V1_ShortRead_SV_Calling.yml"
  shell:
    """
    lumpyexpress \
    -B {input.bam} \
    -o {output}
    """ 
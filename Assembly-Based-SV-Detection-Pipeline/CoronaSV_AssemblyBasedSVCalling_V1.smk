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


SampleIDTo_Nanopore_FQ_Dict = {}

SampleIDTo_Illumina_FQ1_Dict = {}
SampleIDTo_Illumina_FQ2_Dict = {}


for idx, row in df.iterrows():
    
    
    i_SampleID = row["SampleID"]
    
        
    i_Nanopore_FQ = row["OxfordNanopore_FQ_PATH"]


    if i_Nanopore_FQ != "None":
        SampleIDTo_Nanopore_FQ_Dict[i_SampleID] = i_Nanopore_FQ



    i_FastQ_Files = row["IlluminaPE_FQs"]
    i_FastQ_Files_List = i_FastQ_Files.split(";")
    
    if len(i_FastQ_Files_List) == 2:
        FQ_1_PATH, FQ_2_PATH = i_FastQ_Files_List

        SampleIDTo_Illumina_FQ1_Dict[i_SampleID] = FQ_1_PATH
        SampleIDTo_Illumina_FQ2_Dict[i_SampleID] = FQ_2_PATH
    

input_SampleIDs_WiNanopore = list(SampleIDTo_Nanopore_FQ_Dict.keys())
input_SampleIDs_WiIllumina = list(SampleIDTo_Illumina_FQ1_Dict.keys())

#print("SampleIDs with ONT data:", input_SampleIDs_WiNanopore)
#print("SampleIDs with Illumina data:", input_SampleIDs_WiIllumina)

#print(SampleIDTo_Nanopore_FQ_Dict)


rule all:
    input:
        expand(output_Dir + "/{sampleID_WiNanopore}/Nanopore/Nanopore_QC/{sampleID_WiNanopore}.NanoPlot/NanoPlot-report.html", sampleID_WiNanopore=input_SampleIDs_WiNanopore),
        expand(output_Dir + "/{sampleID_WiNanopore}/Nanopore/Flye_Assembly/{sampleID_WiNanopore}.Flye.Assembly.fasta.fai", sampleID_WiNanopore=input_SampleIDs_WiNanopore),
        #expand(output_Dir + "/{sampleID_WiNanopore}/Nanopore/Flye_Assembly_Meta/assembly.fasta.fai", sampleID_WiNanopore=input_SampleIDs_WiNanopore),
        #expand(output_Dir + "/{sampleID_WiNanopore}/Nanopore/Flye_Assembly_KeepHaplotypes/assembly.fasta.fai", sampleID_WiNanopore=input_SampleIDs_WiNanopore),
        #expand(output_Dir + "/{sampleID_WiNanopore}/Nanopore/Unicycler_Assembly/assembly.fasta", sampleID_WiNanopore=input_SampleIDs_WiNanopore),
        
        expand(output_Dir + "/{sampleID_WiIllumina}/IlluminaWGS/Unicycler_SPAdes_Assembly/{sampleID_WiIllumina}.SPAdes.Assembly.fasta.fai", sampleID_WiIllumina=input_SampleIDs_WiIllumina),

        expand(output_Dir + "/{sampleID_WiIllumina}/VariantCalling/Minimap2_GC3_PP_AlignTo_H37rv/{sampleID_WiIllumina}_mm2_GC3_PP_AssemblyToH37rv.vcf", sampleID_WiIllumina=input_SampleIDs_WiIllumina),

        expand(output_Dir + "/{sampleID_WiIllumina}/VariantCalling/NucDiff_Analysis_{sampleID_WiIllumina}_V2_WiVCFout/results/NucDiff_{sampleID_WiIllumina}_ref_struct.Filtered.SVs.gff", sampleID_WiIllumina=input_SampleIDs_WiIllumina),




rule nanoplot_QC:
    input:
        ONT_reads_fq = lambda wildcards: SampleIDTo_Nanopore_FQ_Dict[wildcards.sampleID_WiNanopore],
    output:
         output_Dir + "/{sampleID_WiNanopore}/Nanopore/Nanopore_QC/{sampleID_WiNanopore}.NanoPlot/NanoPlot-report.html"
    conda:
        "Envs/nanoplot_128_Conda.yml"
    threads: 2
    shell:
        "NanoPlot -t {threads} --fastq {input.ONT_reads_fq} -o {output_Dir}/{wildcards.sampleID_WiNanopore}/Nanopore/Nanopore_QC/{wildcards.sampleID_WiNanopore}.NanoPlot/"



### A) Assembly with SPAdes through Unicycler
rule unicycler_Assemble_Nanopore_WGS:
    input:
        ONT_reads_fq = lambda wildcards: SampleIDTo_Nanopore_FQ_Dict[wildcards.sampleID_WiNanopore],
    output:
        assembly_GFA = output_Dir + "/{sampleID_WiNanopore}/Nanopore/Unicycler_Assembly/assembly.gfa",
        assembly_fa = output_Dir + "/{sampleID_WiNanopore}/Nanopore/Unicycler_Assembly/assembly.fasta",
    conda:
        "Envs/unicycler_4_8.yml"
    threads: 1
    shell:
        "unicycler --vcf -t {threads}   "
        " -l {input.ONT_reads_fq}"
        " -o {output_Dir}/{wildcards.sampleID_WiNanopore}/Nanopore/Unicycler_Assembly/ "
        
        # " ‑‑mode conservative " # This version of unicycler doesn't support the --mode arguement







###################################################
##### Flye (correct reads & create assembly) ######
###################################################

rule flye_Assemble:
    input:
        ONT_reads_fq = lambda wildcards: SampleIDTo_Nanopore_FQ_Dict[wildcards.sampleID_WiNanopore],
    output:
        ONT_Flye_Assembly_fa = output_Dir + "/{sampleID_WiNanopore}/Nanopore/Flye_Assembly/assembly.fasta",
        ONT_Flye_Assembly_info_txt = output_Dir + "/{sampleID_WiNanopore}/Nanopore/Flye_Assembly/assembly_info.txt"
    conda:
        "Envs/flye_2_8_Conda.yml"
    threads: 2
    shell:
        "flye --nano-raw {input.ONT_reads_fq} --out-dir {output_Dir}/{wildcards.sampleID_WiNanopore}/Nanopore/Flye_Assembly/  --threads {threads}"



rule rename_FlyeAssembly_WithBioAWK:
    input:
        ONT_Flye_Assembly_fa = output_Dir + "/{sampleID_WiNanopore}/Nanopore/Flye_Assembly/assembly.fasta",
    output:
        ONT_Flye_Assembly_Renamed_fa = output_Dir + "/{sampleID_WiNanopore}/Nanopore/Flye_Assembly/{sampleID_WiNanopore}.Flye.Assembly.fasta",

    shell:
        "bioawk -c fastx '{{ print \">{wildcards.sampleID_WiNanopore}_\"$name \"\\n\" $seq }}' {input} > {output}"



rule samtools_faidx_FlyeAssembly:
    input:
        ONT_Flye_Assembly_Renamed_fa = output_Dir + "/{sampleID_WiNanopore}/Nanopore/Flye_Assembly/{sampleID_WiNanopore}.Flye.Assembly.fasta",
    output:
        output_Dir + "/{sampleID_WiNanopore}/Nanopore/Flye_Assembly/{sampleID_WiNanopore}.Flye.Assembly.fasta.fai"
    conda:
        "Envs/samtools_AND_bcftools_200128_Conda.yml"
    threads: 1
    shell: "samtools faidx {input}"


rule flye_Assemble_Meta:
    input:
        ONT_reads_fq = lambda wildcards: SampleIDTo_Nanopore_FQ_Dict[wildcards.sampleID_WiNanopore],
    output:
        ONT_Flye_Assembly_fa = output_Dir + "/{sampleID_WiNanopore}/Nanopore/Flye_Assembly_Meta/assembly.fasta",
        ONT_Flye_Assembly_info_txt = output_Dir + "/{sampleID_WiNanopore}/Nanopore/Flye_Assembly_Meta/assembly_info.txt"
    conda:
        "Envs/flye_2_8_Conda.yml"
    threads: 2
    shell:
        "flye --nano-raw {input.ONT_reads_fq} --out-dir {output_Dir}/{wildcards.sampleID_WiNanopore}/Nanopore/Flye_Assembly_Meta/ --meta --threads {threads}"



rule samtools_faidx_FlyeAssembly_Meta:
    input:
        output_Dir + "/{sampleID_WiNanopore}/Nanopore/Flye_Assembly_Meta/assembly.fasta"
    output:
        output_Dir + "/{sampleID_WiNanopore}/Nanopore/Flye_Assembly_Meta/assembly.fasta.fai"
    conda:
        "Envs/samtools_AND_bcftools_200128_Conda.yml"
    threads: 1
    shell: "samtools faidx {input}"



rule flye_Assemble_KeepHaplotypes:
    input:
        ONT_reads_fq = lambda wildcards: SampleIDTo_Nanopore_FQ_Dict[wildcards.sampleID_WiNanopore],
    output:
        ONT_Flye_Assembly_fa = output_Dir + "/{sampleID_WiNanopore}/Nanopore/Flye_Assembly_KeepHaplotypes/assembly.fasta",
        ONT_Flye_Assembly_info_txt = output_Dir + "/{sampleID_WiNanopore}/Nanopore/Flye_Assembly_KeepHaplotypes/assembly_info.txt"
    conda:
        "Envs/flye_2_8_Conda.yml"
    threads: 2
    shell:
        "flye --nano-raw {input.ONT_reads_fq} --out-dir {output_Dir}/{wildcards.sampleID_WiNanopore}/Nanopore/Flye_Assembly_KeepHaplotypes/ --keep-haplotypes --threads {threads}"


rule samtools_faidx_FlyeAssembly_KeepHaplotypes:
    input:
        output_Dir + "/{sampleID_WiNanopore}/Nanopore/Flye_Assembly_KeepHaplotypes/assembly.fasta"
    output:
        output_Dir + "/{sampleID_WiNanopore}/Nanopore/Flye_Assembly_KeepHaplotypes/assembly.fasta.fai"
    conda:
        "Envs/samtools_AND_bcftools_200128_Conda.yml"
    threads: 1
    shell: "samtools faidx {input}"







#SampleIDTo_Illumina_FQ1_Dict[i_SampleID] = FQ_1_PATH
#        SampleIDTo_Illumina_FQ2_Dict[i_SampleID] = FQ_2_PATH
    
# lambda wildcards: SampleIDTo_Nanopore_FQ_Dict[wildcards.sampleID_WiNanopore],



# adapter list from: https://github.com/stephenturner/adapters/blob/master/adapters_combined_256_unique.fasta
# can move adapter list to config file later
rule trimmomatic_Illumina_PE_Trimming:
    input:
        r1 = lambda wildcards: SampleIDTo_Illumina_FQ1_Dict[wildcards.sampleID_WiIllumina],
        r2 = lambda wildcards: SampleIDTo_Illumina_FQ2_Dict[wildcards.sampleID_WiIllumina],
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
    threads: 8
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
    threads: 2
    shell:
        "unicycler --vcf -t {threads}   "
        " -1 {input.fq1_trimmed} -2 {input.fq2_trimmed} "
        " -o {output_Dir}/{wildcards.sampleID_WiIllumina}/IlluminaWGS/Unicycler_SPAdes_Assembly/ "
        
        # " ‑‑mode conservative " # This version of unicycler doesn't support the --mode arguement


rule rename_SPAdesAssembly_WithBioAWK:
    input:
        Ill_SPAdes_assembly_fa = output_Dir + "/{sampleID_WiIllumina}/IlluminaWGS/Unicycler_SPAdes_Assembly/assembly.fasta",
    output:
        Ill_SPAdes_Assembly_Renamed_fa = output_Dir + "/{sampleID_WiIllumina}/IlluminaWGS/Unicycler_SPAdes_Assembly/{sampleID_WiIllumina}.SPAdes.Assembly.fasta",

    shell:
        "bioawk -c fastx '{{ print \">{wildcards.sampleID_WiIllumina}_\"$name \"\\n\" $seq }}' {input} > {output}"



rule samtools_faidx_SPAdesAssembly:
    input:
        Ill_SPAdes_Assembly_Renamed_fa = output_Dir + "/{sampleID_WiIllumina}/IlluminaWGS/Unicycler_SPAdes_Assembly/{sampleID_WiIllumina}.SPAdes.Assembly.fasta",
    output:
        output_Dir + "/{sampleID_WiNanopore}/IlluminaWGS/Unicycler_SPAdes_Assembly/{sampleID_WiIllumina}.SPAdes.Assembly.fasta.fai"
    conda:
        "Envs/samtools_AND_bcftools_200128_Conda.yml"
    threads: 1
    shell: "samtools faidx {input}"











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
        MM2_Ill_SPAdes_Assembly_To_Ref_SAM = output_Dir + "/{sampleID_WiIllumina}/VariantCalling/Minimap2_GC3_PP_AlignTo_H37rv/{sampleID_WiIllumina}_mm2_GC3_PP_AssemblyToH37rv.sam",
        MM2_Ill_SPAdes_Assembly_To_Ref_BAM = output_Dir + "/{sampleID_WiIllumina}/VariantCalling/Minimap2_GC3_PP_AlignTo_H37rv/{sampleID_WiIllumina}_mm2_GC3_PP_AssemblyToH37rv.bam",
        MM2_Ill_SPAdes_Assembly_To_Ref_bai = output_Dir + "/{sampleID_WiIllumina}/VariantCalling/Minimap2_GC3_PP_AlignTo_H37rv/{sampleID_WiIllumina}_mm2_GC3_PP_AssemblyToH37rv.bam.bai",
        MM2_Ill_SPAdes_Assembly_To_Ref_VCF = output_Dir + "/{sampleID_WiIllumina}/VariantCalling/Minimap2_GC3_PP_AlignTo_H37rv/{sampleID_WiIllumina}_mm2_GC3_PP_AssemblyToH37rv.vcf",
        MM2_Ill_SPAdes_Assembly_To_Ref_PAF = output_Dir + "/{sampleID_WiIllumina}/VariantCalling/Minimap2_GC3_PP_AlignTo_H37rv/{sampleID_WiIllumina}_mm2_GC3_PP_AssemblyToH37rv.paf",

    #conda:"envs/PacBio_Software_py27_Conda.yml"

    threads: 1
    shell: 
        "minimap2 -ax asm10 --cs {input.Ref_FA} {input.Ill_SPAdes_Assembly_Renamed_fa} | awk '$5 != 0' > {output.MM2_Ill_SPAdes_Assembly_To_Ref_SAM} \n"
        "samtools view -bS {output.MM2_Ill_SPAdes_Assembly_To_Ref_SAM} | samtools sort - > {output.MM2_Ill_SPAdes_Assembly_To_Ref_BAM} \n"
        "samtools index {output.MM2_Ill_SPAdes_Assembly_To_Ref_BAM} \n"
        "minimap2 -cx asm10 --cs {input.Ref_FA} {input.Ill_SPAdes_Assembly_Renamed_fa} | awk '$5 != 0' | sort -k6,6 -k8,8n | paftools.js call -s {wildcards.sampleID_WiIllumina} -L {MinAlnLen_ForVariantCalling} -l {MinAlnLen_ForCoverage} -f {input.Ref_FA} - > {output.MM2_Ill_SPAdes_Assembly_To_Ref_VCF} \n"
        "minimap2 -cx asm10 --cs {input.Ref_FA} {input.Ill_SPAdes_Assembly_Renamed_fa} | awk '$5 != 0' | sort -k6,6 -k8,8n | paftools.js call -s {wildcards.sampleID_WiIllumina} -L {MinAlnLen_ForVariantCalling} -l {MinAlnLen_ForCoverage} - > {output.MM2_Ill_SPAdes_Assembly_To_Ref_PAF} \n"






rule NucDiff_Analysis_G3_PP_vs_H37rv_WithVCFoutput:
    input:
        Ill_SPAdes_Assembly_Renamed_fa = output_Dir + "/{sampleID_WiIllumina}/IlluminaWGS/Unicycler_SPAdes_Assembly/{sampleID_WiIllumina}.SPAdes.Assembly.fasta",
        Ref_FA = refGenome_FA_PATH,
    output:
        NucDiff_SmallVariants_GFF = output_Dir + "/{sampleID_WiIllumina}/VariantCalling/NucDiff_Analysis_{sampleID_WiIllumina}_V2_WiVCFout/results/NucDiff_{sampleID_WiIllumina}_ref_snps.gff",
        NucDiff_SmallVariants_VCF = output_Dir + "/{sampleID_WiIllumina}/VariantCalling/NucDiff_Analysis_{sampleID_WiIllumina}_V2_WiVCFout/results/NucDiff_{sampleID_WiIllumina}_ref_snps.vcf",
        NucDiff_Struct_GFF = output_Dir + "/{sampleID_WiIllumina}/VariantCalling/NucDiff_Analysis_{sampleID_WiIllumina}_V2_WiVCFout/results/NucDiff_{sampleID_WiIllumina}_ref_struct.gff",
        NucDiff_Struct_FilteredSVs_GFF = output_Dir + "/{sampleID_WiIllumina}/VariantCalling/NucDiff_Analysis_{sampleID_WiIllumina}_V2_WiVCFout/results/NucDiff_{sampleID_WiIllumina}_ref_struct.Filtered.SVs.gff",
    conda:
        "Envs/nucdiff_2_0_3_Conda.yml"
    threads: 1
    shell:
        "nucdiff {input.Ref_FA} {input.Ill_SPAdes_Assembly_Renamed_fa} {output_Dir}/{wildcards.sampleID_WiIllumina}/VariantCalling/NucDiff_Analysis_{wildcards.sampleID_WiIllumina}_V2_WiVCFout NucDiff_{wildcards.sampleID_WiIllumina} --vcf yes \n"
        ' grep -v "reshuffling" {output.NucDiff_Struct_GFF} > {output.NucDiff_Struct_FilteredSVs_GFF} \n'















_____# This is a snakemake file/script for SV Detection within the SARS-COV2 genome using an assembly based approach

### Input Data Types:

### Maximillian Marin (mgmarin@g.harvard.edu)


# Define PATH to the reference genome to be used: 
refGenome_FA_PATH = config["RefGenome_FA_PATH"]
refGenome_GFF_PATH = config["RefGenome_GFF_PATH"]




# Define PATH of OUTPUT directory

output_Dir = config["output_dir"]



# Read in meta data regarding input 
df = pd.read_csv( config["inputSampleData_TSV"], sep=',')






SampleIDTo_Nanopore_FQ_Dict = {}

SampleIDTo_Illumina_FQ1_Dict = {}
SampleIDTo_Illumina_FQ2_Dict = {}



for idx, row in df.iterrows():
    
    
    i_SampleID = row["SampleID"]
    
        
    i_Nanopore_FQ = row["OxfordNanopore_FQ_PATH"]

    SampleIDTo_Nanopore_FQ_Dict = i_Nanopore_FQ


    i_FastQ_Files = row["IlluminaPE_FQs"]
    i_FastQ_Files_List = FastQ_Files.split(";")
    
    if len(i_FastQ_Files_List) == 2:
        FQ_1_PATH, FQ_2_PATH = i_FastQ_Files_List

        SampleIDTo_FQ1_Dict[SampleID_i] = FQ_1_PATH
        SampleIDTo_FQ2_Dict[SampleID_i] = FQ_2_PATH
    



input_SampleIDs_WiNanopore = list(SampleIDTo_Nanopore_FQ_Dict.keys())
input_SampleIDs_WiIllumina = list(SampleIDTo_FQ1_Dict.keys())




rule all:
    input:



rule nanoplot_QC:
    input:
        lambda wildcards: SampleIDTo_Nanopore_FQ_Dict[wildcards.RunID],
    output:
        output_Dir + "{sampleID}/Nanopore/Nanoplot_QC/{sampleID}.subreads.NanoPlot/NanoPlot-report.html"
    conda:
        "envs/nanoplot_128_Conda.yml"
    threads: 4
    shell: "NanoPlot -t {threads} --fastq {input} -o {output_Dir}{wildcards.sampleID}/Nanopore/Nanoplot_QC/{wildcards.sampleID}.subreads.NanoPlot/"



###################################################
##### Flye (correct reads & create assembly) ######
###################################################

rule flye_Assemble:
    input:
        ONT_reads_fq = lambda wildcards: SampleIDTo_Nanopore_FQ_Dict[wildcards.RunID],
    output:
        assembly_fa = output_Dir + "{sampleID}/Nanopore/Flye_Assembly/assembly.fasta",
        assembly_info_txt = output_Dir + "{sampleID}/Nanopore/Flye_Assembly/assembly_info.txt"
    conda:
        "envs/flye_2_8_Conda.yml"
    threads: 10
    shell:
        "flye --nano-raw {input.ONT_reads_fq} --out-dir {output_Dir}{wildcards.sampleID}/Nanopore/Flye_Assembly/ --genome-size 5m --threads {threads}"


rule samtools_faidx_FlyeAssembly:
    input:
        output_Dir + "{sampleID}/Nanopore/Flye_Assembly/{sampleID}.flyeassembly.fixstart.fasta"
    output:
        output_Dir + "{sampleID}/Nanopore/Flye_Assembly/{sampleID}.flyeassembly.fixstart.fasta.fai"
    conda:
        "envs/PacBio_Software_py27_Conda.yml"
    threads: 1
    shell: "samtools faidx {input}"













# adapter list from: https://github.com/stephenturner/adapters/blob/master/adapters_combined_256_unique.fasta
# can move adapter list to config file later
rule trimmomatic_Illumina_PE_Trimming:
    input:
         r1 = output_Dir + "{sampleID_WiIll}/IlluminaWGS/FASTQs/{sampleID_WiIll}_1.fastq",
         r2 = output_Dir + "{sampleID_WiIll}/IlluminaWGS/FASTQs/{sampleID_WiIll}_2.fastq",
    output:
        r1 = output_Dir + "{sampleID_WiIll}/IlluminaWGS/FASTQs/{sampleID_WiIll}_1_trimmed.fastq",
        r2 = output_Dir + "{sampleID_WiIll}/IlluminaWGS/FASTQs/{sampleID_WiIll}_2_trimmed.fastq",
        # reads where trimming entirely removed the mate
        r1_unpaired = output_Dir + "{sampleID_WiIll}/IlluminaWGS/FASTQs/{sampleID_WiIll}_1_trimmed.unpaired.fastq",
        r2_unpaired = output_Dir + "{sampleID_WiIll}/IlluminaWGS/FASTQs/{sampleID_WiIll}_2_trimmed.unpaired.fastq",
    log:
        "logs/trimmomatic/{sampleID_WiIll}.log"
    params:
        # list of trimmers (see manual)
        trimmer=["ILLUMINACLIP:'references/CustomTrimmoatic_IlluminaWGS_AdapterList.fasta':2:30:10:2:true SLIDINGWINDOW:4:20 MINLEN:75"],
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
        fq1_trimmed = output_Dir + "{sampleID_WiIll}/IlluminaWGS/FQs_Trimmomatic_Trimming/{sampleID_WiIll}_1_trimmed.fastq",
        fq2_trimmed = output_Dir + "{sampleID_WiIll}/IlluminaWGS/FQs_Trimmomatic_Trimming/{sampleID_WiIll}_2_trimmed.fastq",
    output:
        assembly_GFA = output_Dir + "{sampleID_WiIll}/IlluminaWGS/Unicycler_SPAdesAssembly/assembly.gfa",
        assembly_fa = output_Dir + "{sampleID_WiIll}/IlluminaWGS/Unicycler_SPAdesAssembly/assembly.fasta",
    conda:
        "envs/unicycler_4_8.yml"
    threads: 8
    shell:
        "unicycler --vcf -t {threads}   "
        " -1 {input.fq1_trimmed} -2 {input.fq2_trimmed} "
        " -o {output_Dir}{wildcards.sampleID_WiIll}/IlluminaWGS/Unicycler_SPAdesAssembly/ "
        
        # " ‑‑mode conservative " # This version of unicycler doesn't support the --mode arguement















##############################################################################
############ MM2: Assembly To H37rv Alignment & Variant Calling ##############
##############################################################################


MinAlnLen_ForCoverage = 5000  # The minimum length of an alignment to call a 
MinAlnLen_ForVariantCalling = 5000


rule Minimap2_GC3_PP_AlignTo_H37rv:
    input:
        Q3Assembly_PilonPolished_FA = output_Dir+ "{sampleID_WiIll}/GC3_IlluminaPolishing/pilon_IllPE_Polishing_Q3_Assembly_ChangeSNPsINDELsOnly/{sampleID_WiIll}.IllPE.Q3Assembly.fasta",
        H37rv_FA = refGenome_FA_PATH,
    output:
        MM2_GC3_PP_To_H37rv_SAM = output_Dir + "{sampleID_WiIll}/VariantCallingVersusH37Rv/Minimap2_GC3_PP_AlignTo_H37rv/{sampleID_WiIll}_mm2_GC3_PP_AssemblyToH37rv.sam",
        MM2_GC3_PP_To_H37rv_BAM = output_Dir + "{sampleID_WiIll}/VariantCallingVersusH37Rv/Minimap2_GC3_PP_AlignTo_H37rv/{sampleID_WiIll}_mm2_GC3_PP_AssemblyToH37rv.bam",
        MM2_GC3_PP_To_H37rv_bai = output_Dir + "{sampleID_WiIll}/VariantCallingVersusH37Rv/Minimap2_GC3_PP_AlignTo_H37rv/{sampleID_WiIll}_mm2_GC3_PP_AssemblyToH37rv.bam.bai",
        MM2_GC3_PP_To_H37rv_VCF = output_Dir + "{sampleID_WiIll}/VariantCallingVersusH37Rv/Minimap2_GC3_PP_AlignTo_H37rv/{sampleID_WiIll}_mm2_GC3_PP_AssemblyToH37rv.vcf",
        MM2_GC3_PP_To_H37rv_PAF = output_Dir + "{sampleID_WiIll}/VariantCallingVersusH37Rv/Minimap2_GC3_PP_AlignTo_H37rv/{sampleID_WiIll}_mm2_GC3_PP_AssemblyToH37rv.paf",
    conda:
        "envs/PacBio_Software_py27_Conda.yml"
    threads: 1
    shell: 
        "minimap2 -ax asm10 --cs {input.H37rv_FA} {input.Q3Assembly_PilonPolished_FA} | awk '$5 != 0' > {output.MM2_GC3_PP_To_H37rv_SAM} \n"
        "samtools view -bS {output.MM2_GC3_PP_To_H37rv_SAM} | samtools sort - > {output.MM2_GC3_PP_To_H37rv_BAM} \n"
        "samtools index {output.MM2_GC3_PP_To_H37rv_BAM} \n"
        "minimap2 -cx asm10 --cs {input.H37rv_FA} {input.Q3Assembly_PilonPolished_FA} | awk '$5 != 0' | sort -k6,6 -k8,8n | paftools.js call -s {wildcards.sampleID_WiIll} -L {MinAlnLen_ForVariantCalling} -l {MinAlnLen_ForCoverage} -f {input.H37rv_FA} - > {output.MM2_GC3_PP_To_H37rv_VCF} \n"
        "minimap2 -cx asm10 --cs {input.H37rv_FA} {input.Q3Assembly_PilonPolished_FA} | awk '$5 != 0' | sort -k6,6 -k8,8n | paftools.js call -s {wildcards.sampleID_WiIll} -L {MinAlnLen_ForVariantCalling} -l {MinAlnLen_ForCoverage} - > {output.MM2_GC3_PP_To_H37rv_PAF} \n"






rule NucDiff_Analysis_G3_PP_vs_H37rv_WithVCFoutput:
    input:
        input_Coronavirus_Flye_Assembly_FA = output_Dir + "{sampleID_WiIll}/____/{sampleID_WiIll}.fasta",
        Reference_FA = refGenome_FA_PATH
    output:
        NucDiff_SmallVariants_GFF = output_Dir + "{sampleID_WiIll}/VariantCallingVersusH37Rv/NucDiff_Analysis_{sampleID_WiIll}_V2_WiVCFout/results/NucDiff_{sampleID_WiIll}_ref_snps.gff",
        NucDiff_SmallVariants_VCF = output_Dir + "{sampleID_WiIll}/VariantCallingVersusH37Rv/NucDiff_Analysis_{sampleID_WiIll}_V2_WiVCFout/results/NucDiff_{sampleID_WiIll}_ref_snps.vcf",
        NucDiff_Struct_GFF = output_Dir + "{sampleID_WiIll}/VariantCallingVersusH37Rv/NucDiff_Analysis_{sampleID_WiIll}_V2_WiVCFout/results/NucDiff_{sampleID_WiIll}_ref_struct.gff",
        NucDiff_Struct_FilteredSVs_GFF = output_Dir + "{sampleID_WiIll}/VariantCallingVersusH37Rv/NucDiff_Analysis_{sampleID_WiIll}_V2_WiVCFout/results/NucDiff_{sampleID_WiIll}_ref_struct.Filtered.SVs.gff",
    conda:
        "envs/nucdiff_2_0_3_Conda.yml"
    threads: 1
    shell:
        "nucdiff {input.Reference_FA} {input.input_Coronavirus_Flye_Assembly_FA} {output_Dir}{wildcards.sampleID_WiIll}/VariantCallingVersusH37Rv/NucDiff_Analysis_{wildcards.sampleID_WiIll}_V2_WiVCFout NucDiff_{wildcards.sampleID_WiIll} --vcf yes \n"
        ' grep -v "reshuffling" {output.NucDiff_GC3_PP_To_H37rv_Struct_GFF} > {output.NucDiff_GC3_PP_To_H37rv_Struct_FilteredSVs_GFF} \n'















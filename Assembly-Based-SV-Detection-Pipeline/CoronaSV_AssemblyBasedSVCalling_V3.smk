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

print(len(input_SampleIDs_WiIllumina))
print(len(input_SampleIDs_WiNanopore))


rule all:
    input:
        expand(output_Dir + "/Nanopore_FQs/{sampleID_WiNanopore}.fastq", sampleID_WiNanopore=input_SampleIDs_WiNanopore),
        expand(output_Dir + "/Illumina_PE_FQs/{sampleID_WiIllumina}_1.fastq", sampleID_WiIllumina=input_SampleIDs_WiIllumina),
        expand(output_Dir + "/Illumina_PE_FQs/{sampleID_WiIllumina}_2.fastq", sampleID_WiIllumina=input_SampleIDs_WiIllumina),
        #expand(output_Dir + "/{sampleID_WiIllumina}/IlluminaWGS/Unicycler_SPAdes_Assembly/{sampleID_WiIllumina}.SPAdes.Assembly.fasta.fai", sampleID_WiIllumina=input_SampleIDs_WiIllumina),
        #expand(output_Dir + "/{sampleID_WiIllumina}/VariantCalling/Minimap2_GC3_PP_AlignTo_H37rv/{sampleID_WiIllumina}_mm2_GC3_PP_AssemblyToH37rv.vcf", sampleID_WiIllumina=input_SampleIDs_WiIllumina),
        #expand(output_Dir + "/{sampleID_WiIllumina}/VariantCalling/NucDiff_Analysis_{sampleID_WiIllumina}_V2_WiVCFout/results/NucDiff_{sampleID_WiIllumina}_ref_struct.Filtered.SVs.gff", sampleID_WiIllumina=input_SampleIDs_WiIllumina),



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















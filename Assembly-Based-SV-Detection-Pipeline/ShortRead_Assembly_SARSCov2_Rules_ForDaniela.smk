



################################################################################
######### UniCycler (Using SPAdes Internally) for short read assembly  #########
################################################################################

# Adding SPADes assembly with UNIcycler to analysis

### A) Assembly with SPAdes through Unicycler
rule unicycler_SPAdes_Assemble_IlluminaWGS:
    input:
        fq1_trimmed = __________
        fq2_trimmed = __________
    output:
        assembly_GFA = "results/Assemblies/Unicycler_SPAdes_Assembly_{run}/assembly.gfa",
        Ill_SPAdes_assembly_fa = "results/Assemblies/Unicycler_SPAdes_Assembly_{run}/assembly.fasta",
    conda:
        "Envs/unicycler_4_8.yml"
    threads: 4
    shell:
        "unicycler -t {threads} " 
        " -1 {input.fq1_trimmed} -2 {input.fq2_trimmed} "
        " -o results/Assemblies/Unicycler_SPAdes_Assembly_{wildcards.run}/ "


rule rename_SPAdesAssembly_WithBioAWK:
    input:
        Ill_SPAdes_assembly_fa = "results/Assemblies/Unicycler_SPAdes_Assembly_{run}/assembly.fasta",
    output:
        Ill_SPAdes_Assembly_Renamed_fa = "results/Assemblies/Unicycler_SPAdes_Assembly_{run}/{run}.SPAdes.Assembly.fasta",
    conda:
        "Envs/SV_Hack_V1.yml"
    shell:
        "bioawk -c fastx '{{ print \">{wildcards.run}_\"$name \"\\n\" $seq }}' {input} > {output}"




rule samtools_faidx_SPAdesAssembly:
    input:
        Ill_SPAdes_Assembly_Renamed_fa = "results/Assemblies/Unicycler_SPAdes_Assembly_{run}/{run}.SPAdes.Assembly.fasta",
    output:
        "results/Assemblies/Unicycler_SPAdes_Assembly_{run}/{run}.SPAdes.Assembly.fasta.fai"
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
        Ill_SPAdes_Assembly_Renamed_fa = "results/Assemblies/Unicycler_SPAdes_Assembly_{run}/{run}.SPAdes.Assembly.fasta",
        Ref_FA = refGenome_FA_PATH,
    output:
        MM2_Ill_SPAdes_Assembly_To_Ref_SAM = "results/AssemblyToRef_Variants/Minimap2_Alignment/{run}.minimap2.paftools.sam",
        MM2_Ill_SPAdes_Assembly_To_Ref_BAM = "results/AssemblyToRef_Variants/Minimap2_Alignment/{run}.minimap2.paftools.bam",
        MM2_Ill_SPAdes_Assembly_To_Ref_bai = "results/AssemblyToRef_Variants/Minimap2_Alignment/{run}.minimap2.paftools.bam.bai",
        MM2_Ill_SPAdes_Assembly_To_Ref_VCF = "results/AssemblyToRef_Variants/Minimap2_Alignment/{run}.minimap2.paftools.vcf",
        MM2_Ill_SPAdes_Assembly_To_Ref_PAF = "results/AssemblyToRef_Variants/Minimap2_Alignment/{run}.minimap2.paftools.paf",
    threads: 1
    conda:
        "Envs/SV_Hack_V1.yml"
    shell: 
        "minimap2 -ax asm10 --cs {input.Ref_FA} {input.Ill_SPAdes_Assembly_Renamed_fa} | awk '$5 != 0' > {output.MM2_Ill_SPAdes_Assembly_To_Ref_SAM} \n"
        "samtools view -bS {output.MM2_Ill_SPAdes_Assembly_To_Ref_SAM} | samtools sort - > {output.MM2_Ill_SPAdes_Assembly_To_Ref_BAM} \n"
        "samtools index {output.MM2_Ill_SPAdes_Assembly_To_Ref_BAM} \n"
        "minimap2 -cx asm10 --cs {input.Ref_FA} {input.Ill_SPAdes_Assembly_Renamed_fa} | awk '$5 != 0' | sort -k6,6 -k8,8n | paftools.js call -s {wildcards.run} -L {MinAlnLen_ForVariantCalling} -l {MinAlnLen_ForCoverage} -f {input.Ref_FA} - > {output.MM2_Ill_SPAdes_Assembly_To_Ref_VCF} \n"
        "minimap2 -cx asm10 --cs {input.Ref_FA} {input.Ill_SPAdes_Assembly_Renamed_fa} | awk '$5 != 0' | sort -k6,6 -k8,8n | paftools.js call -s {wildcards.run} -L {MinAlnLen_ForVariantCalling} -l {MinAlnLen_ForCoverage} - > {output.MM2_Ill_SPAdes_Assembly_To_Ref_PAF} \n"





rule NucDiff_Ill_SPAdes_Assembly_AlignTo_Ref_WithVCFoutput:
    input:
        Ill_SPAdes_Assembly_Renamed_fa = "results/Assemblies/Unicycler_SPAdes_Assembly_{run}/{run}.SPAdes.Assembly.fasta",
        Ref_FA = refGenome_FA_PATH,
    output:
        NucDiff_SmallVariants_GFF = "results/AssemblyToRef_Variants/NucDiff_{run}/results/{run}.NucDiff_ref_snps.gff",
        NucDiff_SmallVariants_VCF = "results/AssemblyToRef_Variants/NucDiff_{run}/results/{run}.NucDiff_ref_snps.vcf",
        NucDiff_Struct_GFF = "results/AssemblyToRef_Variants/NucDiff_{run}/results/{run}.NucDiff_ref_struct.gff",
    conda:
        "Envs/nucdiff_2_0_3_Conda.yml"
    threads: 1
    shell:
        "nucdiff {input.Ref_FA} {input.Ill_SPAdes_Assembly_Renamed_fa} results/AssemblyToRef_Variants/NucDiff_{wildcards.run} {wildcards.run}.NucDiff --vcf yes \n"


rule SVanalyzer_SVrefine_Ill_SPAdes_Assembly_AlignTo_Ref:
    input:
        input_MUMmer_Delta_Aln_File = "results/AssemblyToRef_Variants/{run}.NucDiff.delta",
        Ill_SPAdes_Assembly_Renamed_fa = "results/Assemblies/Unicycler_SPAdes_Assembly_{run}/{run}.SPAdes.Assembly.fasta",
        Ref_FA = refGenome_FA_PATH,
    output:
        SVrefine_Ill_SPAdes_Assembly_VCF = "results/AssemblyToRef_Variants/SVanalyzer_SVrefine/{run}.SVrefine.vcf"
    threads: 1
    conda:
        "Envs/SV_Hack_V1.yml"
    shell:
        "SVrefine --delta {input.input_MUMmer_Delta_Aln_File} "
        "--ref_fasta {input.Ref_FA} "
        "--query_fasta {input.Ill_SPAdes_Assembly_Renamed_fa} "
        "--outvcf {output.SVrefine_Ill_SPAdes_Assembly_VCF} "
        "--refname 'NC_045512.2' --samplename {wildcards.run} "







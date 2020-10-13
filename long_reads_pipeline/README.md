# Command lines for Long Reads SV detection

> Deleted all generated intermediate files cause they are big

> I used the SRS6189915 for command line example, but it seems SRS6189915.fastq is too small to use.

## Reads Quality Control

### Get the stat plot

```
NanoPlot -t 2 --fastq raw_input/SRS6189915.fastq --plots hex dot -o 00.stats -p SRS6189915
```

> only stats of SRS6189915 after trimming was run with nanoplot

### Trim the reads

```
NanoFilt -q 10 --headcrop 50 raw_input/SRS6189915.fastq > 01.filter/SRS6189915.trimmed.fastq
```

## Mapping reads to reference genome

### minimap2 mapping

```
minimap2 -ax map-ont reference_genome/SARS_2CoV_2.fasta 01.filter/SRS6189915.trimmed.fastq > 02.mapping/SRS6189915.sam
```
### covert .sam to .bam and sort

```
samtools view 02.mapping/SRS6189915.sam | samtools sort > 02.mapping/SRS6189915.sorted.bam
```

### add MD tag to bam file (which need by Sniffles)

```
samtools calmd 02.mapping/SRS6189915.sorted.bam | samtools view -bS > 02.mapping/SRS6189915.final.bam
```

### indexing the result bam files

```
samtools index 02.mapping/SRS6189915.final.bam
```

## SV Calling

### Sniffles

```
sniffles -m 02.mapping/SRS6189915.final.bam -v 03.results/01.sniffles/SRS6189915.vcf
```

### SVIM

```
svim alignment 03.results/02.svim/SRS6189915 reference_genome/SARS_2CoV_2.fasta
```

### nanoSV

I skipped it because nanoSV only reports breakpoints, rather than SV types. So we can't merge with other vcfs effeciently

### cuteSV

```
mkdir -p tmp/SRS6189915; cuteSV 02.mapping/SRS6189915.final.bam 03.results/03.cuteSV/SRS6189915.vcf tmp/SRS6189915
```

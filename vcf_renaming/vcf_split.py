# This script is used for split vcf files by SV type and size
# Output vcf is named in the following format
# [RunAccession]_[SampleAccession]_[DUP/DEL/INV]_[small|large]_[ASM/LR/SR]_[RNAseq/Amplicon/WGS]_[VariantCaller].vcf

import os
import vcf
import csv

SV_SIZE_CUTOFF = 50

# File structure on DNAnexus is vcf_files/[ASM/LR/SR]/[VariantCaller]/[RunAccession].vcf
# This script is written to accommodate such file structure

vcf_files_dir = "vcf_files"
sample_run_mapping = "meta.csv"

# read metadata for RunAccession, SampleAccession, and Assay_Type
run2sample = dict()
run2assay = dict()
with open(sample_run_mapping, "r") as csvfile:
    csvreader = csv.reader(csvfile, delimiter='\t')
    next(csvreader)
    for row in csvreader:
        RunAccession = row[2]
        SampleAccession = row[5]
        Assay_Type = row[6]
        Platform = row[7]
        run2sample[RunAccession] = SampleAccession
        run2assay[RunAccession] = Assay_Type

read_types = ["short-reads", "long-reads", "asm"]

read_types_map = {"short-reads": "SR",
                  "long-reads": "LR",
                  "asm": "ASM"}

# fix naming error in the current dataset
renameing_dict = {"SRS6189924":"SRR11140744",
                  "SRS6189920":"SRR11140745",
                  "SRS6189919":"SRR11140746",
                  "SRS6189917":"SRR11140747",
                  "SRS6189918":"SRR11140748",
                  "SRS6189916":"SRR11140749",
                  "SRS6189914":"SRR11140750",
                  "SRS6189915":"SRR11140751"}

for read_type in read_types:
    read_type_dir = os.path.join(vcf_files_dir, read_type)
    if os.path.exists(read_type_dir):
        callers = os.listdir(read_type_dir)
        for caller in callers:
            for vcf_file in os.listdir(os.path.join(read_type_dir, caller)):
                
                RunAccession = vcf_file.split(".")[0]
                
                try:
                    SampleAccession = run2sample[RunAccession]
                except KeyError:
                    SampleAccession = RunAccession
                    RunAccession = renameing_dict[SampleAccession]
                
                Assay_Type = run2assay[RunAccession]
                
                vcf_reader = vcf.Reader(
                    filename=os.path.join(os.path.join(read_type_dir, caller, vcf_file)))
        
                del_small = []
                del_large = []
                inv_small = []
                inv_large = []
                dup_small = []
                dup_large = []
                
                for record in vcf_reader:
                    SvType = record.INFO["SVTYPE"]
                    if SvType in ["DEL","INV","DUP"]:
                        try:
                            SvLen = abs(int(record.INFO["SVLEN"][0]))
                        except TypeError:
                            SvLen = abs(record.INFO["SVLEN"])
                        except KeyError:
                            SvLen = abs(int(record.INFO["END"]) - int(record.POS))

                        # only include DEL/INV/DUP
                        if SvType == "DEL" and SvLen <= SV_SIZE_CUTOFF:
                            del_small.append(record)
                        if SvType == "DEL" and SvLen > SV_SIZE_CUTOFF:
                            del_large.append(record)
                        if SvType == "INV" and SvLen <= SV_SIZE_CUTOFF:
                            inv_small.append(record)
                        if SvType == "INV" and SvLen > SV_SIZE_CUTOFF:
                            inv_large.append(record)                           
                        if SvType == "DUP" and SvLen <= SV_SIZE_CUTOFF:
                            dup_small.append(record)
                        if SvType == "DUP" and SvLen > SV_SIZE_CUTOFF:
                            dup_large.append(record)                       
                
                if len(del_small) >= 0:
                    vcf_writer = vcf.Writer(open(f'vcf_split/{RunAccession}_{SampleAccession}_DEL_small_{read_types_map[read_type]}_{Assay_Type}_{caller}.vcf', 'w'),
                                            vcf_reader)
                    for record in del_small:
                        vcf_writer.write_record(record)
                    vcf_writer.close()
                        
                if len(del_large) >= 0:
                    vcf_writer = vcf.Writer(open(f'vcf_split/{RunAccession}_{SampleAccession}_DEL_large_{read_types_map[read_type]}_{Assay_Type}_{caller}.vcf', 'w'),
                                                  vcf_reader)
                    for record in del_large:
                        vcf_writer.write_record(record)
                    vcf_writer.close()
                        
                if len(inv_small) >= 0:
                    vcf_writer = vcf.Writer(open(f'vcf_split/{RunAccession}_{SampleAccession}_INV_small_{read_types_map[read_type]}_{Assay_Type}_{caller}.vcf', 'w'),
                                                  vcf_reader)
                    for record in inv_small:
                        vcf_writer.write_record(record)
                    vcf_writer.close()
                        
                if len(inv_large) >= 0:
                    vcf_writer = vcf.Writer(open(f'vcf_split/{RunAccession}_{SampleAccession}_INV_large_{read_types_map[read_type]}_{Assay_Type}_{caller}.vcf', 'w'),
                                                  vcf_reader)
                    for record in inv_large:
                        vcf_writer.write_record(record)
                    vcf_writer.close()
                        
                if len(dup_small) >= 0:
                    vcf_writer = vcf.Writer(open(f'vcf_split/{RunAccession}_{SampleAccession}_DUP_small_{read_types_map[read_type]}_{Assay_Type}_{caller}.vcf', 'w'),
                                                  vcf_reader)
                    for record in dup_small:
                        vcf_writer.write_record(record)
                    vcf_writer.close()
                        
                if len(dup_large) >= 0:
                    vcf_writer = vcf.Writer(open(f'vcf_split/{RunAccession}_{SampleAccession}_DUP_large_{read_types_map[read_type]}_{Assay_Type}_{caller}.vcf', 'w'),
                                                  vcf_reader)
                    for record in dup_large:
                        vcf_writer.write_record(record)
                    vcf_writer.close()


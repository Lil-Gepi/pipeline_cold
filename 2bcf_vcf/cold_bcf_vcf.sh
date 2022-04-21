#!/bin/sh

## 4. Filter mpileup to keep only SNPs
#
# input file, can also be "-" if you omit 3.c
in_bcf="/Users/ychen/COLD/result/bcf_piled/cold_piled.bcf"
# BED file of regions to exclude (must match the reference)
# norepeats_bed="/path/to/norepeats.bed"
# filter expression removin non-SNPs and positions of dubious
# coverage depth (adjust as needed; here, ~2200x coverage was
# expected across all samples in the mpileup)
flt_expr='TYPE = "snp" & INFO/DP > 1500 & INFO/DP < 7000'
# output path (bgzipped vcf is the most widely used format)
out_vcf="/Users/ychen/COLD/result/vcf/cold_piled_flt.vcf.gz"

# 4.a filter SNPs close to INDELs or in repeat regions
# -g 5       ... filter out SNPs that are this close to INDELs
# -T bedfile ... only output positions described by the provided BED file
# bcftools filter -g 5 -T "$norepeats_bed" -Ou "$in_bcf" |\
bcftools filter -g 5 -Ou "$in_bcf" |\
  # 4.b keep only SNPs and adjust coverage depth thresholds as appropriate
  # (here ~2200x was expected across all samples)
  bcftools filter -i "$flt_expr" -Ou - | \
  # 4.c make positions single-allelic for easier manipulation
  bcftools norm -m- -Oz - > "$out_vcf"

exit

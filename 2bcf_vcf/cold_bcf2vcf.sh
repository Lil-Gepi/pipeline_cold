
## 3. Combine multiple single-sample pilups into a multi-sample pileup
#
# folder where only the annotated pileups (output of step 2) to be combined are found
# "ls $in_path/*.bcf" shows the files in the same order in which you get the columns in the mpileup
in_folder="/Users/ychen/COLD/result/bcf_single"
out_bcf="/Users/ychen/COLD/result/bcf_piled/cold_piled.bcf"
# 3.a merge multiple raw pileups into multi-sample mpileup
bcftools merge --no-index --merge both --threads 8 -Ov -o - "$in_folder"/*.bcf |\
  # 3.b annotate with FORMAT/AF (observed allele frequencies) and FORMAT/XF (expected allele frequencies
  # under a multinomial sampling model) tags; also add  FORMAT/SAC (sum of allele counts) tag for convenience
  # finally, change REF of positions with a reference count of 0 to ALT{1} and modify all affected tags accordingly
  # tag such positions with the INFO/RMOD flag
  ../rupert_pipe_v3/post-merging.awk |\
  # directly pipe into the step 4.a or
  # 3.c save VCF to compressed format
  bcftools view --no-version -Ob > "$out_bcf"

## 4. Filter mpileup to keep only SNPs
#
# input file, can also be "-" if you omit 3.c
in_bcf="/Users/ychen/COLD/result/bcf_piled/cold_piled.bcf"
# BED file of regions to exclude (must match the reference)
# norepeats_bed="/path/to/norepeats.bed"
# filter expression removin non-SNPs and positions of dubious
# coverage depth (adjust as needed; here, ~2200x coverage was
# expected across all samples in the mpileup)
flt_expr='TYPE = "snp" & INFO/DP > 1000 & INFO/DP < 4500'
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

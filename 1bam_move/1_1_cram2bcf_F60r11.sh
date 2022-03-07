## 1. Prepare the raw single-sample pileup

# adjust these variables as needed
# read length
read_length=150
# output path (I use BCF, the compressed binary version of VCF, here)
out_bcf="/Users/ychen/COLD/data/coldF60_r11/coldF60_r11_raw.bcf"
# reference genome in FASTA format
ref_fa="/Users/ychen/COLD/reference/dsimM252v1.2+microbiome.fa"
# max. depth to consider, set to 5x expected autosomal coverage
max_cov=400
# input path (if your samples have been sequenced on multiple lanes, the CRAMs
# for a single sample need to be samtools-merged into one file instead, see below)
in_folder_329="/Users/ychen/COLD/data/coldF60r11_1"
in_folder_356="/Users/ychen/COLD/data/coldF60r11_2"

# 1.a convert to SAM (needed by awk scripts)
# samtools view -h "$in_cram" |\
# or with multiple CRAMs for the same sample
samtools merge -O SAM -o - "$in_folder_329/coldF60r11_1_a.bam" "$in_folder_329/coldF60r11_1_b.bam" "$in_folder_356/coldF60r11_2_a.bam" |\
  # 1.b flag short fragments  (replace $rl with the actual read length)
  /Users/ychen/COLD/pipeline_cold/1bam_move/rupert_pipe_v2/flag-short.awk -v READ_LENGTH=$read_length |\
  # 1.c recalculate mapping quality
  /Users/ychen/COLD/pipeline_cold/1bam_move/rupert_pipe_v2/remapq.awk -v NORD=7.5 |\
  # 1.d make single-sample pileup
  # --incl-flags 0x2 ... include only proper pairs
  # -q 10 ... include only reads with a mapping quality that apsses the threshold
  #          (may need to be adjusted if remapq step is omitted or NORD!=7.5 is used)
  #          (if in doubt, check a histogram of the MAPQ values in $in_cram)
  # -Q 20 ... count only bases with min(base quality,BAQ) passing the threshold
  # -D    ... enable BAQ calculation for all positions (check bcftools manpage for details)
  # -a DP,AD ... annotate with additional fields (AD = allele depth for each allele is needed below!)
  bcftools mpileup -f "$ref_fa" -d $max_cov --incl-flags 0x2 -q 10 -Q 20 -D -a DP,AD,QS,SCR -Ob - > "$out_bcf"

## 2. Annotate raw single-sample pileup with more FORMAT tags (to save the information from INFO)
# input path (a pileup generated in step 1 in BCF, VCF, or bgzipped VCF format, cannot be a pipe!)
in_bcf="/Users/ychen/COLD/data/coldF60_r11/coldF60_r11_raw.bcf"
# output path (use - to write uncompressed BCF to stdout, or leave empty to write VCF to stdout)
out_bcf="/Users/ychen/COLD/data/coldF60_r11/coldF60_r11_pre.bcf"
./rupert_pipe_v2/pre-merging.sh "$in_bcf" "$out_bcf"

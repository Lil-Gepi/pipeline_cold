## 1. Prepare the raw single-sample pileup

# adjust these variables as needed

cold_folder="/Users/ychen/COLD"
# read length
read_length_329=150
read_length_356=125
# output path (I use BCF, the compressed binary version of VCF, here)
out_bcf="${cold_folder}/data/coldF60_r11/coldF60_r11_piped_raw.bcf"
# reference genome in FASTA format
ref_fa="${cold_folder}/reference/dsimM252v1.2+microbiome.fa"
# max. depth to consider, set to 5x expected autosomal coverage
max_cov=400
# input path (if your samples have been sequenced on multiple lanes, the CRAMs
# for a single sample need to be samtools-merged into one file instead, see below)
in_folder_329="/Volumes/Data/329"
in_folder_356="/Volumes/Data/356"

# 1.a convert to SAM (needed by awk scripts)
# samtools view -h "$in_cram" |\
# or with multiple CRAMs for the same sample
samtools merge -O SAM -o - \
<(samtools view -T "$ref_fa" -h  "$in_folder_329/a/LB_329.D707.cram" | ./rupert_pipe_v2/flag-short.awk -v READ_LENGTH=$read_length_329 ) \
<(samtools view -T "$ref_fa" -h  "$in_folder_329/b/LB_329.D707.cram" | ./rupert_pipe_v2/flag-short.awk -v READ_LENGTH=$read_length_329 ) \
<(samtools view -T "$ref_fa" -h  "$in_folder_356/a/LB_356.D707+D507.cram" | ./rupert_pipe_v2/flag-short.awk -v READ_LENGTH=$read_length_356 )|\
  # 1.b flag short fragments  (replace $rl with the actual read length)
  # ./rupert_pipe_v2/flag-short.awk -v READ_LENGTH=$read_length |\
  # 1.c recalculate mapping quality
  ./rupert_pipe_v2/remapq.awk -v NORD=7.5 |\
  # 1.d make single-sample pileup
  # --incl-flags 0x2 ... include only proper pairs
  # -q 10 ... include only reads with a mapping quality that apsses the threshold
  #          (may need to be adjusted if remapq step is omitted or NORD!=7.5 is used)
  #          (if in doubt, check a histogram of the MAPQ values in $in_cram)
  # -Q 20 ... count only bases with min(base quality,BAQ) passing the threshold
  # -D    ... enable BAQ calculation for all positions (check bcftools manpage for details)
  # -a DP,AD ... annotate with additional fields (AD = allele depth for each allele is needed below!)
  bcftools mpileup -f "$ref_fa" -d $max_cov --skip-all-unset 0x2 -q 10 -Q 20 -D -a DP,AD,QS,SCR -Ob - > "$out_bcf"

## 2. Annotate raw single-sample pileup with more FORMAT tags (to save the information from INFO)
# input path (a pileup generated in step 1 in BCF, VCF, or bgzipped VCF format, cannot be a pipe!)
in_bcf="${cold_folder}/data/coldF60_r11/coldF60_r11_piped_raw.bcf"
# output path (use - to write uncompressed BCF to stdout, or leave empty to write VCF to stdout)
out_bcf="${cold_folder}/data/coldF60_r11/coldF60_r11_piped_pre.bcf"
./rupert_pipe_v2/pre-merging.sh "$in_bcf" "$out_bcf"

exit

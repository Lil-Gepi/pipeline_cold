## 1. Prepare the raw single-sample pileup

# adjust these variables as needed

cold_dir="/Users/ychen/COLD"
# read length
rl=120
#rl_356=125
# output path (I use BCF, the compressed binary version of VCF, here)
out_bcf="${cold_dir}/data/F40_r18/coldF40_r18_raw.bcf"
# reference genome in FASTA format
ref_fa="${cold_dir}/reference/dsimM252v1.2+microbiome.fa"
# max. depth to consider, set to 5x expected autosomal coverage
max_cov=400
# input path (if your samples have been sequenced on multiple lanes, the CRAMs
# for a single sample need to be samtools-merged into one file instead, see below)
in_dir_1="/Volumes/Data/185"


# 1.a convert to SAM (needed by awk scripts)
# samtools view -h "$in_cram" |\
# or with multiple CRAMs for the same sample
samtools merge -O SAM -o - \
<(samtools view -T "$ref_fa" -h -F 0x400 "$in_dir_1/a/LB_185.D712+D505.cram" | ./rupert_pipe_v2/flag-short.awk -v READ_LENGTH=$rl) \
<(samtools view -T "$ref_fa" -h -F 0x400 "$in_dir_1/b/LB_185.D712+D505.cram" | ./rupert_pipe_v2/flag-short.awk -v READ_LENGTH=$rl) \
<(samtools view -T "$ref_fa" -h -F 0x400 "$in_dir_1/c/LB_185.D712+D505.cram" | ./rupert_pipe_v2/flag-short.awk -v READ_LENGTH=$rl) \
<(samtools view -T "$ref_fa" -h -F 0x400 "$in_dir_1/d/LB_185.D712+D505.cram" | ./rupert_pipe_v2/flag-short.awk -v READ_LENGTH=$rl) |\
# testing the time it takes for auto detecting read length, so no read length was provided to flag-short.awk
samtools collate -Ouf - | samtools fixmate -ru - - | samtools view -f 0x2 -h - |\
  # collate sort the SAM files by order of qname, fixmate pairs only primary pairs, the following view command keeps only primary pairs
  # Non-primary pairs flag 0x100, is an additional line in the CRAM for a read that already has a primary alignment.
  # It is independent of the "proper pair" flag 0x2, so a "proper pair" can be a primary alignment of one read and a secondary of its mate.
  # Sometimes, collate -f leaves broken pairs, i.e. reads that have a 0x2 flag, but no mate, so we need the fixmate step to rewrite those flags,
  # and the view -f 0x2 to get only the real remaining pairs
  # So with fixmate and view -f 0x2 guarantees that you only have proper primary pairs for remapq

  # 1.b flag short fragments  (replace $rl with the actual read length)
  # ./rupert_pipe_v2/flag-short.awk -v READ_LENGTH=$read_length |\
  # 1.c recalculate mapping quality
  # remapq wants name-sorted, that's why we are collating the SAM based on qname
  ./rupert_pipe_v2/remapq.awk -v NORD=7.5 |\
  # and before entering bcftools mpileup, the input should be sorted back to coordinate order
  samtools sort -@3 -m 4G -u -T /tmp/sortsam.$(basename ${out_bcf%.bcf}).$tag - |\
  # 1.d make single-sample pileup
  # --incl-flags 0x2 ... include only proper pairs
  # -q 10 ... include only reads with a mapping quality that apsses the threshold
  #          (may need to be adjusted if remapq step is omitted or NORD!=7.5 is used)
  #          (if in doubt, check a histogram of the MAPQ values in $in_cram)
  # -Q 20 ... count only bases with min(base quality,BAQ) passing the threshold
  # -D    ... enable BAQ calculation for all positions (check bcftools manpage for details)
  # -a DP,AD ... annotate with additional fields (AD = allele depth for each allele is needed below!)
  bcftools mpileup -f "$ref_fa" -d $max_cov --skip-all-unset 0x2 -q 10 -Q 20 -D -a DP,AD,QS,SCR -Ob --threads 3 - > "$out_bcf"

## 2. Annotate raw single-sample pileup with more FORMAT tags (to save the information from INFO)
# input path (a pileup generated in step 1 in BCF, VCF, or bgzipped VCF format, cannot be a pipe!)
in_bcf="${cold_dir}/data/F40_r18/coldF40_r18_raw.bcf"
# output path (use - to write uncompressed BCF to stdout, or leave empty to write VCF to stdout)
out_bcf="${cold_dir}/data/F40_r18/coldF40_r18_pre.bcf"
#./rupert_pipe_v2/pre-merging.sh "$in_bcf" "$out_bcf"
# the pre-merging.sh has the "mktemp" syntax difference between linux and mac,
# so the temprary files ann.hdr and ann.txt.gz are not created, instead use the following command to do the same job.
bcftools annotate -Ov -x INFO/IMF,INFO/VDB,INFO/RPBZ,INFO/MQBZ,INFO/MQSBZ,INFO/SCBZ,INFO/BQBZ,INFO/FS,INFO/SGB,INFO/I16,INFO/QS,INFO/MQ0F,FORMAT/PL "$in_bcf" |\
 ./rupert_pipe_v2/info2fmt.awk -v tags=DP > "$out_bcf"
# these long options just remove a bunch of tags from the INFO field to make the files much smaller

exit

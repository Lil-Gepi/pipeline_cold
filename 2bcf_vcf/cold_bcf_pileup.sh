#!/bin/sh

## 3. Combine multiple single-sample pilups into a multi-sample pileup
#
# folder where only the annotated pileups (output of step 2) to be combined are found
# "ls $in_path/*.bcf" shows the files in the same order in which you get the columns in the mpileup
in_folder="/Users/ychen/COLD/result/bcf_single"
out_bcf="/Users/ychen/COLD/result/bcf_piled/cold_piled.bcf"

#!/bin/bash
# NOTE : Quote it else use array to avoid problems #

for bcf_files in /Users/ychen/COLD/result/bcf_single/*.bcf
do
  bcftools index -f --threads 14 $bcf_files
done

# 3.a merge multiple raw pileups into multi-sample mpileup
bcftools merge --no-index --merge both --threads 14 -Ov -o - "$in_folder"/*.bcf |\
  # 3.b annotate with FORMAT/AF (observed allele frequencies) and FORMAT/XF (expected allele frequencies
  # under a multinomial sampling model) tags; also add  FORMAT/SAC (sum of allele counts) tag for convenience
  # finally, change REF of positions with a reference count of 0 to ALT{1} and modify all affected tags accordingly
  # tag such positions with the INFO/RMOD flag
  /Users/ychen/COLD/pipeline_cold/rupert_pipe_v3/post-merging.awk |\
  # directly pipe into the step 4.a or
  # 3.c save VCF to compressed format
  bcftools view --thread 14 --no-version -Ob -o "$out_bcf"

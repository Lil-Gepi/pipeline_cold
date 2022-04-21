#!/bin/bash

in_dir="/Users/ychen/COLD/result/bcf_single_raw"
out_dir="/Users/ychen/COLD/result/bcf_single"

for bcf_pre_file in $(find $in_dir -type f -maxdepth 1 -name "*_raw.bcf" -exec basename {} \;)
do
  sample_id=$(echo $bcf_pre_file | cut -f1,2 -d '_')
  bcftools index -f --threads 16 "${in_dir}/${bcf_pre_file}"
  bcftools view --threads 16 --regions X,2L,2R,3L,3R,4 --no-version -Ob "${in_dir}/${bcf_pre_file}"|
  bcftools annotate -Ov -x INFO/IMF,INFO/VDB,INFO/RPBZ,INFO/MQBZ,INFO/MQSBZ,INFO/SCBZ,INFO/BQBZ,INFO/FS,INFO/SGB,INFO/I16,INFO/QS,INFO/MQ0F,FORMAT/PL |\
  /Users/ychen/COLD/pipeline_cold/rupert_pipe_v3/info2fmt.awk -v tags=DP | bcftools view --no-version -Ob -o "${out_dir}/${sample_id}_pre.bcf"

done

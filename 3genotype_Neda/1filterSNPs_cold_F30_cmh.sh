#!/bin/sh
# to filter the SNPs piled up in our mapped data based on Neda's SNP set.
# this awk command line is from Sheng-kai
for repl in 11 12 13 14
do
  awk 'BEGIN{OFS="\t"} NR==FNR {f1[$1$2] = $0; next} ($1$2 in f1) {print $0}' \
  ~/RS/reference/Dsim_base_F60_bam.sorted.rmdup.filter.mpileup.Q20.polymorphic_indel_repeat_Ytranslocremoved_masked_MajorChr.sync_furthermasked.cmh \
  ~/RS/result/cold_F30_cmh/cold_F30_r${repl}_all_chr_mq20_bq20.sync > \
  ~/RS/result/cold_F30_cmh/cold_F30_r${repl}_all_chr_mq20_bq20_Neda_filtered.sync
done

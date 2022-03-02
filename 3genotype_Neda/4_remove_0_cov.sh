#!/bin/bash

grep -v "[[:blank:]]0:0:0:0" cold_F30_r11_all_chr_mq20_bq20_Neda_filtered_base_duplicated.sync \
> cold_F30_r11_all_chr_mq20_bq20_Neda_filtered_base_duplicated_nozerocov.sync

grep -v "[[:blank:]]0:0:0:0" cold_F30_r12_all_chr_mq20_bq20_Neda_filtered_base_duplicated.sync \
> cold_F30_r12_all_chr_mq20_bq20_Neda_filtered_base_duplicated_nozerocov.sync

grep -v "[[:blank:]]0:0:0:0" cold_F30_r13_all_chr_mq20_bq20_Neda_filtered_base_duplicated.sync \
> cold_F30_r13_all_chr_mq20_bq20_Neda_filtered_base_duplicated_nozerocov.sync

grep -v "[[:blank:]]0:0:0:0" cold_F30_r14_all_chr_mq20_bq20_Neda_filtered_base_duplicated.sync \
> cold_F30_r14_all_chr_mq20_bq20_Neda_filtered_base_duplicated_nozerocov.sync

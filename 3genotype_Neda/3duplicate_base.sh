#!/bin/bash

awk '{ $4=$4 OFS $4 OFS $4 OFS $4 OFS $4}1' OFS='\t' cold_F30_r11_all_chr_mq20_bq20_Neda_filtered.sync > cold_F30_r11_all_chr_mq20_bq20_Neda_filtered_base_duplicated.sync

awk '{ $4=$4 OFS $4 OFS $4 OFS $4 OFS $4}1' OFS='\t' cold_F30_r12_all_chr_mq20_bq20_Neda_filtered.sync > cold_F30_r12_all_chr_mq20_bq20_Neda_filtered_base_duplicated.sync

awk '{ $4=$4 OFS $4 OFS $4 OFS $4 OFS $4}1' OFS='\t' cold_F30_r13_all_chr_mq20_bq20_Neda_filtered.sync > cold_F30_r13_all_chr_mq20_bq20_Neda_filtered_base_duplicated.sync

awk '{ $4=$4 OFS $4 OFS $4 OFS $4 OFS $4}1' OFS='\t' cold_F30_r14_all_chr_mq20_bq20_Neda_filtered.sync > cold_F30_r14_all_chr_mq20_bq20_Neda_filtered_base_duplicated.sync

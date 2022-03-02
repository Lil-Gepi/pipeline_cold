#!/bin/sh

cat \
~/RS/result/grand_PCA_all_generations_2L.sync ~/RS/result/grand_PCA_all_generations_2R.sync \
~/RS/result/grand_PCA_all_generations_3L.sync ~/RS/result/grand_PCA_all_generations_3R.sync \
~/RS/result/grand_PCA_all_generations_4.sync ~/RS/result/grand_PCA_all_generations_X.sync > \
~/RS/result/grand_PCA_all_generations_all_chr_mq20_bq20.sync

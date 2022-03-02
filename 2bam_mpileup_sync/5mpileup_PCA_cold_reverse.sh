#!/bin/sh

for chr in 2L 2R 3L 3R 4 X
do
  samtools mpileup -B -q 20 -Q 0 -f ~/RS/reference/dsimM252v1.2+microbiome.fa \
  ~/RS/data/coldF40r11/coldF40r11_${chr}.bam ~/RS/data/coldF40r12/coldF40r12_${chr}.bam \
  ~/RS/data/coldF40r13/coldF40r13_${chr}.bam ~/RS/data/coldF40r14/coldF40r14_${chr}.bam \
  ~/RS/data/F30r11sub1/F30r11sub1_${chr}.bam ~/RS/data/F30r11sub2/F30r11sub2_${chr}.bam \
  ~/RS/data/F30r11sub3/F30r11sub3_${chr}.bam ~/RS/data/F30r11sub4/F30r11sub4_${chr}.bam \
  ~/RS/data/F30r11sub5/F30r11sub5_${chr}.bam \
  ~/RS/data/F30r12sub1/F30r12sub1_${chr}.bam ~/RS/data/F30r12sub2/F30r12sub2_${chr}.bam \
  ~/RS/data/F30r12sub3/F30r12sub3_${chr}.bam ~/RS/data/F30r12sub4/F30r12sub4_${chr}.bam \
  ~/RS/data/F30r12sub5/F30r12sub5_${chr}.bam \
  ~/RS/data/F30r13sub1/F30r13sub1_${chr}.bam ~/RS/data/F30r13sub2/F30r13sub2_${chr}.bam \
  ~/RS/data/F30r13sub3/F30r13sub3_${chr}.bam ~/RS/data/F30r13sub4/F30r13sub4_${chr}.bam \
  ~/RS/data/F30r13sub5/F30r13sub5_${chr}.bam \
  ~/RS/data/F30r14sub1/F30r14sub1_${chr}.bam ~/RS/data/F30r14sub2/F30r14sub2_${chr}.bam \
  ~/RS/data/F30r14sub3/F30r14sub3_${chr}.bam ~/RS/data/F30r14sub4/F30r14sub4_${chr}.bam \
  ~/RS/data/F30r14sub5/F30r14sub5_${chr}.bam | \
  # need to figure out how to specifz the Pop and Gen instead of * all bam files.
  # maybe partition the gen? and R? so that I can put them into the code.
  java -jar ~/RS/pipeline/popoolation2/mpileup2sync.jar --input /dev/stdin \
  --output ~/RS/result/cold_F30_PCA/cold_F30_${chr}.sync --fastq-type sanger \
  --min-qual 20 --threads 16
done


cat \
~/RS/result/cold_F30_PCA/cold_F30_2L.sync ~/RS/result/cold_F30_PCA/cold_F30_2R.sync \
~/RS/result/cold_F30_PCA/cold_F30_3L.sync ~/RS/result/cold_F30_PCA/cold_F30_3R.sync \
~/RS/result/cold_F30_PCA/cold_F30_4.sync ~/RS/result/cold_F30_PCA/cold_F30_X.sync > \
~/RS/result/cold_F30_PCA/cold_F30_all_chr_mq20_bq20.sync

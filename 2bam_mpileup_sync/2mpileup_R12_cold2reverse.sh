#!/bin/sh

# do the job from bam to cram to sync for each subreplicate
# bamfile are like F30r11sub1 etc.
for bamfile in coldF40r12 F30r12sub1 F30r12sub2 F30r12sub3 F30r12sub4 F30r12sub5
do
  samtools index -b ~/RS/data/${bamfile}/${bamfile}.bam
  for chr in 2L 2R 3L 3R 4 X
  do
    samtools view -b ~/RS/data/${bamfile}/${bamfile}.bam ${chr} > \
    ~/RS/data/${bamfile}/${bamfile}_${chr}.bam
    samtools index ~/RS/data/${bamfile}/${bamfile}_${chr}.bam
  done
done

for chr in 2L 2R 3L 3R 4 X
do
  samtools mpileup -B -q 20 -Q 0 -f ~/RS/reference/dsimM252v1.2+microbiome.fa \
  ~/RS/data/coldF40r12/coldF40r12_${chr}.bam ~/RS/data/F30r12sub1/F30r12sub1_${chr}.bam \
  ~/RS/data/F30r12sub2/F30r12sub2_${chr}.bam ~/RS/data/F30r12sub3/F30r12sub3_${chr}.bam \
  ~/RS/data/F30r12sub4/F30r12sub4_${chr}.bam ~/RS/data/F30r12sub5/F30r12sub5_${chr}.bam | \
  # need to figure out how to specifz the Pop and Gen instead of * all bam files.
  # maybe partition the gen? and R? so that I can put them into the code.
  java -jar ~/RS/pipeline/popoolation2/mpileup2sync.jar --input /dev/stdin \
  --output ~/RS/result/cold_F30_cmh/cold_F30_r12_${chr}.sync --fastq-type sanger \
  --min-qual 20 --threads 12

done


cat \
~/RS/result/cold_F30_cmh/cold_F30_r12_2L.sync ~/RS/result/cold_F30_cmh/cold_F30_r12_2R.sync \
~/RS/result/cold_F30_cmh/cold_F30_r12_3L.sync ~/RS/result/cold_F30_cmh/cold_F30_r12_3R.sync \
~/RS/result/cold_F30_cmh/cold_F30_r12_4.sync ~/RS/result/cold_F30_cmh/cold_F30_r12_X.sync > \
~/RS/result/cold_F30_cmh/cold_F30_r12_all_chr_mq20_bq20.sync

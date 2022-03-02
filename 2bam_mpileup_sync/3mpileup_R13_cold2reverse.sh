#!/bin/sh

# do the job from bam to cram to sync for each subreplicate
# bamfile are like F30r13sub1 etc.
for bamfile in coldF40r13 F30r13sub1 F30r13sub2 F30r13sub3 F30r13sub4 F30r13sub5
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
  ~/RS/data/coldF40r13/coldF40r13_${chr}.bam ~/RS/data/F30r13sub1/F30r13sub1_${chr}.bam \
  ~/RS/data/F30r13sub2/F30r13sub2_${chr}.bam ~/RS/data/F30r13sub3/F30r13sub3_${chr}.bam \
  ~/RS/data/F30r13sub4/F30r13sub4_${chr}.bam ~/RS/data/F30r13sub5/F30r13sub5_${chr}.bam | \
  # need to figure out how to specifz the Pop and Gen instead of * all bam files.
  # maybe partition the gen? and R? so that I can put them into the code.
  java -jar ~/RS/pipeline/popoolation2/mpileup2sync.jar --input /dev/stdin \
  --output ~/RS/result/cold_F30_cmh/cold_F30_r13_${chr}.sync --fastq-type sanger \
  --min-qual 20 --threads 12

done


cat \
~/RS/result/cold_F30_cmh/cold_F30_r13_2L.sync ~/RS/result/cold_F30_cmh/cold_F30_r13_2R.sync \
~/RS/result/cold_F30_cmh/cold_F30_r13_3L.sync ~/RS/result/cold_F30_cmh/cold_F30_r13_3R.sync \
~/RS/result/cold_F30_cmh/cold_F30_r13_4.sync ~/RS/result/cold_F30_cmh/cold_F30_r13_X.sync > \
~/RS/result/cold_F30_cmh/cold_F30_r13_all_chr_mq20_bq20.sync

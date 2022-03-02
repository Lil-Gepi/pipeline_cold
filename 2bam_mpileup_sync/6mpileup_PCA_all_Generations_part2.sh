#!/bin/sh

for chr in 2L 2R 3L 3R 4 X
do
  samtools mpileup -B -q 20 -Q 0 -f ~/RS/reference/dsimM252v1.2+microbiome.fa \
  $(ls -R ~/RS/data/*/*_${chr}.bam) | \
  java -jar ~/RS/pipeline/popoolation2/mpileup2sync.jar --input /dev/stdin \
  --output ~/RS/result/grand_PCA_all_generations_${chr}.sync --fastq-type sanger \
  --min-qual 20 --threads 12
done

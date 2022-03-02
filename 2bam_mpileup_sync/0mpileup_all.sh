#!/bin/sh

# do the job from bam to cram to sync for each subreplicate
# bamfile are like F30r11sub1 etc.
for bamfile in $(ls ~/RS/data/)
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
  $(ls -R /Users/ychen/RS/data/*/*_${chr}.bam) | \
  # need to figure out how to specifz the Pop and Gen instead of * all bam files.
  # maybe partition the gen? and R? so that I can put them into the code.
  java -jar ~/RS/pipeline/popoolation2/mpileup2sync.jar --input /dev/stdin \
  --output ~/RS/result/reverse_${chr}.sync --fastq-type sanger \
  --min-qual 20 --threads 12

done


cat \
~/RS/result/F30_2L.sync ~/RS/result/F30_2R.sync \
~/RS/result/F30_3L.sync ~/RS/result/F30_3R.sync \
~/RS/result/F30_4.sync ~/RS/result/F30_X.sync > \
~/RS/result/F30_all_chr_mq20_bq20.sync




# old codes :
# samtools view -b /Users/ychen/RS/flight_test/data/F30/LB_311.TI02.bam 2L > ../data/F30/Pool_311.F30r11sub1_2L.bam &&\
# samtools view -b /Users/ychen/RS/flight_test/data/F30/LB_311.TI02.bam 2R > ../data/F30/Pool_311.F30r11sub1_2R.bam &&\
# samtools view -b /Users/ychen/RS/flight_test/data/F30/LB_311.TI02.bam 3L > ../data/F30/Pool_311.F30r11sub1_3L.bam &&\
# samtools view -b /Users/ychen/RS/flight_test/data/F30/LB_311.TI02.bam 3R > ../data/F30/Pool_311.F30r11sub1_3R.bam &&\
# samtools view -b /Users/ychen/RS/flight_test/data/F30/LB_311.TI02.bam 4 > ../data/F30/Pool_311.F30r11sub1_4.bam &&\
# samtools view -b /Users/ychen/RS/flight_test/data/F30/LB_311.TI02.bam X > ../data/F30/Pool_311.F30r11sub1_X.bam &&\
#
# samtools index ../data/F30/Pool_311.F30r11sub1_2L.bam &&\
# samtools index ../data/F30/Pool_311.F30r11sub1_2R.bam &&\
# samtools index ../data/F30/Pool_311.F30r11sub1_3L.bam &&\
# samtools index ../data/F30/Pool_311.F30r11sub1_3R.bam &&\
# samtools index ../data/F30/Pool_311.F30r11sub1_4.bam &&\
# samtools index ../data/F30/Pool_311.F30r11sub1_X.bam &&\
#
# samtools mpileup -B -q 20 -Q 0 -f ~/RS/flight_test/data/reference/dsimM252v1.2+microbiome.fa ../data/F30/Pool_311.F30r11sub1_2L.bam > ../data/F30/Pool_311.F30r11sub1_2L.pileup &&\
# samtools mpileup -B -q 20 -Q 0 -f ~/RS/flight_test/data/reference/dsimM252v1.2+microbiome.fa ../data/F30/Pool_311.F30r11sub1_2R.bam > ../data/F30/Pool_311.F30r11sub1_2R.pileup &&\
# samtools mpileup -B -q 20 -Q 0 -f ~/RS/flight_test/data/reference/dsimM252v1.2+microbiome.fa ../data/F30/Pool_311.F30r11sub1_3L.bam > ../data/F30/Pool_311.F30r11sub1_3L.pileup &&\
# samtools mpileup -B -q 20 -Q 0 -f ~/RS/flight_test/data/reference/dsimM252v1.2+microbiome.fa ../data/F30/Pool_311.F30r11sub1_3R.bam > ../data/F30/Pool_311.F30r11sub1_3R.pileup &&\
# samtools mpileup -B -q 20 -Q 0 -f ~/RS/flight_test/data/reference/dsimM252v1.2+microbiome.fa ../data/F30/Pool_311.F30r11sub1_4.bam > ../data/F30/Pool_311.F30r11sub1_4.pileup &&\
# samtools mpileup -B -q 20 -Q 0 -f ~/RS/flight_test/data/reference/dsimM252v1.2+microbiome.fa ../data/F30/Pool_311.F30r11sub1_X.bam > ../data/F30/Pool_311.F30r11sub1_X.pileup &&\
#
# # concate all the files into one
# cd ../data/F30/
# cat Pool_311.F30r11sub1_2L.pileup Pool_311.F30r11sub1_2R.pileup Pool_311.F30r11sub1_3L.pileup Pool_311.F30r11sub1_3R.pileup Pool_311.F30r11sub1_4.pileup Pool_311.F30r11sub1_X.pileup > Pool_311.F30r11sub1_allChr.pileup

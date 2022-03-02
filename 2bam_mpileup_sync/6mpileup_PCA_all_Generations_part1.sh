#!/bin/sh

# do the job from bam to cram to sync for each subreplicate
# bamfile are like F30r11sub1 etc.
bamfile_id_list=(F10r12sub1 F10r12sub2 F10r13sub1 F10r13sub2 F10r14sub1 F10r14sub2 F20r11sub1 F20r11sub2 F20r12sub1 F20r12sub2 F20r13sub1 F20r13sub2 F20r14sub1 F20r14sub2 F40r11sub1 F40r11sub2 F40r12sub1 F40r12sub2 F40r13sub1 F40r13sub2 F40r14sub1 F40r14sub2 F50r11sub1 F50r11sub2 F50r12sub1 F50r12sub2 F50r13sub1 F50r13sub2 F50r14sub1 F50r14sub2 F60r11sub1 F60r11sub2 F60r12sub1 F60r12sub2 F80r11sub1 F80r11sub2 F80r12sub1 F80r12sub2 F80r13sub1 F80r13sub2 F80r14sub1 F80r14sub2)
for bamfile in "${bamfile_id_list[@]}"
do
  samtools index -b ~/RS/data/${bamfile}/${bamfile}.bam
  for chr in 2L 2R 3L 3R 4 X
  do
    echo "samtools view -b ~/RS/data/${bamfile}/${bamfile}.bam ${chr} > \
    ~/RS/data/${bamfile}/${bamfile}_${chr}.bam; samtools index ~/RS/data/${bamfile}/${bamfile}_${chr}.bam"
  done | parallel -j 8
done

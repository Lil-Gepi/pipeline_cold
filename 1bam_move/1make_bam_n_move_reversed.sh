#!/bin/bash
samtools faidx ~/RS/reference/dsimM252v1.2+microbiome.fa
# It could be better if we intake a txt containing .cram file name and coordinate biological meaning, so that he bam file name make sense

for pool in 267 311 312 488 564 565 566
do
  for cramfile in $(find /Volumes/Data/${pool}/a/ -type f -maxdepth 1 -name "*.cram" -exec basename {} \;)
  do
    cramfilename=${cramfile%.*}
    bamfilename=$(awk -v var="$cramfilename" '$1==var{print $2}' ~/RS/pipeline/1bam_move/cram2bam_name.txt)
    if [ ! -d ~/RS/data/${bamfilename}/ ]; then
      mkdir ~/RS/data/${bamfilename}/
    fi
  done
done

for pool in 267 311 312 488 564 565 566
do
  for cramfile in $(find /Volumes/Data/${pool}/a/ -type f -maxdepth 1 -name "*.cram" -exec basename {} \;)
  do
    cramfilename=${cramfile%.*}
    bamfilename=$(awk -v var="$cramfilename" '$1==var{print $2}' ~/RS/pipeline/1bam_move/cram2bam_name.txt)
    echo "samtools view -b -q 20 -F 0x400 -T ~/RS/reference/dsimM252v1.2+microbiome.fa -o ~/RS/data/${bamfilename}/${bamfilename}.bam /Volumes/Data/${pool}/a/${cramfile}"
  done
done | parallel -j 10

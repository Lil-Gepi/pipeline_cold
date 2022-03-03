#!/bin/bash
for pool in 126_2 132 133 134 185 313 329 356 357
do
  for cramfile in $(find /Volumes/Data/${pool}/a/ -type f -maxdepth 1 -name "*.cram" -exec basename {} \;)
  do
    cramfilename=${cramfile%.*}
    bamfilename=$(awk -v var="$cramfilename" '$1==var{print $2}' ~/COLD/pipeline/1bam_move/cram2bam_name.txt)
    if [ ! -d ~/COLD/data/${bamfilename}/ ]; then
      mkdir ~/COLD/data/${bamfilename}/
    fi
  done
done

for pool in 126_2 132 133 134 185 313 329 356 357
do
  for x in {a..j}
  do
    if [ -d /Volumes/Data/${pool}/${x}/ ]; then
      for cramfile in $(find /Volumes/Data/${pool}/${x}/ -type f -maxdepth 1 -name "*.cram" -exec basename {} \;)
      do
        cramfilename=${cramfile%.*}
        bamfilename=$(awk -v var="$cramfilename" '$1==var{print $2}' ~/COLD/pipeline/1bam_move/cram2bam_name.txt)
        samtools view -b -q 20 -F 0x400 -T ~/COLD/reference/dsimM252v1.2+microbiome.fa -o ~/COLD/data/${bamfilename}/${bamfilename}_${x}.bam /Volumes/Data/${pool}/${x}/${cramfile}
      done
    fi
  done
  samtools merge ~/COLD/data/${bamfilename}/${bamfilename}.bam ~/COLD/data/${bamfilename}/${bamfilename}_*.bam
done

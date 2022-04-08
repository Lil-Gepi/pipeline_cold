#!/bin/bash
samtools --version
bcftools --version
which gawk
which awk

ref_fa="folder/to/dsimM252v1.2+microbiome.fa"

samtools view -T "$ref_fa" -h -F 0x400 "blablabla/a/LB_132.TI02.cram" |\
 head -n 1000000 | cut -f 10 | perl -ne 'chomp;print length($_) . "\n"' |\
  sort | uniq -c

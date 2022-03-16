#!/bin/bash
samtools --version
bcftools --version
which gawk
which awk

for i in {0,10,20,30,40,50,60}; do
  for j in {11..20}; do
    mkdir "/Users/ychen/COLD/data/F${i}_r${j}/"
  done
done

#!/bin/bash
# Specify range in a for loop

cd /data/scratch/DMP/UCEC/UROTRBIO/slise/ROSETREES/DATA/UKB/VCF
dir=/data/scratch/DMP/UCEC/UROTRBIO/slise/DATA/UKBioBank/HAP
for chr in {1..22}
do
  target=${dir}/ukb_hap_chr${chr}_v2.vcf.gz
  link_name=chr${chr}.UKB.vcf.gz
  ls $target
  ln -s $target $link_name
  target=${dir}/ukb_hap_chr${chr}_v2.vcf.gz.tbi
  link_name=chr${chr}.UKB.vcf.gz.tbi
  ls $target
  ln -s $target $link_name
done

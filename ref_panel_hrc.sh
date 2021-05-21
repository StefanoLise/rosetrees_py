#!/bin/bash
# Specify range in a for loop

cd /data/scratch/DMP/UCEC/UROTRBIO/slise/ROSETREES/DATA/HRC/VCF
dir=/data/scratch/DMP/UCEC/UROTRBIO/slise/DATA/HRC/EGAF00001381867
target=${dir}/HRC.r1-1.EGA.GRCh37.chr1.haplotypes.noIBD.vcf.gz
link_name=chr1.HRC.vcf.gz
ls $target
ln -s $target $link_name
dir=/data/scratch/DMP/UCEC/UROTRBIO/slise/DATA/HRC/EGAF00001381868
target=${dir}/HRC.r1-1.EGA.GRCh37.chr1.haplotypes.noIBD.vcf.gz.tbi
link_name=chr1.HRC.vcf.gz.tbi
ls $target
ln -s $target $link_name

num=1381872
for chr in {2..22}
do
  dir=/data/scratch/DMP/UCEC/UROTRBIO/slise/DATA/HRC/EGAF0000${num}
  target=${dir}/HRC.r1-1.EGA.GRCh37.chr${chr}.haplotypes.vcf.gz
  link_name=chr${chr}.HRC.vcf.gz
  ls $target
  ln -s $target $link_name
  ((num++))
  dir=/data/scratch/DMP/UCEC/UROTRBIO/slise/DATA/HRC/EGAF0000${num}
  target=${dir}/HRC.r1-1.EGA.GRCh37.chr${chr}.haplotypes.vcf.gz.tbi
  link_name=chr${chr}.HRC.vcf.gz.tbi
  ls $target
  ln -s $target $link_name
  ((num++))
  ((num++))
  ((num++))
  ((num++))
done

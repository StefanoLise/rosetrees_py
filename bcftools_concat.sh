#!/bin/bash

patID=TR019
fileID=TR019.4x
array=UKB
chrList=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22);

cd  /data/scratch/DMP/UCEC/UROTRBIO/slise/ROSETREES/${patID}/SHAPEIT4
vcfFiles=( ${chrList[@]/#/${fileID}.${array}.chr} );
vcfFiles=( ${vcfFiles[@]/%/.Output.vcf.gz} );
concat_file=${fileID}.${array}.Output.vcf.gz

for x in ${vcfFiles[@]}; do
    echo $x
    if ! [ -e "${x}" ]; then
        echo "Can not find $x"
        exit
    fi
done

if [ -e "${concat_file}" ]; then
    echo "${concat_file} exists already, are you sure you want to overwrite it?"
    exit
fi

bcftools concat -O z -o ${concat_file} ${vcfFiles[@]} ;
tabix -fp vcf ${concat_file} ;
rm ${fileID}.${array}.chr*.Output.vcf.gz*
rm ${fileID}.${array}.chr*.Output.log

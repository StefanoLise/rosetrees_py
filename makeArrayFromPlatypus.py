#!/usr/bin/env python3

import re, gzip, os, subprocess

patID='TR081'
array='UKB'
gqMin=20
phasingMethod='SHAPEIT4'

#############################################################################

def platypusGT():
    sampleID=f'{patID}_GL'
    gtFile=f'/data/scratch/DMP/UCEC/UROTRBIO/slise/ROSETREES/{patID}/VCF/{sampleID}.{array}.vcf.gz'
    try:
        fhand=gzip.open(gtFile,'rt')
    except:
        print('File cannot be opened:', gtFile)
        exit()

    global snp
    snp=[]
    nCalls=nCallsGQ=nGT=0
    for line in fhand:
        line = line.rstrip()
        if line.startswith('#'): continue
        fields = line.split()
        (chrom,pos,varID,ref,alt,qual,varFilter,info,varFormat,sampleGL)=fields[:10] ## double check this !!!!
        chrom=re.sub('^chr','',chrom)
        if (len(ref) > 1) or (len(alt) >1): continue
        (gt,gl,gof,gq,nr,nv)=sampleGL.split(':')
        nCalls += 1
        gq=int(gq)
        if (gq < gqMin): continue
        nCallsGQ += 1
        gt=gt.replace("|","/")
        if (gt == '1/0'): gt='0/1'
        if not ((gt == '0/0') or (gt == '0/1') or (gt == '1/1')): continue
        if ((gt == '0/1') or (gt == '1/1')):
            if (varFilter != 'PASS'): continue
        snp.append([chrom,pos,ref,alt,gt,gq])
        nGT += 1;
    fhand.close()
    print(nCalls,nCallsGQ,nGT)

########################################################################

def outputVCF():
    vcfFile=f'/data/scratch/DMP/UCEC/UROTRBIO/slise/ROSETREES/{patID}/{phasingMethod}/{patID}.{array}.Input.vcf'
    fout=open(vcfFile,'w')
    fout.write('##fileformat=VCFv4.1\n')
    fout.write('##FILTER=<ID=PASS,Description="All filters passed">\n')
    fout.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    fout.write('##FORMAT=<ID=GQ,Number=.,Type=Integer,Description="Genotype quality as phred score">\n')
    fout.write(f'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{patID}\n')
    for var in snp:
        (chrom,pos,ref,alt,gt,gq)=var
        sample=gt+':'+str(gq)
        fout.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (chrom,pos,'.',ref,alt,'.','.','.','GT:GQ',sample))
    fout.close()
    os.system(f'bgzip -f {vcfFile}')
    os.system(f'tabix -fp vcf {vcfFile}.gz')

#######################################################################

platypusGT()
outputVCF()

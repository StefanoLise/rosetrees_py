#!/usr/bin/env python3

import re, gzip, os, subprocess

patID='TR064'
cov='4x'
array='HRC'
afTH=0.005
gtProbTH=0.99
chromList=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22]

glimpseGT=[]

#####################################################################

def run_command(command):
    p = subprocess.Popen(command,
                         shell=True,
                         text=True,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    return iter(p.stdout.readline, b'')

######################################################################

def read_glimpseGT():

    input_dir=f'/data/scratch/DMP/UCEC/UROTRBIO/slise/ROSETREES/{patID}/GLIMPSE'
    for chrom in chromList:
        imputed_file=f'{input_dir}/{patID}.{cov}.chr{chrom}.bcf'
        if not os.path.isfile(imputed_file): exit(f'File not found: {imputed_file}')
        command = f"bcftools view -H -v snps -m2 -M2 {imputed_file} -r {chrom}"
        for line in run_command(command):
            if line == '':  break
            line = line.rstrip()
            if line.startswith('#'): continue
            fields = line.split()
            if (len(fields) != 10): exit()
            (chr,pos,id,ref,alt,qual,filter,info,format,sample)=fields[:]
            (gt,ds,gp,hs)=sample.split(':')
            gtProb=gp.split(',')
            raf=None
            for infoField in info.split(';'):
                matchObj=re.search('^RAF=(.*)',infoField)
                if matchObj: raf=matchObj.group(1)
                break
            if raf == None : exit('Raf not found'+line+'\n')
            raf=float(raf)
#            if (raf < afTH) or (raf > 1-afTH): continue
            if ((len(ref) > 1) or (len(alt) >1)): continue
            if not re.search('^[ACGT]$',ref) and re.search('^[ACGT]$',alt): continue
            if (gt == '1/0'): gt='0/1'
            if not ((gt == '0/0') or (gt == '0/1') or (gt == '1/1')): continue
            iP=0 if (gt == '0/0') else 1 if (gt == '0/1') else 2
            if float(gtProb[iP]) < gtProbTH: continue 
            glimpseGT.append([chrom,pos,ref,alt,gt])
                
##########################################################################################

def outputVCF():
    output_dir=f'/data/scratch/DMP/UCEC/UROTRBIO/slise/ROSETREES/{patID}/SHAPEIT4'
    vcf_file=f'{output_dir}/{patID}.{cov}.{array}.Input.vcf'
    try:
        fout=open(vcf_file,'w')
    except:
        exit("Can not open: "+vcf_file)
    fout.write('##fileformat=VCFv4.1\n')
    fout.write('##FILTER=<ID=PASS,Description="All filters passed">\n')
    fout.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    fout.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'+patID+'\n')
    for var in glimpseGT:
        (chrom,pos,ref,alt,gt)=var
        fout.write('\t'.join(str(p) for p in (chrom,pos,'.',ref,alt,'.','.','.','GT',gt))+'\n')
#        fout.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (chrom,pos,'.',ref,alt,'.','.','.','GT',gt))
    fout.close()
    os.system(f'bgzip -f {vcf_file}')
    os.system(f'tabix -fp vcf {vcf_file}.gz')

#############################################################################################

read_glimpseGT()
outputVCF()






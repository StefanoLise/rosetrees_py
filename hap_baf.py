#!/usr/bin/env python3

import re, gzip, os, subprocess, statistics, sys

#patID='TR064'                   # These is the patient ID  
patID=str(sys.argv[1])
covGL=4                          # This determine if the haplotype is from low pass or high pass 
covSample=None
tc=None
samples=('BL.8x','PD1','PD2')               # These are the samples for which we will calculate the HAP BAF   
array='HRC'
phasingMethod='SHAPEIT4'
readCountMethod='bcftools'       # could be either platypus or bcftools
afRefPanelTH=0.005
binType='by_snp'        # could be 'by_snp' or 'by_bp'
if binType == 'by_snp':
#    n_snp_bin=2000
#    n_snp_min=1500
    n_snp_bin=int(sys.argv[2])
    n_snp_min=int(3*n_snp_bin/4)
    max_interval_size=150000000
elif  binType == 'by_bp':
    pass
else:
    exit(f"Bin type '{binType}' not defined: it's either 'by_snp' or 'by_bp'")

chrList=(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22)

allele_freq=dict()
phased_hap=[]
chrom_bins=[]
snp_count=[]

##############################################################

def run_command(command):
    p = subprocess.Popen(command,
                         shell=True,
                         text=True,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    return iter(p.stdout.readline, b'')

##############################################################

def read_allele_freq(chrom):
    vcfFile='/data/scratch/DMP/UCEC/UROTRBIO/slise/ROSETREES/DATA/SNPs/HRC/snps.vcf.gz'
    allele_freq={}
    command = f"bcftools query -f'%POS\t%INFO/AF\n' {vcfFile} -r {chrom}"
    for line in run_command(command):
        if line == '':  break
        line = line.rstrip()
        if line.startswith('#'): continue
        fields = line.split()
        (pos,af)=fields[:]
        pos=int(pos)
        af=float(af)
        allele_freq[pos]=af
    return allele_freq


##############################################################

def centromere_positions(chrom):
    centromere=dict()
    centromere[1]=[121535434,124535434]
    centromere[2]=[92326171,95326171]
    centromere[3]=[90504854, 93504854]
    centromere[4]=[49660117, 52660117]
    centromere[5]=[46405641, 49405641]
    centromere[6]=[58830166, 61830166]
    centromere[7]=[58054331, 61054331]
    centromere[8]=[43793000, 46857000]
    centromere[9]=[47367679, 50367679]
    centromere[10]=[39254935, 42254935]
    centromere[11]=[51644205, 54644205]
    centromere[12]=[34856694, 37856694]
    centromere[13]=[16000000, 19000000]
    centromere[14]=[16000000, 19000000]
    centromere[15]=[17000000, 20000000]
    centromere[16]=[35335801, 38335801]
    centromere[17]=[22263006, 25263006]
    centromere[18]=[15460898, 18460898]
    centromere[19]=[24681782, 27681782]
    centromere[20]=[26369569, 29369569]
    centromere[21]=[11288129, 14288129]
    centromere[22]=[13000000, 16000000]
    return centromere[chrom]

###########################################################################################

def read_phased_hap(chrom):
    phased_hap_file=f'/data/scratch/DMP/UCEC/UROTRBIO/slise/ROSETREES/{patID}/{phasingMethod}/' \
                     f'{patID}.HRC.UKB.Output.vcf.gz';
    if covGL:
        phased_hap_file=f'/data/scratch/DMP/UCEC/UROTRBIO/slise/ROSETREES/{patID}/{phasingMethod}/' \
                         f'{patID}.{covGL}x.HRC.UKB.Output.vcf.gz';
    if not os.path.isfile(phased_hap_file): exit(f'File not found: {phased_hap_file}')
    phased_hap=[]
    command = f"bcftools view -H -v snps -m2 -M2 {phased_hap_file} -r {chrom}"
    for line in run_command(command):
        if line == '':  break
        line = line.rstrip()
        if line.startswith('#'): continue
        fields = line.split()
        if (len(fields) != 10): exit()
        (chr,pos,id,ref,alt,qual,filter,info,format,sample)=fields[:]
        pos=int(pos)
        af=allele_freq.get(pos,0)
        if af < afRefPanelTH or af > 1-afRefPanelTH: continue
        gt=sample
        if ((len(ref) > 1) or (len(alt) >1)): continue
        if not re.search('^[ACGT]$',ref) and re.search('^[ACGT]$',alt): continue
        if ((gt == '0|1') or (gt == '1|0')): 
            phased_hap.append([pos,gt])
    return phased_hap

###############################################################

def split_chromosome_into_bins(chrom):
    (centro_start,centro_end)=centromere_positions(chrom)
    n_snp_tot=len(phased_hap)
    chrom_bins=[]
    i_start=0
    while (i_start < n_snp_tot):
        pos=phased_hap[i_start][0]
        while centro_start <= pos <= centro_end:
            i_start += 1
            pos=phased_hap[i_start][0]
        i_max=i_start+n_snp_bin
        if (i_max > n_snp_tot): i_max=n_snp_tot
        pos_max=phased_hap[-1][0]
        if phased_hap[i_start][0] < centro_start:
            pos_max=centro_start
        for i in (range(i_start,i_max)):
            i_end=i
            if i_end == n_snp_tot-1: break
            pos=phased_hap[i_end+1][0]
            if pos > pos_max: break
        chrom_bins.append([i_start,i_end])
        i_start = i_end + 1
    return chrom_bins

###############################################################

def read_count(sample_id,chrom):
    vcf_file=f'/data/scratch/DMP/UCEC/UROTRBIO/slise/ROSETREES/{patID}/VCF/{sample_id}.{array}.vcf.gz'
    if not os.path.isfile(vcf_file): exit(f'File not found: {vcf_file}')
    phased_hap_file=f'/data/scratch/DMP/UCEC/UROTRBIO/slise/ROSETREES/{patID}/{phasingMethod}/' \
                     f'{patID}.HRC.UKB.Output.vcf.gz';
    if covGL:
        phased_hap_file=f'/data/scratch/DMP/UCEC/UROTRBIO/slise/ROSETREES/{patID}/{phasingMethod}/' \
                         f'{patID}.{covGL}x.HRC.UKB.Output.vcf.gz';
    if not os.path.isfile(phased_hap_file): exit(f'File not found: {phased_hap_file}')
    snp_count=[]
    if readCountMethod == 'platypus': 
        command = f'bcftools annotate -a {vcf_file} {phased_hap_file} -c FORMAT/NR,FORMAT/NV -r {chrom} | ' \
                   f'bcftools view -g het -H'
    elif readCountMethod == 'bcftools':
        command = f'bcftools annotate -a {vcf_file} {phased_hap_file} -x INFO,^FORMAT/GT -c INFO/DP4 -r {chrom} | ' \
                   f'bcftools view -g het -H'
    for line in run_command(command):
        if line == '':  break
        line = line.rstrip()
        if line.startswith('#'): continue
        fields = line.split()
        if (len(fields) != 10): exit()
        (pos,id,ref,alt,qual,filter,info,format,format_values)=fields[1:]
        pos=int(pos)
        af=allele_freq.get(pos,0)
        if af < afRefPanelTH or af > 1-afRefPanelTH: continue
        if readCountMethod == 'platypus':
            if (len(format.split(':')) == 3):
                (gt,nr,nv)=format_values.split(':')
                if not nr.isnumeric(): nr=0
                if not nv.isnumeric(): nv=0
            elif (len(format.split(':')) == 1):
                (gt,nr,nv)=(format_values,0,0)
            else:
                exit(f'Something went wrong with the annotation of {sample_id} at {chrom}:pos')
            n_ref=int(nr)-int(nv)
            n_alt=int(nv)
        elif readCountMethod == 'bcftools': 
            gt=format_values
            (refF,refR,altF,altR)=(0,0,0,0)
            dp4=re.search("^DP4=(\d+),(\d+),(\d+),(\d+)$",info)
            if dp4:
                (refF,refR,altF,altR)=map(int,dp4.group(1,2,3,4))
            n_ref=refF+refR
            n_alt=altF+altR        
        if not ((gt == '0|1') or (gt == '1|0')): continue
        if (gt == '0|1'):
            snp_count.append([pos,gt,n_ref,n_alt])
        elif (gt == '1|0'):
            snp_count.append([pos,gt,n_alt,n_ref])
    return snp_count

#####################################################################

def sanity_checks():
    if (len(phased_hap) != len(snp_count)):
        print (f'The length of phased_hap is:',len(phased_hap))
        print (f'The length of snp_count is:',len(snp_count))
        exit("The arrays 'phased_hap' and 'snp_count' shouldn't have different lengths\n")

#####################################################################

def calculate_hap_baf(chrom):
    results=[]
    for bin in chrom_bins:
        (i_start,i_end)=bin
        n_snp=i_end+1-i_start
        pos_start=phased_hap[i_start][0]
        pos_end=phased_hap[i_end][0]
        int_size=pos_end + 1 - pos_start
        n_l=n_r=n_snp=0
        for i in (range(i_start,i_end+1)):
            n_snp +=1
            n_l += snp_count[i][2]
            n_r += snp_count[i][3]
        af_l=n_l/(n_l+n_r)
        af_r=n_r/(n_l+n_r)
        daf_bin=af_l-af_r
        af_l=f"{af_l:.4f}"
        af_r=f"{af_r:.4f}"
        daf_bin=f"{daf_bin:.4f}"
        results.append([chrom,pos_start,pos_end,n_snp,int_size,af_l,af_r,daf_bin])
    return results

##############################################################

def output_file(sample_id,results):
    out_dir=f'/data/scratch/DMP/UCEC/UROTRBIO/slise/ROSETREES/{patID}/HAP_BAF'
    if binType == 'by_snp':
        if n_snp_bin >= 1000:
            int_size_f=f"{n_snp_bin//1000:d}Ksnp"
        else:
            int_size_f=f"{n_snp_bin}snp"
    elif  binType == 'by_bp':
        pass
    else:
        exit('Problems')
    out_file=f'{out_dir}/{sample_id}.{array}.{phasingMethod}.{int_size_f}.tsv';
    if covGL:
        out_file=f'{out_dir}/{sample_id}.{array}.GLIMPSE.{int_size_f}.tsv';

    try:
        fout= open(out_file,'w')
    except:
        exit(f"Can not open: {out_file}")
    print('chrom','bin_start','bin_end','n_snp','int_size','af_l','af_r','daf',sep="\t",file=fout)
    hap_bin_size=[]
    daf_bin_list=[]
    i_bin=0
    for out in results:
        (chrom,pos_start,pos_end,n_snp,int_size,af_l,af_r,daf_bin)=out
        print(*out,sep="\t",file=fout)
        if (n_snp > n_snp_min) and (int_size < max_interval_size):
            hap_bin_size.append(int_size)
            daf_bin=float(daf_bin)
            daf_bin_list.append(daf_bin)
#            if daf_bin > 0.03 or daf_bin < -0.03:
#                print(*[i_bin,*out],sep="\t")
            i_bin += 1
    fout.close()
    print(statistics.mean(hap_bin_size),statistics.stdev(hap_bin_size),max(hap_bin_size),
          min(hap_bin_size),len(hap_bin_size))
    print(statistics.quantiles(hap_bin_size,n=10))
    print(statistics.mean(daf_bin_list),statistics.stdev(daf_bin_list),max(daf_bin_list),
          min(daf_bin_list),len(daf_bin_list))
    print(statistics.quantiles(daf_bin_list,n=10))

###############################################################

for s_t in samples:
    sample_id=f'{patID}_{s_t}'
    if covSample :
        sample_id = sample_id + f'.{covSample}x' 
    if tc :
        sample_id = sample_id + f'.TC{tc}'
    print(sample_id,f"{n_snp_bin=}",f"{n_snp_min=}")
    results=[]
    for chrom in chrList:
        print(chrom)
        allele_freq=read_allele_freq(chrom)
        phased_hap=read_phased_hap(chrom)
        chrom_bins=split_chromosome_into_bins(chrom)
        snp_count=read_count(sample_id,chrom)
        sanity_checks()
        results.extend(calculate_hap_baf(chrom))
    output_file(sample_id,results)
   

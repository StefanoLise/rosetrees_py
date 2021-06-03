#!/usr/bin/env python3


import re, gzip, os, subprocess, statistics, random

patID='V5322'

target_cov=60
target_tc=0.50
cov_t=37.9                  # coverage of the original tumour bam
cov_n=25.1                  # coverage of the original GL bam
tc=0.70                     # tumour content of the original tumour bam
k_2=1.104                   # fraction of reads (N_READS_GL/N_READS_BL) in diploid regions in both samples
target_tc_f=int(100*target_tc)
target_tc_f=f'{target_tc_f:02d}'


bam_dir=f'/data/scratch/DMP/UCEC/UROTRBIO/slise/ROSETREES/{patID}/BAM'
t_bam=f'{bam_dir}/{patID}_BL.bam'
n_bam=f'{bam_dir}/{patID}_GL.bam'
t_bam_tmp=f'{bam_dir}/{patID}_BL.{target_cov}x.TC{target_tc_f}.tmp.bam'
n_bam_tmp=f'{bam_dir}/{patID}_GL.{target_cov}x.TC{target_tc_f}.tmp.bam'
out_bam=f'{bam_dir}/{patID}_BL.{target_cov}x.TC{target_tc_f}.bam'

tc_r=target_tc/tc
f_t=target_cov*(k_2*tc_r)/((1-tc_r)*cov_n+k_2*tc_r*cov_t)
if float(f_t) > 1:
    exit(f"f_t = {f_t} - it needs to be < 1")
f_t=int(100000*f_t)
f_t=f'{f_t:05d}'
f_n=target_cov*(1-tc_r)/((1-tc_r)*cov_n+k_2*cov_t*tc_r)
if float(f_n) > 1:
    exit(f"f_n = {f_n} - it needs to be < 1")
f_n=int(100000*f_n)
f_n=f'{f_n:05d}'

print(f'Sample {patID}_BL.{target_cov}x.TC{target_tc_f}')
print(f'Target coverage: {target_cov}x')
print(f'Target tumor content: {target_tc}')
print(f'Fraction of reads kept in GL: 0.{f_n}')
print(f'Fraction of reads kept in BL: 0.{f_t}')
iseed=random.randint(1,1000000);
s = f"{iseed}.{f_t}";
print(f'Paramater -s for subsampling BL: {s}')
os.system(f'samtools view -bh {t_bam} -s {s} -o {t_bam_tmp}')
os.system(f'samtools index {t_bam_tmp}')

iseed=random.randint(1,1000000);
s = f"{iseed}.{f_n}";
print(f'Paramater -s for subsampling GL: {s}')
os.system(f'samtools view -bh {n_bam} -s {s} -o {n_bam_tmp}')
os.system(f'samtools index {n_bam_tmp}')

os.system(f'samtools merge -f {out_bam}  {n_bam_tmp} {t_bam_tmp}')
os.system(f'samtools index {out_bam}')

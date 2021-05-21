#!/usr/bin/env python3

# This code calculates the read count from a bam file at pre-defined {sites}.
# The original bam file can be downsampled to a target coverage 
#

import gzip, re, os, subprocess, random

pat_id='TR081'
s_t='BL'
target_cov=8
tc=None
sample_id=f'{pat_id}_{s_t}'
array='HRC'
cov_bam={'PT23_GL':35.7,'PT23_BL':86,  'PT23_PD':96.1,
         'PT36_GL':42.6,'PT36_BL':111, 'PT36_PD':113,
         'PT78_GL':36.1,'PT78_BL':74.7,'PT78_PD':78.5,
         'V5322_GL':35.5,'V5322_BL':59.3,'V5322_PD':'???',
         'TR081_GL':25.1,'TR081_BL':37.9,'TR081_PD1':5.8,'TR081_PD2':5.4,
         'TR019_GL':24.1,'TR019_BL':26.3,'TR019_PD1':4,'TR019_PD2':4,
         'TR064_GL':22.6,'TR064_BL':28.1,'TR064_PD1':4,'TR064_PD2':4}

bam_dir=f'/data/scratch/DMP/UCEC/UROTRBIO/slise/ROSETREES/{pat_id}/BAM'
bam_file=f'{bam_dir}/{sample_id}.bam'
sites=f'/data/scratch/DMP/UCEC/UROTRBIO/slise/ROSETREES/DATA/SNPs/{array}/snps'
if not os.path.isdir(f'{bam_dir}'): exit(f'Directory not found: {bam_dir}')
if not os.path.isfile(f'{bam_file}'): exit(f'File not found: {bam_file}')
if not os.path.isfile(f'{sites}.vcf.gz'): exit(f'File not found: {sites}.vcf.gz')
if not os.path.isfile(f'{sites}.tsv.gz'): exit(f'File not found: {sites}.tsv.gz')

read_count_bam=bam_file
if target_cov < cov_bam[sample_id]:
    read_count_bam=f'{bam_dir}/{sample_id}.{target_cov}x.bam'
    if tc :
        read_count_bam=f'{bam_dir}/{sample_id}.{target_cov}x.TC{tc}.bam'
    if os.path.isfile(f'{read_count_bam}'): 
        print(f'File exists already: {read_count_bam} - moving on to read counts')
    elif tc :
        exit(f'Can not create a diluted bam file on the fly')
    else:
        print(f'Sub-sampling {sample_id}.bam')
        iseed=random.randint(1,1000000);
        frac=int(100000*target_cov/cov_bam[sample_id])
        frac=f"{frac:05d}"
        s = f'{iseed}.{frac}'
        print(s)
        os.system(f'samtools view -bh {bam_file} -s {s} -o {read_count_bam}')
        os.system(f'samtools index {read_count_bam}')
        print('Done with subsampling')

out_dir=f'/data/scratch/DMP/UCEC/UROTRBIO/slise/ROSETREES/{pat_id}/VCF'
out_file=f'{out_dir}/{sample_id}.{array}.vcf.gz'
if target_cov < cov_bam[sample_id]:
    out_file=f'{out_dir}/{sample_id}.{target_cov}x.{array}.vcf.gz'
    if tc :
        out_file=f'{out_dir}/{sample_id}.{target_cov}x.TC{tc}.{array}.vcf.gz'

ref_file='/data/scratch/DMP/UCEC/UROTRBIO/slise/DATA/GENOMES/HS37D5/hs37d5.fa'
os.system(f"bcftools mpileup {read_count_bam} --ignore-RG -f {ref_file} -I -E -T {sites}.vcf.gz -d 10000 -Ou |" \
          + f" bcftools call -Aim -C alleles -T {sites}.tsv.gz -Oz -o {out_file}")
os.system(f"tabix -fp vcf {out_file}");
print('Done')

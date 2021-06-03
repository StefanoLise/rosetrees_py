#!/usr/bin/env python3

import re, gzip, os, subprocess

pat_id='TR019'
cov_gl=4                          # This determine if the haplotype is from low pass or high pass
bl_sample=f'{pat_id}_BL.8x'
sample_test=['BL.8x','PD1','PD2']
array='HRC'
phasing_method='SHAPEIT4'
read_count_method_bl='bcftools'   # could be either platypus or bcftools
read_count_method_test='bcftools'   # could be either platypus or bcftools
af_ref_panel_th=0.005
bin_type='by_snp'        # could be 'by_snp' or 'by_bp'
if bin_type == 'by_snp':
    n_snp_bin=1000
    n_snp_min=750
    max_interval_size=1500000
    gl_stdev=0.02
    if n_snp_bin >= 1000:
        int_size_f=f"{n_snp_bin//1000:d}Ksnp"
    else:
        int_size_f=f"{n_snp_bin}snp"
elif bin_type == 'by_bp':
    pass
else:
    exit(f"Bin type '{bin_type}' not defined: it's either 'by_snp' or 'by_bp'")

loh_regions=dict()
allele_freq=dict()
phased_hap=[]
chrom_bins=[]


################################################################

def run_command(command):
    p = subprocess.Popen(command,
                         shell=True,
                         text=True,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    return iter(p.stdout.readline, b'')

##################################################

def read_regions_and_baf(bl_sample):
    if bl_sample.startswith('PT23_BL'):
        loh_regions['5:80000000-180000000']=0.14
        loh_regions['6:63000000-166000000']=0.21
        loh_regions['8:1-43000000']=0.21
        loh_regions['13:20000000-77000000']=0.21
        loh_regions['16:49000000-90000000']=0.21
    elif bl_sample.startswith('PT23_PD'):
        loh_regions['5:80000000-180000000']=0.19
        loh_regions['6:63000000-166000000']=0.27
        loh_regions['8:1-43000000']=0.27
        loh_regions['13:20000000-77000000']=0.27
        loh_regions['16:49000000-90000000']=0.27
    elif bl_sample.startswith('PT36_BL'):
        pass
    elif bl_sample.startswith('PT36_PD'):
        pass
    elif bl_sample.startswith('PT78_BL'):
        pass
    elif bl_sample.startswith('PT78_PD'):
        pass
    elif bl_sample.startswith('V5322_BL'):
        loh_regions['1:69000000-89000000']=0.06
        loh_regions['4:80000000-175000000']=0.34
        loh_regions['6:1-29000000']=0.34
        loh_regions['8:1-43000000']=0.06
        loh_regions['9:1-35000000']=0.06
        loh_regions['13:41000000-70000000']=0.11
        loh_regions['16:53500000-90000000']=0.27
        loh_regions['18:35500000-78000000']=0.06
    elif bl_sample.startswith('TR081_BL'):
        loh_regions['8:1-43000000']=0.15
        loh_regions['12:1-35000000']=0.11
        loh_regions['17:1-22000000']=0.15
    elif bl_sample.startswith('TR064_BL'):
#        loh_regions['8:1-40000000']=0.465
        loh_regions['13:49000000-90500000']=0.42
    elif bl_sample.startswith('TR019_BL'):
        loh_regions['6:84000000-122000000']=0.465
        loh_regions['8:8000000-43000000']=0.465
        loh_regions['11:98000000-130000000']=0.435
    else:
        exit(f'Can not recognize sample {bl_sample}')

##################################################

def parse_region(region):
    region_obj=re.search("^(\d+):(\d+)-(\d+)",region)
    if not region_obj: exit(f"Not an acceptable region: {region}")
    (chrom,region_start,region_end)=map(int,region_obj.group(1,2,3))
    return chrom,region_start,region_end

##################################################
def read_allele_freq(chrom):
    vcfFile='../DATA/SNPs/HRC/snps.vcf.gz'
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

################################################################################

def read_phased_hap(chrom):
    phased_hap_file=f'../{pat_id}/{phasing_method}/{pat_id}.HRC.UKB.Output.vcf.gz';
    if cov_gl:
        phased_hap_file=f'../{pat_id}/{phasing_method}/{pat_id}.{cov_gl}x.HRC.UKB.Output.vcf.gz';
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
        if af < af_ref_panel_th or af > 1-af_ref_panel_th: continue
        gt=sample
        if ((len(ref) > 1) or (len(alt) >1)): continue
        if not re.search('^[ACGT]$',ref) and re.search('^[ACGT]$',alt): continue
        if ((gt == '0|1') or (gt == '1|0')):
            phased_hap.append([pos,ref,alt,gt])
    return phased_hap

#################################################

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

#################################################################

def read_count(sample_id,chrom):
    vcf_file=f'../{pat_id}/VCF/{sample_id}.{array}.vcf.gz'
    if not os.path.isfile(vcf_file): exit(f'File not found: {vcf_file}')
    phased_hap_file=f'../{pat_id}/{phasing_method}/{pat_id}.HRC.UKB.Output.vcf.gz';
    if cov_gl:
        phased_hap_file=f'../{pat_id}/{phasing_method}/{pat_id}.{cov_gl}x.HRC.UKB.Output.vcf.gz';
    if not os.path.isfile(phased_hap_file): exit(f'File not found: {phased_hap_file}')
    allele_count=[]
    if read_count_method_bl == 'platypus':
        command = f'bcftools annotate -a {vcf_file} {phased_hap_file} -c FORMAT/NR,FORMAT/NV -r {chrom} | ' \
                   f'bcftools view -g het -H'
    elif read_count_method_bl == 'bcftools':
        command = f'bcftools annotate -a {vcf_file} {phased_hap_file} -x INFO,^FORMAT/GT -c INFO/DP4 -r {chrom} | ' \
                   f' bcftools view -g het -H'
    for line in run_command(command):
        if line == '':  break
        line = line.rstrip()
        if line.startswith('#'): continue
        fields = line.split()
        if (len(fields) != 10): exit()
        (pos,id,ref,alt,qual,filter,info,format,format_values)=fields[1:]
        pos=int(pos)
        af=allele_freq.get(pos,0)
        if af < af_ref_panel_th or af > 1-af_ref_panel_th: continue
        if read_count_method_bl == 'platypus':
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
        elif read_count_method_bl == 'bcftools':
            gt=format_values
            (refF,refR,altF,altR)=(0,0,0,0)
            dp4=re.search("^DP4=(\d+),(\d+),(\d+),(\d+)$",info)
            if dp4:
                (refF,refR,altF,altR)=map(int,dp4.group(1,2,3,4))
            n_ref=refF+refR
            n_alt=altF+altR
        if not ((gt == '0|1') or (gt == '1|0')): continue
        if (gt == '0|1'):
            allele_count.append([pos,gt,n_ref,n_alt])
        elif (gt == '1|0'):
            allele_count.append([pos,gt,n_alt,n_ref])
    return allele_count

#####################################################################

def sanity_checks():
    if (len(phased_hap) != len(allele_count_bl)):
        print (f'The length of phased_hap is:',len(phased_hap))
        print (f'The length of allele_count is:',len(allele_count_bl))
        exit("The arrays 'phased_hap' and 'allele_count' shouldn't have different lengths\n")

#####################################################################

def switch_gt(gt):
    if (gt == '0|1'):
        gt_new = '1|0'
    elif (gt == '1|0'):
        gt_new = '0|1'
    else:
        exit(f'{gt} is not heterozygous')
    return gt_new

####################################################################

def calculate_bin_baf(chrom_bins):
    baf_region=loh_regions[region]
    baf_error=3*gl_stdev
    baf_th=baf_region+baf_error
    chrom_bins_copy=chrom_bins.copy()
    for i_bin in range(len(chrom_bins_copy)):
        (i_start,i_end)=chrom_bins_copy[i_bin]
        n_l=n_r=n_snp=0
        for i in (range(i_start,i_end+1)):
            n_snp +=1
            n_l += allele_count_bl[i][2]
            n_r += allele_count_bl[i][3]
        af_l=n_l/(n_l+n_r)
        flag  ='+' if af_l <= baf_th else '-' if af_l >= (1-baf_th) else '?'
        if n_snp < n_snp_min: flag='?'
        chrom_bins[i_bin].extend([n_l,n_r,flag])

#################################################################

def hap_building(region):
    (chrom,region_start,region_end)=parse_region(region)
    out_dir=f'../{pat_id}/HAP_BUILDING'
    out_vcf_file=f'{out_dir}/{bl_sample}.chr{region}.{int_size_f}.vcf'
    try:
        fout_vcf=open(out_vcf_file,'w')
    except:
        exit(f"Can not open: {out_vcf_file}")
    fout_vcf.write('##fileformat=VCFv4.1\n')
    fout_vcf.write('##FILTER=<ID=PASS,Description="All filters passed">\n')
    fout_vcf.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    fout_vcf.write(f'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{pat_id}\n')
    n_snp_region=n_snp_hap=0
    bin_baf=[]
    for bin in chrom_bins:
        (i_start,i_end,n_l,n_r,flag)=bin
        pos_start=phased_hap[i_start][0]
        pos_end=phased_hap[i_end][0]
        if (pos_end < region_start or pos_start > region_end): continue
        n_snp = (i_end + 1 - i_start)
        n_snp_region += n_snp
        int_size=pos_end + 1 - pos_start
        bin_baf.append([chrom,pos_start,pos_end,n_snp,int_size,n_l,n_r,flag])
        if flag == '?' :continue
        n_snp_hap += n_snp
        for i in (range(i_start,i_end+1)):
            (pos,ref,alt,gt)=phased_hap[i]
            if flag == '-': gt=switch_gt(gt)
            fout_vcf.write('\t'.join(str(p) for p in ([chrom,pos,'.',ref,alt,'.','PASS','.','GT',gt]))+'\n')
    fout_vcf.close()
    os.system(f'bgzip -f {out_vcf_file}')
    os.system(f'tabix -fp vcf {out_vcf_file}.gz')
    print(region,n_snp_region,n_snp_hap)

    out_baf_file=f'{out_dir}/{bl_sample}.chr{region}.{int_size_f}.baf.tsv'
    try:
        fout_baf=open(out_baf_file,'w')
    except:
        exit(f"Can not open: {out_baf_file}")
    print('chrom','bin_start','bin_end','n_snp','int_size','n_l','n_r','af_l','af_r','daf',sep="\t",file=fout_baf)
    for bin in bin_baf:
        (chrom,pos_start,pos_end,n_snp,int_size,n_l,n_r,flag)=bin
        af_l=n_l/(n_l+n_r)
        af_r=n_r/(n_l+n_r)
        daf_bin=af_l-af_r
        af_l=f"{af_l:.4f}"
        af_r=f"{af_r:.4f}"
        daf_bin=f"{daf_bin:.4f}"
        print(*[chrom,pos_start,pos_end,n_snp,int_size,n_l,n_r,af_l,af_r,daf_bin,flag],sep="\t",file=fout_baf)
    fout_baf.close()

#################################################################

def tc_est(region,sample_id):
    (chrom,region_start,region_end)=parse_region(region)
    hap_file=f'../{pat_id}/HAP_BUILDING/{bl_sample}.chr{region}.{int_size_f}.vcf.gz'
    vcf_file=f'../{pat_id}/VCF/{sample_id}.{array}.vcf.gz'
    if not os.path.isfile(hap_file): exit(f'File not found: {hap_file}')
    if not os.path.isfile(vcf_file): exit(f'File not found: {vcf_file}')

    n_l=n_r=n_snp=0;
    if read_count_method_test == 'platypus':
        command = f'bcftools annotate -a {vcf_file} {hap_file} -c FORMAT/NR,FORMAT/NV -r {region} |' \
                   f'bcftools view -g het -H'
    elif read_count_method_test == 'bcftools':
        command = f'bcftools annotate -a {vcf_file} {hap_file} -x INFO,^FORMAT/GT -c INFO/DP4 -r {chrom} | ' \
                   f'bcftools view -g het -H'
    for line in run_command(command):
        if line == '':  break
        line = line.rstrip()
        if line.startswith('#'): continue
        fields = line.split()
        if (len(fields) != 10): exit()
        (pos,var_id,ref,alt,qual,var_filter,info,format_fields,format_values)=fields[1:]
        if ((len(ref) > 1) or (len(alt) >1)): continue
        if not (re.search('^[ACGT]$',ref) and re.search('^[ACGT]$',alt)):
            print(line)
            continue
        if read_count_method_test == 'platypus':
            if (len(format_fields.split(':')) == 3):
                (gt,nr,nv)=format_values.split(':')
                if not nr.isnumeric(): nr=0
                if not nv.isnumeric(): nv=0
            elif (len(format_values.split(':')) == 1):
                (gt,nr,nv)=(format_values,0,0)
            else:
                exit(f'Something went wrong with the annotation of {sample_id} at {chrom}:{pos}')
            n_ref=int(nr)-int(nv)
            n_alt=int(nv)
        elif read_count_method_test == 'bcftools':
            gt=format_values
            (ref_f,ref_r,alt_f,alt_r)=(0,0,0,0)
            dp4=re.search("^DP4=(\d+),(\d+),(\d+),(\d+)$",info)
            if dp4:
                (ref_f,ref_r,alt_f,alt_r)=map(int,dp4.group(1,2,3,4))
            n_ref=ref_f+ref_r
            n_alt=alt_f+alt_r
        if not ((gt == '0|1') or (gt == '1|0')): continue
        n_snp += 1
        if (gt == '0|1'):
            n_l += n_ref
            n_r += n_alt
        elif (gt == '1|0'):
            n_l += n_alt
            n_r += n_ref
    if (n_l+n_r):
        af_l=f"{n_l/(n_l+n_r):.4f}"
        af_r=f"{n_r/(n_l+n_r):.4f}"
        daf=f"{(n_r-n_l)/(n_l+n_r):.4f}"
    else:
        af_l='NA'
        af_r='NA'
        daf='NA'
    print(f'{sample_id}',region,n_snp,n_l,n_r,n_l+n_r,af_l,af_r,daf,sep='\t')

#########################################################################################################

read_regions_and_baf(bl_sample)
for region in loh_regions:
    (chrom,region_start,region_end)=parse_region(region)
    allele_freq=read_allele_freq(chrom)
    phased_hap=read_phased_hap(chrom)
    chrom_bins=split_chromosome_into_bins(chrom)
    allele_count_bl=read_count(bl_sample,chrom)
    sanity_checks()
    calculate_bin_baf(chrom_bins)
    hap_building(region)
    for s_t in sample_test:
        sample_id=f'{pat_id}_{s_t}'
        tc_est(region,sample_id)

#!/usr/bin/env python3

# This script analyses the phasing accuracya of a reconstructing haplotype, e.g. using GLIMPSE for
#   imputation and SHAPEIT4 for phasing. In particular it determines the phasing acciracy over bins 
#   of either (i) a predetermined genomic size {hapSize} or (ii) with a predetrmined number of het 
#   snps.
#

import gzip, re, os, subprocess, statistics

sampleID='NA12878'
confRegionsFlag=None
afRefPanelTH=0.005
bam='GIAB'
array='HRC'
covGL=4
binType='by_snp'        # could be 'by_snp' or 'by_bp'
if binType == 'by_snp':
    n_snp_bin=2000
    n_snp_min=1500
    max_interval_size=150000000
elif  binType == 'by_bp':
    hapSize=5000000 
else:
    exit(f"Bin type '{binType}' not defined: it's either 'by_snp' or 'by_bp'")

chrList=(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22)

allele_freq=dict()
giab_hap=dict()
phased_hap=dict()


centromere=dict()

#####################################################################

def run_command(command):
    p = subprocess.Popen(command,
                         shell=True,
                         text=True,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    return iter(p.stdout.readline, b'')

###############################################################################

def centromere_positions():
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

###############################################################################

def read_allele_freq(chrom):
    vcfFile='../REF_PANEL/HRC/VCF/sites.vcf.gz' 
    allele_freq={}

    n=nCommon=0
    command = f"bcftools query -f'%POS\t%INFO/AF\n' {vcfFile} -r {chrom}"
    for line in run_command(command):
        if line == '':  break
        line = line.rstrip()
        if line.startswith('#'): continue
        fields = line.split()
        (pos,af)=fields[:]
        af=float(af)
        allele_freq[int(pos)]=af
        if afRefPanelTH <= af <= (1-afRefPanelTH): nCommon += 1
    n=len(allele_freq)
    print("Number of genotyped positions in the reference panel: {0:d}, of which {1:d} common".format(n,nCommon))
    return allele_freq

#############################################################################

def read_giab_hap(chrom):
    hapFile=f'/scratch/DMP/BIOCHLEU/slise/LowPassWGS/NA12878/VCF/GIAB/snps.{array}.vcf.gz'
    if not os.path.isfile(hapFile): exit(f'File not found: f{hapFile}')

    n=nHet=0
    giab_hap=dict()
    command = f"bcftools view -H -v snps -m2 -M2 {hapFile} -r {chrom}"
    if confRegionsFlag:
        confRegions=f'../NA12878/VCF/GIAB/confRegions.chr{chrom}.bed'
        if not os.path.isfile(confRegions): exit(f'File not found: {confRegions}')
        command = f"bcftools view -H -v snps -m2 -M2 {hapFile} -R {confRegions}"
    for line in run_command(command):
        if line == '':  break
        line = line.rstrip()
        if line.startswith('#'): continue
        fields = line.split()
        if (len(fields) != 10): exit()
        (chr,pos,id,ref,alt,qual,filter,info,format,gt)=fields[:]
        pos=int(pos)
        if ((len(ref) > 1) or (len(alt) >1)): continue
        if not (re.search('^[ACGT]$',ref) and re.search('^[ACGT]$',alt)): continue 
        if not ((gt == '0|0') or (gt == '0|1') or (gt == '1|0') or (gt == '1|1')): continue
        af=allele_freq.get(pos,0)
        if af < afRefPanelTH or af > 1-afRefPanelTH: continue
        n += 1
        if ((gt == '0|1') or (gt == '1|0')): nHet +=1
        giab_hap[pos]=gt
    print("Number of common genotyped positions in GIAB: {0:d}, of which {1:d} heterozygous".format(n,nHet))

#############################################################################

def read_phased_hap(chrom):
    phased_hap=dict()
    command = f"bcftools view -H -v snps -m2 -M2 {phasedFile} -r {chrom}"
    if confRegionsFlag:
        confRegions=f'../NA12878/VCF/GIAB/confRegions.chr{chrOfInterest}.bed'
        if not os.path.isfile(confRegions): exit(f'File not found: {confRegions}')
        command = f"bcftools view -H -v snps -m2 -M2 {phasedFile} -R {confRegions}"
 
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
        if pos not in giab_hap: continue
        if not ((gt == '0|1') or (gt == '1|0')): continue
        phased_hap[pos]=gt
    print("Number of het positions (present also in GIAB) in phased file:",len(phased_ap))
    return phased_hap

#############################################################################################

def phased_hap_accuracy_by_genomic_size_interval(chrOfInterest):

    gtHapCounts=[]
    for i in range(3000):
        gtHapCounts.append([0,0,0,0,0])

    hapSizef=f"{hapSize//1000000:d}Mb"
    outFile=f'../{sampleID}/SHAPEIT4/PHASING_ANALYSIS/' \
             + f'{sampleID}.{covGL}x.{array}.UKB.{hapSizef}.chr{chrOfInterest}.tsv';
    try:
        fout= open(outFile,'w')
    except:
        exit(f"Can not open: {outFile}")
    for pos in sorted(list(phasedHap.keys())):
        bin=pos//hapSize
        gtHapCounts[bin][0] += 1
        if giabHap[pos] == '0|0':
            gtHapCounts[bin][1] += 1 
        elif giabHap[pos] == '1|1':
            gtHapCounts[bin][2] += 1
        if not (giabHap[pos] == '0|1' or giabHap[pos] == '1|0'): continue
        if giabHap[pos] == phasedHap[pos]:
            gtHapCounts[bin][3] += 1
        else:
            gtHapCounts[bin][4] += 1
        binMax=bin

    centromere_bins=[x//hapSize for x in centromere[chrOfInterest]]
    nBin=nBinCorrect=0
    for bin in range(binMax):
        if not gtHapCounts[bin][0]: continue
        flag='+'
        nCorrect=gtHapCounts[bin][3]
        if (gtHapCounts[bin][4] > gtHapCounts[bin][3]):
            nCorrect=gtHapCounts[bin][4]
            flag='-'
        fracPhased=f"{nCorrect/gtHapCounts[bin][0]:.3f}"
        print(chrOfInterest,int((bin+0.5)*hapSize),*gtHapCounts[bin],fracPhased,flag,sep="\t",file=fout)
        if centromere_bins[0] <= bin <= centromere_bins[1]:continue 
        if gtHapCounts[bin][0] > n_snp_min: 
            nBin+=1
            if float(fracPhased)>=0.985: nBinCorrect +=1
    fout.close()
    return nBin,nBinCorrect

##############################################################################

def phased_hap_accuracy_by_number_of_snps_interval(chrOfInterest):
    max_interval_size=15000000
    n_snp_min=3000
    (centro_start,centro_end)=centromere[chrom]
    snp= sorted(list(phasedHap.keys()))

    n_snps_bin_f=f"{n_snps_bin//1000:d}Ksnp"
    outFile=f'../{sampleID}/SHAPEIT4/PHASING_ANALYSIS/' \
             + f'{sampleID}.{covGL}x.{array}.UKB.{n_snps_bin_f}.chr{chrOfInterest}.tsv';
    try:
        fout= open(outFile,'w')
    except:
        exit(f"Can not open: {outFile}")
    nBin=nBinCorrect=0
    hap_size=[]
    n_snp_tot=len(snp)
    i_start=0
    while (i_start < n_snp_tot):
        pos=snp[i_start]
        while centro_start <= pos <= centro_end:
            i_start += 1
            pos=snp[i_start]
        i_max=i_start+n_snps_bin
        if (i_max > n_snp_tot): i_max=n_snp_tot
        pos_max=snp[-1]
        if snp[i_start] < centro_start:
            pos_max=centro_start
        gtHapCounts=[0,0,0,0,0]
        for i in (range(i_start,i_max)):
            gtHapCounts[0] += 1
            pos=snp[i]
            if giabHap[pos] == '0|0':
                gtHapCounts[1] += 1 
            elif giabHap[pos] == '1|1':
                gtHapCounts[2] += 1
            if not (giabHap[pos] == '0|1' or giabHap[pos] == '1|0'): continue
            if giabHap[pos] == phasedHap[pos]:
                gtHapCounts[3] += 1
            else:
                gtHapCounts[4] += 1
            i_end=i
            if i_end == n_snp_tot-1: break
            pos=snp[i_end+1]
            if pos > pos_max: break
        pos_start=snp[i_start]
        pos_end=snp[i_end]
        int_size=pos_end + 1 - pos_start
        flag='+'
        nCorrect=gtHapCounts[3]
        if (gtHapCounts[4] > gtHapCounts[3]):
            nCorrect=gtHapCounts[4]
            flag='-'
        fracPhased=f"{nCorrect/gtHapCounts[0]:.3f}"
        print(chrOfInterest,pos_start,pos_end,int_size,*gtHapCounts,fracPhased,flag,sep="\t",file=fout)
        if gtHapCounts[0] > n_snp_min: 
            nBin+=1
            hap_size.append(int_size)
            if float(fracPhased)>=0.985: nBinCorrect +=1
        i_start = i_end + 1
#        pos=snp[i_start]
#        while centro_start <= pos <= centro_end:
#            i_start += 1
#            pos=snp[i_start]
    fout.close()
    return nBin,nBinCorrect,hap_size



##############################################################################

centromere_positions()
n_bin_genome=n_bin_correct_genome=0
hap_size_genomewide=[]
for chrom in chrList:
    print(chrom)
    phasedFile=f'../{sampleID}/SHAPEIT4/{sampleID}.{covGL}x.{array}.UKB.chr{chrom}.Output.vcf.gz';
    if not os.path.isfile(phasedFile): exit(f'File not found: {phasedFile}')
    read_allele_freq(chrom)
    read_giab_hap(chrom)
    read_phased_hap(chrom)
#    n_bin_chrom,n_bin_correct_chrom=phased_hap_accuracy_by_genomic_size_interval(chrom)
    n_bin_chrom,n_bin_correct_chrom,hap_size_chrom=phased_hap_accuracy_by_number_of_snps_interval(chrom)
    hap_size_genomewide.extend(hap_size_chrom)
    print(chrom,n_bin_chrom,n_bin_correct_chrom,sep='\t')
    n_bin_genome += n_bin_chrom
    n_bin_correct_genome += n_bin_correct_chrom

frac_bin=f"{n_bin_correct_genome/n_bin_genome:.3f}"
print(n_bin_genome,n_bin_correct_genome,frac_bin)
print(statistics.mean(hap_size_genomewide),statistics.stdev(hap_size_genomewide),
      max(hap_size_genomewide),min(hap_size_genomewide),len(hap_size_genomewide))
print(statistics.quantiles(hap_size_genomewide,n=10))

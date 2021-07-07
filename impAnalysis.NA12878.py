#!/usr/bin/env python3

import gzip,re, os, subprocess
import numpy as np

sampleID='NA12878'
chrOfInterest=8
confRegions=None
#confRegions=f'../NA12878/VCF/GIAB/confRegions.chr{chrOfInterest}.bed'
if confRegions and not os.path.isfile(confRegions): exit('File not found: '+ confRegions)
cov='LowCov'
array='HRC'
gtProbTH=0.99
afTH=0.005

#platypusFile=f'../{sampleID}/VCF/Platypus/{sampleID}.GIAB.30x.chr{chrOfInterest}.{array}.vcf.gz'
#if not os.path.isfile(platypusFile): exit('File not found: '+ platypusFile)

imputedFile=f'../{sampleID}/GLIMPSE/{cov}/{sampleID}.{cov}.chr{chrOfInterest}.bcf'
if not os.path.isfile(imputedFile): exit('File not found: '+ imputedFile)

#####################################################################

def run_command(command):
    p = subprocess.Popen(command,
                         shell=True,
                         text=True,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    return iter(p.stdout.readline, b'')

###############################################################################

def read_alleleFreq():
    global alleleFreq
    alleleFreq=dict()

    vcfFile='../DATA/REF_PANEL/HRC/VCF/sites.vcf.gz'

    n=nCommon=0
    command = f"bcftools query -f'%POS\t%INFO/AF\n' {vcfFile} -r {chrOfInterest}"
    for line in run_command(command):
        if line == '':  break
        line = line.rstrip()
        if line.startswith('#'): continue
        fields = line.split()
        (pos,af)=fields[:]
        af=float(af) 
        alleleFreq[int(pos)]=af 
        if af >= afTH and af <= 1-afTH: nCommon += 1
    n=len(alleleFreq)
    print("Number of genotyped positions in the reference panel: {0:d}, of which {1:d} common".format(n,nCommon))

#############################################################################

def read_giabGT():
    global giabGT

    gtFile=f'../NA12878/VCF/GIAB/snps.{array}.vcf.gz'
    if not os.path.isfile(gtFile): exit('File not found: '+gtFile)

    n=nHet=nHomAlt=0
    nCommon=nCommonHet=nCommonHomAlt=0
    giabGT=dict()
    command = f"bcftools view -H -v snps -m2 -M2 {gtFile} -r {chrOfInterest}"
    if confRegions:
        command = f"bcftools view -H -v snps -m2 -M2 {gtFile} -R {confRegions}"
    for line in run_command(command):
        if line == '':  break
        line = line.rstrip()
        if line.startswith('#'): continue
        fields = line.split()
        if (len(fields) != 10): exit()
        (chrom,pos,id,ref,alt,qual,filter,info,format,gt)=fields[:]
        pos=int(pos)  
        if ((len(ref) > 1) or (len(alt) >1)): continue
        if not (re.search('^[ACGT]$',ref) and re.search('^[ACGT]$',alt)): continue 
        gt=gt.replace("|","/")
        if (gt == '1/0'): gt='0/1'
        if not ((gt == '0/0') or (gt == '0/1') or (gt == '1/1')): continue
        n += 1
        if (gt == '0/1'): 
            nHet +=1
        elif (gt == '1/1'):
            nHomAlt +=1
        af=alleleFreq.get(pos,0)
        if af < afTH or af > 1-afTH: continue
        nCommon += 1
        if (gt == '0/1'): 
            nCommonHet +=1
        elif (gt == '1/1'):
            nCommonHomAlt +=1
#        matchObj=re.search('^AF=(.*)',info) 
#        af=float(matchObj.group(1))
        giabGT[pos]=[ref,alt,gt,af]

    print("Number of genotyped positions in GIAB: {0:d}, of which {1:d} heterozygous and {2:d} homozygous alt".format(n,nHet,nHomAlt))
    print("Number of common genotyped positions in GIAB: {0:d}, of which {1:d} heterozygous and {2:d} homozygous alt".
          format(nCommon,nCommonHet,nCommonHomAlt))


#############################################################################

def read_platypusGT():
    global platypusGT
    gqMin=20


    nHet=0;nHomAlt=0;nHomRef=0;
    platypusGT=dict()
    command = f"bcftools view -H -v snps -m2 -M2 {platypusFile} -r {chrOfInterest}"
    for line in run_command(command):
        if line == '':  break
        line = line.rstrip()
        if line.startswith('#'): continue
        fields = line.split()
        if (len(fields) != 10): exit()
        (chr,pos,id,ref,alt,qual,filter,info,format,sample)=fields[:]
        chr=chr.strip('^chr')
        if (len(ref) > 1) or (len(alt) >1): continue
        if not re.search('^[ACGT]$',ref) and re.search('^[ACGT]$',alt): continue
        (gt,gl,gof,gq,nr,nv)=sample.split(':')
        gq=int(gq)
        if (gq < gqMin): continue
        gt=gt.replace("|","/")
        if (gt == '1/0'): gt='0/1'
        if not ((gt == '0/0') or (gt == '0/1') or (gt == '1/1')): continue
        if ((gt == '0/1') or (gt == '1/1')):
            if (filter != 'PASS'): continue
        platypusGT[int(pos)]=[ref,alt,gt,gq]    
        if (gt == '0/0'):
            nHomRef +=1
        elif (gt == '0/1'):
            nHet +=1
        elif (gt == '1/1'):
            nHomAlt +=1
        
    print("Number of genotyped positions in platypus file: %d, of which %d heterozygous and %d homozygous alt" % (len(platypusGT),nHet,nHomAlt))

#############################################################################

def read_bcftoolsGT():
    global bcftoolsGT

    #    bcftoolsFile='../../Tools/GLIMPSE/tutorial_hg19/GLIMPSE_validation/NA12878.chr22.vcf.gz'
    #    bcftoolsFile='../../Tools/GLIMPSE/tutorial/GLIMPSE_validation/EUR.validation.NA12878.chr22.bcf'
    bcftoolsFile='~/Tools/GLIMPSE/tutorial/GLIMPSE_validation/EUR.mock.vcf.gz'

    nHet=0;nHomAlt=0;nHomRef=0;
    bcftoolsGT=dict()
    command = f"bcftools view -H -v snps -m2 -M2 {bcftoolsFile} -r {chrOfInterest}"
    for line in run_command(command):
        if line == '':  break
        line = line.rstrip()
        if line.startswith('#'): continue
        fields = line.split()
        if (len(fields) != 10): exit()
        (chr,pos,id,ref,alt,qual,filter,info,format,sample)=fields[:]
        (gt,pl,dp)=sample.split(':')
        if (gt == '.'): continue
        plPhred=pl.split(',')
        for x in plPhred:
            if x =='.':
                print(line)
                exit()
        plProb=[10**(-float(x)/10) for x in plPhred]
        probTot=sum(plProb)
        plProb=[x/probTot for x in plProb]
        if ((len(ref) > 1) or (len(alt) >1)): continue
        if not re.search('^[ACGT]$',ref) and re.search('^[ACGT]$',alt): continue
        if (gt == '1/0'): gt='0/1'
        if not ((gt == '0/0') or (gt == '0/1') or (gt == '1/1')): continue
        if (gt == '0/0'):
            iP=0
            nHomRef +=1
        elif (gt == '0/1'):
            iP=1 
            nHet +=1
        elif (gt == '1/1'):
            iP=2
            nHomAlt +=1
        bcftoolsGT[int(pos)]=[ref,alt,gt,int(dp)]+ plProb

    print("Number of genotyped positions in bcftools file: %d, of which %d heterozygous and %d homozygous alt" % (len(bcftoolsGT),nHet,nHomAlt))


##############################################################################

def read_imputedGT():
    global imputedGT

    n=nHet=nHomAlt=0
    nFilter=nFilterHet=nFilterHomAlt=0
    imputedGT=dict()
    command = f"bcftools view -H -v snps -m2 -M2 {imputedFile} -r {chrOfInterest}"
    if confRegions:
        command = f"bcftools view -H -v snps -m2 -M2 {imputedFile} -R {confRegions}"
    for line in run_command(command):
        if line == '':  break                       
        line = line.rstrip()
        if line.startswith('#'): continue
        fields = line.split()
        if (len(fields) != 10): exit()
        (chr,pos,id,ref,alt,qual,filter,info,format,sample)=fields[:]
        pos=int(pos)
        if ((len(ref) > 1) or (len(alt) >1)): continue
        if not re.search('^[ACGT]$',ref) and re.search('^[ACGT]$',alt): continue
        (gt,ds,gp,hs)=sample.split(':')
        if (gt == '1/0'): gt='0/1'
        if not ((gt == '0/0') or (gt == '0/1') or (gt == '1/1')): continue
        iP=0 if (gt == '0/0') else 1 if (gt == '0/1') else 2
        gtProb=gp.split(',')[iP]
        gtProb=float(gtProb)
#        raf=None
#        for infoField in info.split(';'):
#            matchObj=re.search('^RAF=(.*)',infoField)
#            if matchObj: raf=matchObj.group(1)
#            break
#        if raf == None : exit('Raf not found'+line+'\n')
        n += 1
        if (gt == '0/1'): 
            nHet +=1
        elif (gt == '1/1'):
            nHomAlt +=1
        af=alleleFreq.get(pos,0)
        if af < afTH or af > 1-afTH: continue
        if gtProb < gtProbTH: continue
        nFilter += 1
        if (gt == '0/1'): 
            nFilterHet +=1
        elif (gt == '1/1'):
            nFilterHomAlt +=1
        imputedGT[pos]=[ref,alt,gt,gtProb,af]

    print("Number of genotyped positions in imputed file: {0:d}, of which {1:d} heterozygous and {2:d} homozygous alt".
          format(n,nHet,nHomAlt))
    print("Number of filtered genotyped positions in imputed file: {0:d}, of which {1:d} heterozygous and {2:d} homozygous alt".
          format(nFilter,nFilterHet,nFilterHomAlt))

###############################################################################

def gt_concordance():
    
    
    ds_bcftools=[]
    ds_imputed=[]
    afMin,afMax=(0.3,0.41)  # 0.00000 0.00100 0.00200 0.00500 0.01000 0.05000 0.10000 0.20000 0.50000
    for pos in bcftoolsGT:
        pos=int(pos)
        dp=bcftoolsGT[pos][3]
        
        if (dp < 8): continue
        if  pos not in imputedGT: continue
        if (imputedGT[pos][3] < gtProbTH): continue
        af=alleleFreq.get(pos,0)
        swap=(af>0.5)
        if swap: 
            af=1-af
        if not (afMin <= af < afMax): continue 
        bcftoolsDS=bcftoolsGT[pos][5]+2*bcftoolsGT[pos][6]
#        imputedDS=float(imputedGT[pos][7])+2*float(imputedGT[pos][6])
        imputedDS=imputedGT[pos][4]
        if swap:
            bcftoolsDS=bcftoolsGT[pos][5]+2*bcftoolsGT[pos][4]
            imputedDS=float(imputedGT[pos][7])+2*float(imputedGT[pos][6])

#        if bcftoolsGT[pos][2] == '0/0':
#            if (af < 0.5):
#                iBCF = 0
#            else:
#                iBCF = 0
#        elif bcftoolsGT[pos][2] == '0/1':
#            iBCF = 1
#        elif bcftoolsGT[pos][2] == '1/1':
#            if (af < 0.5):
#                iBCF = 2
#            else:
#                iBCF = 2 
#        else: 
#            exit("Something wrong with giabGT genotype at "+pos+" :"+bcftoolsGT[pos])
        if imputedGT[pos][2] == '0/0':
            if (af < 0.5):
                iImp = 0
            else:
                iImp = 0
        elif  imputedGT[pos][2] == '0/1':
            iImp = 1
        elif  imputedGT[pos][2] == '1/1':
            if (af < 0.5):
                iImp = 2
            else:
                iImp = 2
        else:
            exit("Something wrong with imputedGT genotype at "+pos+" : "+imputedGT[pos])
        ds_bcftools.append(bcftoolsDS)
        ds_imputed.append(imputedDS)
    
    print(len(ds_imputed)) 
    r = np.corrcoef(ds_bcftools,ds_imputed)
    print(r[0][1]**2)    
         
##############################################################################

def gtAccuracy():
    concordMatrix=[[0,0,0],[0,0,0],[0,0,0]]
    notInGIAB=notInGIABHet=0
    for pos in imputedGT:
        if pos not in giabGT:
            notInGIAB += 1
            if imputedGT[pos][2] == '0/1': notInGIABHet += 1
            continue
        if imputedGT[pos][2] == '0/0':
            iImp = 0
        elif  imputedGT[pos][2] == '0/1':
            iImp = 1
        elif  imputedGT[pos][2] == '1/1':
            iImp = 2
        else:
            exit("Something wrong with imputedGT genotype at "+pos+" : "+imputedGT[pos])
        if giabGT[pos][2] == '0/0':
            iGIAB = 0
        elif  giabGT[pos][2] == '0/1':
            iGIAB = 1
        elif  giabGT[pos][2] == '1/1':
            iGIAB = 2
        else: 
            exit("Something wrong with giabGT genotype at "+pos+" :"+giabGT[pos])
        concordMatrix[iImp][iGIAB] += 1

    tot=acc=0
    for i in range(3):
        for j in range(3):
            tot += concordMatrix[i][j]
            if i==j: acc += concordMatrix[i][j]
        print(concordMatrix[i])

    acc /= tot
    print("Overall accuracy of {0:.5f} over {1:d} overlapping positions".format(acc,tot))
    print("Number of imputed calls not in GIAB: {0:d}, of which {1:d} heterozygous".format(notInGIAB,notInGIABHet))

    notImputed=notImputedHet=0
    for pos in giabGT:
        if (pos in imputedGT): continue
        notImputed +=1
        if (giabGT[pos][2] == '0/1'): notImputedHet +=1
    print("Number of GIAB calls not imputed: {0:d}, of which {1:d} heterozygous".format(notImputed,notImputedHet))


###############################################################################

def imputationAccuracy():

    nNotImputed=nWrongAlt=0
    concordMatrix=[[0,0,0],[0,0,0],[0,0,0]]
    dsGIAB=[];
    dsImputed=[] 
    for pos in giabGT:
        if int(pos) not in imputedGT: 
            nNotImputed +=1
            continue
        if (imputedGT[pos][3] < gtProbTH):
            del imputedGT[pos] 
            continue
        if (giabGT[pos][1] != imputedGT[pos][1]):
            nWrongAlt += 1
            continue
        if giabGT[pos][2] == '0/0':
            iGIAB = 0
        elif  giabGT[pos][2] == '0/1':
            iGIAB = 1
        elif  giabGT[pos][2] == '1/1':
            iGIAB = 2
        else: 
            exit("Something wrong with giabGT genotype at "+pos+" :"+giabGT[pos])
        if imputedGT[pos][2] == '0/0':
            iImp = 0
        elif  imputedGT[pos][2] == '0/1':
            iImp = 1
        elif  imputedGT[pos][2] == '1/1':
            iImp = 2
        else:
            exit("Something wrong with imputedGT genotype at "+pos+" : "+imputedGT[pos])
        concordMatrix[iGIAB][iImp] += 1
        dsGIAB.append(iGIAB)
        dsImputed.append(imputedGT[pos][4])
    
    r = np.corrcoef(dsGIAB,dsImputed)
    print(r[0][1]**2)    
    print(nWrongAlt)

    tot=acc=totNoHomRef=accNoHomRef=0
    for iGIAB in range(3):
        for iImp in range(3):
            tot += concordMatrix[iGIAB][iImp]
            if iGIAB==iImp: acc += concordMatrix[iGIAB][iImp]
            if ((iGIAB != 0) or (iImp != 0)):
               totNoHomRef += concordMatrix[iGIAB][iImp]
               if iGIAB==iImp: accNoHomRef += concordMatrix[iGIAB][iImp]
        print(concordMatrix[iGIAB])

    acc /= tot
    accNoHomRef /= totNoHomRef
    print("Total number of imputed file: %d (excluding true positive homozygous ref sites: %d) " % (tot,totNoHomRef))
    print("Overall acuracy: %f (excluding true positive homozygous ref sites: %f)" % (acc, accNoHomRef))
    print('Number of GIAB SNP not imputed',nNotImputed)


#############################################################################

def outputVCF():
    vcfFile=f'/scratch/DMP/BIOCHLEU/slise/LowPassWGS/{sampleID}/SHAPEIT4/{sampleID}.{cov}.chr{chrOfInterest}.Input.vcf'
    fout=open(vcfFile,'w')
    fout.write('##fileformat=VCFv4.1\n')
    fout.write('##FILTER=<ID=PASS,Description="All filters passed">\n')
    fout.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    fout.write('##FORMAT=<ID=GQ,Number=.,Type=Integer,Description="Genotype quality as phred score">\n')
    fout.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'+sampleID+'\n')
    for pos in sorted(list(imputedGT.keys())):
            (ref,alt,gt,gtProb)=imputedGT[pos][0:4]
            fout.write("%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (chrOfInterest,pos,'.',ref,alt,'.','.','.','GT',gt))
    fout.close()
    subprocess.call(["bgzip","-f",vcfFile])
    subprocess.call(["tabix","-fp","vcf",vcfFile+'.gz'])


##############################################################################

read_alleleFreq()

read_giabGT()
#read_bcftoolsGT()
read_imputedGT()
gtAccuracy()
#read_platypusGT()
#read_alleleFreq()
#gt_concordance()
#imputationAccuracy()
#outputVCF()

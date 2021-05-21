#!/bin/bash
#
#SBATCH --job-name=gtPlatypus
#SBATCH --account=DMPGXAAAM
#SBATCH --mail-user=stefano.lise@icr.ac.uk
#SBATCH --mail-type=END
#SBATCH --output=LOG/gtPlatypus.%A.out
#SBATCH --error=LOG/gtPlatypus.%A.err
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=2G

module load anaconda/3
source /opt/software/applications/anaconda/3/etc/profile.d/conda.sh
conda init
conda activate platypus-variant0.8.1.2

patID=TR064;
sampleID=${patID}_BL
array=HRC;
refGenome=/data/scratch/DMP/UCEC/UROTRBIO/slise/DATA/GENOMES/HS37D5/hs37d5.fa;        # hs37d5.fa or ucsc.hg19.fasta

bamDir=/data/scratch/DMP/UCEC/UROTRBIO/slise/ROSETREES/${patID}/BAM
bamFiles=${bamDir}/${sampleID}.bam
if [[ $array == "HRC" ]]; then
    sourceFile=/data/scratch/DMP/UCEC/UROTRBIO/slise/ROSETREES/DATA/SNPs/${array}/snps.f005.vcf.gz;
elif [[ $array == "UKB" ]]; then
    sourceFile=/data/scratch/DMP/UCEC/UROTRBIO/slise/ROSETREES/DATA/SNPs/${array}/snps.vcf.gz;
fi
outFile=/data/scratch/DMP/UCEC/UROTRBIO/slise/ROSETREES/${patID}/VCF/${sampleID}.${array}.vcf ;

platypus callVariants --refFile=$refGenome --bamFiles=${bamFiles} --output=$outFile --source=$sourceFile --minPosterior=0  --getVariantsFromBAMs=0 --logFileName=$outFile.log --verbosity=1 --trimOverlapping=0
bgzip -f $outFile
tabix -fp vcf $outFile.gz

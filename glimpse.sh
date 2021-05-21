#!/bin/bash
#
#SBATCH --job-name=glimpse
#SBATCH --account=DMPGXAAAM
#SBATCH --mail-user=stefano.lise@icr.ac.uk
#SBATCH --mail-type=END
#SBATCH --output=LOG/glimpse.%A.%a.out
#SBATCH --error=LOG/glimpse.%A.%a.err
#SBATCH --array=1-22
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=3G

PAT_ID=TR064

#
# Skeleton of a pipeline based on the GLIMPSE tutorial 
#            https://odelaneau.github.io/GLIMPSE/tutorial_hg19.html   

module load SAMtools BCFtools HTSlib
REFGEN=/data/scratch/DMP/UCEC/UROTRBIO/slise/DATA/GENOMES/HS37D5/hs37d5.fa
COV=4x
CHR=$SLURM_ARRAY_TASK_ID

# Split bam file by chromosomes

IN_BAM=/data/scratch/DMP/UCEC/UROTRBIO/slise/ROSETREES/${PAT_ID}/BAM/${PAT_ID}_GL.${COV}.bam
OUT_BAM=/data/scratch/DMP/UCEC/UROTRBIO/slise/ROSETREES/${PAT_ID}/GLIMPSE/BAM/${PAT_ID}_GL.${COV}.chr${CHR}.bam
samtools view -bh $IN_BAM -o $OUT_BAM $CHR
samtools index $OUT_BAM

# step 3.2: computing GLs

#
SITES_VCF=/data/scratch/DMP/UCEC/UROTRBIO/slise/ROSETREES/DATA/SNPs/HRC/snps.chr${CHR}.vcf.gz
SITES_TSV=/data/scratch/DMP/UCEC/UROTRBIO/slise/ROSETREES/DATA/SNPs/HRC/snps.chr${CHR}.tsv.gz
OUT_GL=/data/scratch/DMP/UCEC/UROTRBIO/slise/ROSETREES/${PAT_ID}/GLIMPSE/GL/${PAT_ID}.${COV}.chr${CHR}.gl.vcf.gz
BAM=/data/scratch/DMP/UCEC/UROTRBIO/slise/ROSETREES/${PAT_ID}/GLIMPSE/BAM/${PAT_ID}_GL.${COV}.chr${CHR}.bam

#############################################################

# Commands to create the SITES_VCF and SITES_TSV files if not already presentx

#bcftools view -O z -r $CHR -o $SITES_VCF /data/scratch/DMP/UCEC/UROTRBIO/slise/ROSETREES/DATA/SNPs/HRC/snps.vcf.gz
#tabix -fp vcf $SITES_VCF
#bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' $SITES_VCF | bgzip -c > $SITES_TSV
#tabix -s1 -b2 -e2 $SITES_TSV


###################################################################################

bcftools mpileup -f ${REFGEN} -I -E -a 'FORMAT/DP' -T ${SITES_VCF} -r $CHR ${BAM} -Ou | bcftools call -Aim -C alleles -T ${SITES_TSV} -Oz -o ${OUT_GL}
bcftools index -f ${OUT_GL}


# step 4:

GLIMPSE_chunk=~/TOOLS/GLIMPSE/static_bins/GLIMPSE_chunk_static
CHUNKS=/data/scratch/DMP/UCEC/UROTRBIO/slise/ROSETREES/${PAT_ID}/GLIMPSE/CHUNKS/chunks.chr${CHR}.txt

$GLIMPSE_chunk --input $SITES_VCF --region $CHR --window-size 2000000 --buffer-size 200000 --output $CHUNKS

# step 5

GLIMPSE_phase=~/TOOLS/GLIMPSE/static_bins/GLIMPSE_phase_static
IN_GL=/data/scratch/DMP/UCEC/UROTRBIO/slise/ROSETREES/${PAT_ID}/GLIMPSE/GL/${PAT_ID}.${COV}.chr${CHR}.gl.vcf.gz
REF=/data/scratch/DMP/UCEC/UROTRBIO/slise/ROSETREES/DATA/REF_PANEL/HRC/VCF/chr${CHR}.HRC.vcf.gz
MAP=~/TOOLS/GLIMPSE/maps/genetic_maps.b37/chr${CHR}.b37.gmap.gz

while IFS="" read -r LINE || [ -n "$LINE" ];
do
   printf -v ID "%02d" $(echo $LINE | cut -d" " -f1)
   IRG=$(echo $LINE | cut -d" " -f3)
   ORG=$(echo $LINE | cut -d" " -f4)
   OUT=/data/scratch/DMP/UCEC/UROTRBIO/slise/ROSETREES/${PAT_ID}/GLIMPSE/CHUNKS/${PAT_ID}.${COV}.chr${CHR}.imputed.${ID}.bcf
   $GLIMPSE_phase --input ${IN_GL} --reference ${REF} --map ${MAP} --input-region ${IRG} --output-region ${ORG} --output ${OUT}
   bcftools index -f ${OUT}
done < $CHUNKS

#exit 0

# step 6

GLIMPSE_ligate=~/TOOLS/GLIMPSE/static_bins/GLIMPSE_ligate_static
LST=/data/scratch/DMP/UCEC/UROTRBIO/slise/ROSETREES/${PAT_ID}/GLIMPSE/CHUNKS/list.${COV}.chr${CHR}.txt
ls /data/scratch/DMP/UCEC/UROTRBIO/slise/ROSETREES/${PAT_ID}/GLIMPSE/CHUNKS/${PAT_ID}.${COV}.chr${CHR}.imputed.*.bcf > ${LST}
OUT=/data/scratch/DMP/UCEC/UROTRBIO/slise/ROSETREES/${PAT_ID}/GLIMPSE/${PAT_ID}.${COV}.chr${CHR}.bcf

$GLIMPSE_ligate --input ${LST} --output $OUT
bcftools index -f ${OUT}

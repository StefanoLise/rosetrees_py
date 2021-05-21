#!/bin/bash
#
#SBATCH --job-name=replace_read_groups
#SBATCH --account=DMPGXAAAM
#SBATCH --mail-user=stefano.lise@icr.ac.uk
#SBATCH --mail-type=END
#SBATCH --output=LOG/replace_read_groups.%A.out
#SBATCH --error=LOG/replace_read_groups.%A.err
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=1G

module load picard-tools/2.23.8 SAMtools
pat_id=TR081
in_bam=/data/scratch/DMP/UCEC/UROTRBIO/slise/ROSETREES/${pat_id}/BAM/H860003.deDup.bam
out_bam=/data/scratch/DMP/UCEC/UROTRBIO/slise/ROSETREES/${pat_id}/BAM/${pat_id}_PD2.bam
picard AddOrReplaceReadGroups I=$in_bam O=$out_bam  RGID=$pat_id RGLB=$pat_id RGPL=ILLUMINA RGPU=NA RGSM=$pat_id  
samtools index $out_bam

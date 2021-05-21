#!/bin/bash
#
#SBATCH --job-name=shapeit4
#SBATCH --account=DMPGXAAAM
#SBATCH --mail-user=stefano.lise@icr.ac.uk
#SBATCH --mail-type=END
#SBATCH --output=LOG/shapeit4.%A.%a.out
#SBATCH --error=LOG/shapeit4.%A.%a.err
#SBATCH --array=6-6 
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=6
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=8042

# use -n 4 with UKB and -n 6 with HRC

module load anaconda/3
source /opt/software/applications/anaconda/3/etc/profile.d/conda.sh
conda init
conda activate shapeit4.1.3
module load BCFtools

patID=TR019
fileID=${patID}.4x
chr=$SLURM_ARRAY_TASK_ID
mapDir=/data/scratch/DMP/UCEC/UROTRBIO/slise/ROSETREES/DATA/REF_PANEL/MAPS
mapFile=${mapDir}/chr${chr}.b37.gmap.gz ;    

array=UKB;
refDir=/data/scratch/DMP/UCEC/UROTRBIO/slise/ROSETREES/DATA/REF_PANEL/${array}/VCF
refFile=${refDir}/chr${chr}.${array}.vcf.gz
workDir=/data/scratch/DMP/UCEC/UROTRBIO/slise/ROSETREES/${patID}/SHAPEIT4
#inFile=${workDir}/${fileID}.${array}.Input.vcf.gz
inFile=${workDir}/${fileID}.HRC.Input.vcf.gz
outFile=${workDir}/${fileID}.${array}.chr${chr}.Output ;  # Output from Shapeit4. Include both SNP and Indels
#shapeit4 --input $inFile --map $mapFile --region $chr --reference $refFile --output ${outFile}.vcf.gz --thread 4 --log ${outFile}.log
#tabix -fp vcf ${outFile}.vcf.gz

array=HRC
refDir=/data/scratch/DMP/UCEC/UROTRBIO/slise/ROSETREES/DATA/REF_PANEL/${array}/VCF
refFile=${refDir}/chr${chr}.${array}.vcf.gz
workDir=/data/scratch/DMP/UCEC/UROTRBIO/slise/ROSETREES/${patID}/SHAPEIT4
inFile=${workDir}/${fileID}.${array}.Input.vcf.gz
scaffoldFile=${workDir}/${fileID}.UKB.chr${chr}.Output.vcf.gz
outFile=${workDir}/${fileID}.${array}.UKB.chr${chr}.Output ;     # Output from Shapeit4. Include both SNP and Indels
shapeit4 --input $inFile --map $mapFile --region $chr --reference $refFile --output ${outFile}.vcf.gz --thread 6 --log ${outFile}.log --scaffold ${scaffoldFile}
tabix -fp vcf ${outFile}.vcf.gz



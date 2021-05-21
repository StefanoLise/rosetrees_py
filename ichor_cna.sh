#!/bin/bash
#



module load anaconda/3
source /opt/software/applications/anaconda/3/etc/profile.d/conda.sh
conda init


# https://github.com/broadinstitute/ichorCNA/wiki
scratch_dir=/data/scratch/DMP/UCEC/UROTRBIO/slise/ROSETREES

pat_id=TR064
gl_sample=${pat_id}_GL.4x
bl_sample=${pat_id}_PD2
out_dir=${scratch_dir}/${pat_id}/ICHOR_CNA

gl_bam=${scratch_dir}/${pat_id}/BAM/${gl_sample}.bam
bl_bam=${scratch_dir}/${pat_id}/BAM/${bl_sample}.bam
window=1000000
if [[ $window -eq 10000 ]]; then
	window_f=10kb
elif [[ $window -eq 50000 ]]; then
	window_f=50kb
elif [[ $window -eq 500000 ]]; then
	window_f=500kb
elif [[ $window -eq 1000000 ]]; then
	window_f=1000kb	
else
	echo "Window size not acceptable"
	exit 0
fi
gl_wig=${out_dir}/${gl_sample}.${window_f}.wig
bl_wig=${out_dir}/${bl_sample}.${window_f}.wig


min_qual=20
region="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y"

conda activate hmmcopy0.1.1
if [ ! -f $gl_wig ]; then
	echo "Running readCounter on ${gl_sample}"
	readCounter --window $window --quality $min_qual --chromosome $region $gl_bam > $gl_wig
fi
if [ ! -f $bl_wig ]; then
	echo "Running readCounter on ${bl_sample}"
	readCounter --window $window --quality $min_qual --chromosome $region $bl_bam > $bl_wig
fi
conda deactivate

conda activate r-essentials4.0
runIchorCNA="Rscript /home/slise/TOOLS/ichorCNA/scripts/runIchorCNA.R"
gc_wig=/home/slise/TOOLS/ichorCNA/inst/extdata/gc_hg19_${window_f}.wig
map_wig=/home/slise/TOOLS/ichorCNA/inst/extdata/map_hg19_${window_f}.wig
id=${bl_sample}
centromere=/home/slise/TOOLS/ichorCNA/inst/extdata/GRCh37.p13_centromere_UCSC-gapTable.txt
normal=0.9                                # fraction of normal cells           
ploidy="c(2)"
max_cn=4
sc_states="c()"

#$runIchorCNA --id $id --WIG $bl_wig --NORMWIG $gl_wig \
$runIchorCNA --id $id --WIG $bl_wig \
  --ploidy $ploidy --normal $normal --maxCN $max_cn \
  --gcWig $gc_wig --mapWig $map_wig --centromere $centromere \
  --includeHOMD False --chrs "c(1:22)" --chrTrain "c(1:22)" \
  --estimateNormal True --estimatePloidy True --estimateScPrevalence True \
  --scStates $sc_states --txnE 0.9999 --txnStrength 10000 --outDir $out_dir 
conda deactivate


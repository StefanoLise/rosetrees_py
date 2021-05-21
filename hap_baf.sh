#!/bin/bash
#
#SBATCH --job-name=hap_baf
#SBATCH --account=DMPGXAAAM
#SBATCH --mail-user=stefano.lise@icr.ac.uk
#SBATCH --mail-type=END
#SBATCH --output=LOG/hap_baf.%A.out
#SBATCH --error=LOG/hap_baf.%A.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=12:00:00

module load SAMtools BCFtools HTSlib
pat_id=TR019
python hap_baf.py $pat_id 100
python hap_baf.py $pat_id 200
python hap_baf.py $pat_id 400
python hap_baf.py $pat_id 1000
python hap_baf.py $pat_id 2000
python hap_baf.py $pat_id 4000




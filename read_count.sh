#!/bin/bash
#
#SBATCH --job-name=read_count
#SBATCH --account=DMPGXAAAM
#SBATCH --mail-user=stefano.lise@icr.ac.uk
#SBATCH --mail-type=END
#SBATCH --output=LOG/read_count.%A.out
#SBATCH --error=LOG/read_count.%A.err
#SBATCH --ntasks=1
#SBATCH --time=12:00:00
#SBATCH --mem-per-cpu=8042

module load SAMtools BCFtools HTSlib
python read_count.py


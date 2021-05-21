#!/bin/bash
#
#SBATCH --job-name=make_bl_sample
#SBATCH --account=DMPGXAAAM
#SBATCH --mail-user=stefano.lise@icr.ac.uk
#SBATCH --mail-type=END
#SBATCH --output=LOG/make_bl_sample.%A.out
#SBATCH --error=LOG/make_bl_sample.%A.err
#SBATCH --ntasks=1
#SBATCH --time=12:00:00
#SBATCH --mem-per-cpu=2G

module load SAMtools
python make_bl_sample.py


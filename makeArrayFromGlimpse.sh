#!/bin/bash
#
#SBATCH --job-name=makeArrayFromPlatypus
#SBATCH --account=DMPGXAAAM
#SBATCH --mail-user=stefano.lise@icr.ac.uk
#SBATCH --mail-type=END
#SBATCH --output=LOG/makeArrayFromGlimpse.%A.out
#SBATCH --error=LOG/makeArrayFromGlimpse.%A.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=7G
#SBATCH --time=12:00:00

module load SAMtools BCFtools HTSlib
python makeArrayFromGlimpse.py


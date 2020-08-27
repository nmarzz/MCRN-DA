#!/bin/bash
#SBATCH -p batch
#SBATCH -N 1
#SBATCH --time=24:00:00
#SBATCH --mem=2GB

module load matlab/2018b

matlab nodisplay -r "PFdriveswe" > output_process_SWE.log
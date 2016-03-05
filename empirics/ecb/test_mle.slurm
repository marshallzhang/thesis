#!/bin/bash
#SBATCH -n 1               # (Max) number of tasks per job, for R usually 1
#SBATCH -N 1
#SBATCH -o dump/out_mle.txt         # File for the standard output
#SBATCH -e dump/err_mle_%a.txt         # File for the standard error
#SBATCH --open-mode=append # Append to standard output and error files
#SBATCH -p general  # Partition to use
#SBATCH --mem-per-cpu=4096  # Memory required per CPU, in MegaBytes
#SBATCH -t 1:00:00        # Running time of 1 min
#SBATCH -a 1-68

Rscript for_odyssey.R
#!/bin/bash
#SBATCH -n 1               # (Max) number of tasks per job, for R usually 1
#SBATCH -N 1
#SBATCH -o dump/out_%A_%a.txt         # File for the standard output
#SBATCH -e dump/err_%A_%a.txt         # File for the standard error
#SBATCH --open-mode=append # Append to standard output and error files
#SBATCH -p serial_requeue  # Partition to use
#SBATCH --mem-per-cpu=4096  # Memory required per CPU, in MegaBytes
#SBATCH -t 12:00:00        # Running time of 1 min
#SBATCH -a 1-9

Rscript for_odyssey.R
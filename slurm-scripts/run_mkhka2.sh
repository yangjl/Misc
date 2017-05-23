#!/bin/bash
#SBATCH -D /home/jolyang/Documents/Github/Misc
#SBATCH -o /home/jolyang/Documents/Github/Misc/slurm-log/testout-%j.txt
#SBATCH -e /home/jolyang/Documents/Github/Misc/slurm-log/err-%j.txt
#SBATCH -J hka2
#SBATCH --mail-user=yangjl0930@gmail.com
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL #email if fails
set -e
set -u

MLHKA -i largedata/infile2.txt -o largedata/outfile2.txt -s 1234 -c 100000

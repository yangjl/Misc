#!/bin/bash
#SBATCH -D /home/jolyang/Documents/Github/Misc
#SBATCH -o /home/jolyang/Documents/Github/Misc/slurm-log/testout-%j.txt
#SBATCH -e /home/jolyang/Documents/Github/Misc/slurm-log/error-%j.txt
#SBATCH -J pheno-fit
set -e
set -u

R --no-save < profiling/7.pheno/7.1.A_pheno_KRN-GWAS.R

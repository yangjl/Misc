#!/bin/bash
#SBATCH -D /Users/yangjl/Documents/Github/pvpDiallel
#SBATCH -o /home/jolyang/Documents/pvpDiallel/slurm-log/testout-%j.txt
#SBATCH -e /home/jolyang/Documents/pvpDiallel/slurm-log/error-%j.txt
#SBATCH -J tw
set -e
set -u

R --no-save < profiling/7.pheno/7.1.A_pheno_KRN-GWAS.R

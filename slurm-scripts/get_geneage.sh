#!/bin/bash
#SBATCH -D /home/jolyang/dbcenter/phyloBlastDB
#SBATCH -o /home/jolyang/Documents/Github/Misc/slurm-log/testout-%j.txt
#SBATCH -e /home/jolyang/Documents/Github/Misc/slurm-log/error-%j.txt
#SBATCH -J geneage
set -e
set -u

PSmap -i ~/Documents/Github/Misc/largedata/2.geneage/Zea_mays.AGPv3.25_newheader.fa -d phyloBlastDB.fa -p BLAST_AGPv3 -r maize_blast_results -t 30 -a 16

python /home/jolyang/bin/send_email.py -s slurm-scripts/get_geneage.sh

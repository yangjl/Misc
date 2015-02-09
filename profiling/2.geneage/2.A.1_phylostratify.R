### Jinliang Yang


### download maize protein
###  genome build AGPv3.25

wget ftp://ftp.ensemblgenomes.org/pub/plants/release-25/fasta/zea_mays/pep/Zea_mays.AGPv3.25.pep.all.fa.gz
gzip -d Zea_mays.AGPv3.25.pep.all.fa.gz 
rep -o ">" Zea_mays.AGPv3.25.pep.all.fa | wc -w


#####
cd ~/Documents/Github/Misc/largedata/2.geneage/
PSmap -i ~/Documents/Github/Misc/largedata/2.geneage/mays_test.fa -d phyloBlastDB.fa -p BLAST_maize_test -r maize_blast_results -t 30 -a 1


### substitute in vi
%s/pep.*/| [Zea mays] | [ Eukaryota; Viridiplantae; Streptophyta; Streptophytina; Embryophyta; Tracheophyta; Euphyllophyta; Spermatophyta; Magnoliophyta; Liliopsida; commelinids;
Poales; Poaceae; PACMAD clade; Panicoideae; Andropogoneae; Zea ]/g


### sut up slurm run
source("lib/setUpslurm.R")
setUpslurm(slurmsh="slurm-scripts/get_geneage.sh",
           oneline=TRUE,
           codesh="PSmap -i ~/Documents/Github/Misc/largedata/2.geneage/Zea_mays.AGPv3.25_newheader.fa -d phyloBlastDB.fa -p BLAST_AGPv3 -r maize_blast_results -t 30 -a 16",
           wd="/home/jolyang/dbcenter/phyloBlastDB",
           sbatho="/home/jolyang/Documents/Github/Misc/slurm-log/testout-%j.txt",
           sbathe="/home/jolyang/Documents/Github/Misc/slurm-log/error-%j.txt",
           sbathJ="geneage")

###>>> In this path: cd /home/jolyang/dbcenter/phyloBlastDB
###>>> note --ntask=x, 8GB of memory per CPU
###>>> RUN: sbatch -p serial --ntasks 16 slurm-scripts/get_geneage.sh

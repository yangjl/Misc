### Jinliang Yang


### download maize protein
###  genome build AGPv3.25

wget ftp://ftp.ensemblgenomes.org/pub/plants/release-25/fasta/zea_mays/pep/Zea_mays.AGPv3.25.pep.all.fa.gz
gzip -d Zea_mays.AGPv3.25.pep.all.fa.gz 
rep -o ">" Zea_mays.AGPv3.25.pep.all.fa | wc -w


#####
cd ~/Documents/Github/Misc/largedata/2.geneage/
PSmap -i ~/Documents/Github/Misc/largedata/2.geneage/mays_test.fa -d phyloBlastDB.fa -p BLAST_maize_test -r maize_blast_results -t 30 -a 1







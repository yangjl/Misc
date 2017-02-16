### Jinliang
### Feb 16th, 2017

#source("~/Documents/Github/zmSNPtools/Rcodes/dsnp2GenABEL.R")
#library("data.table")
library("GenABEL.data", lib="~/bin/Rlib/")
library("GenABEL", lib="~/bin/Rlib/")

# plink -bfile agpv4_chr3_agpv3_140-160M --snps-only --recode transpose --out agpv4_chr3_agpv3_140-160M
convert.snp.tped(tpedfile="/home/jolyang/dbcenter/HapMap/HapMap3/agpv4_chr3_agpv3_140-160M.tped", 
                tfamfile="/home/jolyang/dbcenter/HapMap/HapMap3/agpv4_chr3_agpv3_140-160M.tfam", 
                outfile="/home/jolyang/dbcenter/HapMap/HapMap3/agpv4_chr3_agpv3_140-160M.raw",
                strand = "u", bcast = 10000)

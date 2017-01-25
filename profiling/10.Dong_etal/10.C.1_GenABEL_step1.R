### Jinliang Yang
### 01-23-2017
### find the common set of the genotype and phenotype data

#### pheno types
getp <- function(){
    file <- list.files(path="largedata/BLUPs", pattern="csv", full.names = TRUE)
    p <- read.csv(file[1])
    names(p)[2] <- gsub(".*/|.csv", "", file[1])
    
    for(i in 2:length(file)){
        tem <- read.csv(file[i])
        names(tem)[2] <- gsub(".*/|.csv", "", file[i])
        p <- merge(p, tem, by="X", all=TRUE)
    }
    
    names(p)[1] <- "id"
    p$id <- toupper(p$id)
}



#### Germplasm GWAS using GenABEL
commonGenoPheno <- function(){
    #### ==>
    pheno <- read.csv("largedata/BLUPs/MainSpikeLength.csv", header=TRUE)
    pheno$X <- toupper(pheno$X)
    #### ==>
    geno <- fread("largedata/10.Dong/ZeaGBSv27_Ames282_agpv3.illumina", sep="\t", header=TRUE)
    geno <- as.data.frame(geno)
    
    names(geno) <- toupper(names(geno))
    idx <- names(geno)[names(geno) %in% pheno$X]
    message(sprintf("###>>> common ids [ %s ]", length(idx)))
    geno <- geno[, c("NAME", "CHR", "POS", idx)]
    names(geno)[1:3] <- c("Name", "Chr", "Pos")
    write.table(geno, "largedata/10.Dong/ZeaGBSv27_278_agpv3.illumina", sep="\t", row.names=FALSE, quote=FALSE)
    
}

library("data.table")
library("GenABEL.data", lib="~/bin/Rlib/")
library("GenABEL", lib="~/bin/Rlib/")

commonGenoPheno()
###>>> common ids [ 264 ]
convert.snp.illumina(infile="largedata/10.Dong/ZeaGBSv27_278_agpv3.illumina", 
                     outfile="largedata/10.Dong/ZeaGBSv27_278_agpv3.raw", strand = "+", bcast = 10000000)






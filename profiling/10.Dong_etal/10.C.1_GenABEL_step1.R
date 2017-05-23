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
    p$sex <- 1
    write.table(p, "largedata/pheno_282_dong_lc.txt", sep="\t", row.names=FALSE, quote=FALSE)
    
    p$id <- toupper(p$id)
    write.table(p, "largedata/pheno_282_dong.txt", sep="\t", row.names=FALSE, quote=FALSE)
}



#### Germplasm GWAS using GenABEL
commonGenoPheno <- function(){
    #### ==>
    pheno <- read.csv("largedata/BLUPs/MainSpikeLength.csv", header=TRUE)
    pheno$X <- toupper(pheno$X)
    #### ==>
    geno <- fread("largedata/10.Dong/ZeaGBSv27_Ames282_agpv3.illumina", data.table=FALSE, sep="\t", header=TRUE)
    
    ### ===> AGPv4 coordinates
    bed <- fread("largedata/output_ZeaGBSv27_Ames282_agpv3.bed", data.table=FALSE)
    newgeno <- merge(bed[, -2], geno[, -2:-3], by.x="V4", by.y="Name")
    newgeno <- subset(newgeno, V1 %in% 1:10)
    names(newgeno)[1:3] <- c("Name", "Chr", "Pos")
    
    names(newgeno) <- toupper(names(newgeno))
    idx <- names(newgeno)[names(newgeno) %in% pheno$X]
    message(sprintf("###>>> common ids [ %s ]", length(idx)))
    newgeno <- newgeno[, c("NAME", "CHR", "POS", idx)]
    names(newgeno)[1:3] <- c("Name", "Chr", "Pos")
    write.table(newgeno, "largedata/10.Dong/ZeaGBSv27_278_agpv4.illumina", sep="\t", row.names=FALSE, quote=FALSE)
    
}

library("data.table")
library("GenABEL.data", lib="~/bin/Rlib/")
library("GenABEL", lib="~/bin/Rlib/")

### up
geno <- fread("largedata/10.Dong/ZeaGBSv27_Ames282_agpv3.illumina", sep="\t", header=TRUE, data.table=FALSE)
bed <- as.data.frame(geno[, 1:3])
bed <- bed[, c("Chr", "Pos", "Pos", "Name")]
bed$Pos <- bed$Pos -1
write.table(bed, "largedata/ZeaGBSv27_Ames282_agpv3.bed", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")

commonGenoPheno()
###>>> common ids [ 264 ]


convert.snp.illumina(infile="largedata/10.Dong/ZeaGBSv27_278_agpv4.illumina", 
                     outfile="largedata/10.Dong/ZeaGBSv27_278_agpv4.raw", strand = "+", bcast = 10000000)






### Jinliang
### May 12th, 2015

#source("~/Documents/Github/zmSNPtools/Rcodes/dsnp2GenABEL.R")
library("data.table")
#library("GenABEL.data", lib="~/bin/Rlib/")
#library("GenABEL", lib="~/bin/Rlib/")

bed2illumina <- function(){
    #######==> GBS data for GenABEL
    gbs <- fread("~/dbcenter/AllZeaGBS/ZeaGBSv27_Ames282_agpv3.bed5", data.table=FALSE)
    
    nms <- names(gbs)
    nms2 <- gsub(":.*", "", nms)
    names(gbs) <- nms2
    
    gbs1 <- gbs
    gbs1 <- gbs1[, c(-2, -5)]
    gbs1 <- gbs1[, c(3,1, 2, 4:ncol(gbs1))]
    names(gbs1)[1:3] <- c("Name", "Chr", "Pos")
    
    gbs1[gbs1 == "A"] <- "AA"
    gbs1[gbs1 == "T"] <- "TT"
    gbs1[gbs1 == "C"] <- "CC"
    gbs1[gbs1 == "G"] <- "GG"
    gbs1[gbs1 == "-"] <- "--"
    gbs1[gbs1 == "N"] <- "00"
    
    message("start to writing ...")
    write.table(gbs1, "largedata/10.Dong/ZeaGBSv27_Ames282_agpv3.illumina", sep="\t", row.names=FALSE, quote=FALSE)
}

bed2illumina()
###>>> Read 509572 rows and 293 (of 293) columns from 0.291 GB file in 00:00:14
###>>> Input [ 509572 ] GBS data, after filtering, [ 306190 ] remaining


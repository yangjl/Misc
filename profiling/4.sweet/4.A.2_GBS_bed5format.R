### Jinliang Yang
### purpose: transform ZeaGBSv2.7 Ames inbred => bed+ format

######=======>
library(data.table, lib = "~/bin/Rlib")

source("lib/gbs2bed_ames282.R")
library(data.table)
gbs2bed_ames(gbsfile="/group/jrigrp4/AllZeaGBSv2.7impV5/ZeaGBSv27_Ames282.hmp.txt",
             outfile="/group/jrigrp4/AllZeaGBSv2.7impV5/ZeaGBSv27_Ames282.bed5")

### log
#>>> Read 955690 rows and 299 (of 299) columns from 0.555 GB file in 00:00:34
#>>> Changed to BED5+ format and start filtering ...
#>>> Remaining [ 509572 ] sites with two variations!
#>>> Start to IUPAC=>N transforming, recoding and writing ...


snpfrq -p /group/jrigrp4/AllZeaGBSv2.7impV5 -i ZeaGBSv27_Ames282.bed5 -s 6 -m "0" -a 0 -b 1 -c 2 -o ZeaGBSv27_Ames282.frq

####### Plot the missing rate and MAF

frq <- fread("/group/jrigrp4/AllZeaGBSv2.7impV5/ZeaGBSv27_Ames282.frq", header=TRUE)
frq <- as.data.frame(frq)
frq$chr <- gsub("_.*", "", frq$snpid)
frq$pos <- as.numeric(as.character(gsub(".*_", "", frq$snpid)))

### plot 
par(mfrow=c(1,2))
hist(frq$MAF, breaks=50, xlab="Minor Allele Frequency", main="MAF of Ames282")
hist(frq$missing, breaks=50, xlab="Missing Rate", main="Missing Rate of Ames282")


dim(subset(frq, MAF > 0.05 & missing < 0.8))
### [1] 319531      5

dim(subset(frq, chr == "S5" & pos > 127464000 & pos < 127468000))
# 3
dim(subset(frq, chr == "S5" & pos > 127000000 & pos < 128000000))
# 135


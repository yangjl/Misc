### Jinliang Yang
### May 12th, 2015
### conduct MLM GWAS

library("data.table", lib="~/bin/Rlib/")
library("GenABEL.data", lib="~/bin/Rlib/")
library("GenABEL", lib="~/bin/Rlib/")

### load the data
gm <- load.gwaa.data(phe="data/pheno_ames282.txt", gen="largedata/4.sweet/geno_ames264.raw", force=T)
head(gm@phdata)
#gm@gtdata

################# QC1 #########################
qc1 <- check.marker(gm, p.level=0, callrate=0.5, ibs.mrk=0, maf=0.05, perid.call=0.5)
# 0 people excluded because too high autosomal heterozygosity (FDR <1%)
# In total, 304009 (100%) markers passed all criteria
# In total, 264 (100%) people passed all criteria

### using the simple model GWAS
frq <- fread("largedata/4.sweet/ZeaGBSv27_Ames264.illumina", header=TRUE)
frq <- as.data.frame(frq)
frq <- frq[, 1:3]
### gene located on 127466k region.
chr5snp <- subset(frq, Chr == 5 &Pos > 122000000 & Pos < 132000000)$Name
#349
gm.qc1 <- gm[qc1$idok, chr5snp]

############### simple model ##########################
par(mfrow=c(1,2)) 
#chr3:166596669-166598742     RefGen_v2 	B73 

res1 <- qtscore(X10KW, data=gm, trait = "gaussian" )
res2 <- qtscore(TKW, data=gm, trait = "gaussian")
res3 <- qtscore(KC, data=gm, trait = "gaussian" )
res4 <- qtscore(CD, data=gm, trait = "gaussian")
res5 <- qtscore(CL, data=gm, trait = "gaussian" )
res6 <- qtscore(CW, data=gm, trait = "gaussian")

plot(res1)
abline(v=127466000)
plot(res2)
abline(v=127466000)

plot(res3)
plot(res5)
plot(res6)
############### kinship ##########################
gkin <- ibs(gm, weight="freq")
res1.eg <- egscore(X10KW, data=gm, kinship.matrix=gkin)
res2.eg <- egscore(TKW, data=gm, kinship.matrix=gkin)
res3.eg <- egscore(KC, data=gm, kinship.matrix=gkin)
res4.eg <- egscore(CD, data=gm, kinship.matrix=gkin)
res5.eg <- egscore(CL, data=gm, kinship.matrix=gkin)
res6.eg <- egscore(CW, data=gm, kinship.matrix=gkin)

plot(res1.eg)
plot(res2.eg) ## cool
res2.eg <- checkpeak(res = res2.eg)
plot(res3.eg)
res3.eg <- checkpeak(res = res3.eg)
plot(res4.eg)
plot(res5.eg)
plot(res6.eg)


checkpeak <- function(res = res3){
    r2 <- results(res)
    nrow(subset(r2, Chromosome==3 & Position > 166500000 & Position < 166700000))
    r2 <- r2[order(r2$chi2.1df, decreasing=TRUE),]
    return(r2)
}
############### MLM #################
h2a1 <- polygenic(X10KW, data=gm, kin=gkin)
res1.mm <- mmscore(h2a1, data=gm)
h2a2 <- polygenic(TKW, data=gm, kin=gkin)
res2.mm <- mmscore(h2a2, data=gm)
h2a3 <- polygenic(KC, data=gm, kin=gkin)
res3.mm <- mmscore(h2a3, data=gm)
h2a4 <- polygenic(CD, data=gm, kin=gkin)
res4.mm <- mmscore(h2a4, data=gm)
h2a5 <- polygenic(CL, data=gm, kin=gkin)
res5.mm <- mmscore(h2a5, data=gm)
h2a6 <- polygenic(CW, data=gm, kin=gkin)
res6.mm <- mmscore(h2a6, data=gm)

plot(res1.mm, main="10 kernel weight", pch=16)
plot(res2.mm, main="total kernel weight", cex=0.5, pch=16)
abline(h=5, col="red", lty=2)
res2.mm <- checkpeak(res = res2.mm)
plot(res3.mm, main="Kernel Count", cex=0.5, pch=16)
abline(h=5, col="red", lty=2)
res3.mm <- checkpeak(res = res3.mm)
plot(res4.mm, main="CD", pch=16)
plot(res5.mm, main="CL", pch=16)
plot(res6.mm, main="CW", pch=16)








save(file="cache/gwas_res.RData", list=c("gm", "res1", "res2", "res1.eg", "res2.eg", "res1.mm", "res2.mm"))




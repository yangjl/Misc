### Jinliang Yang
### May 12th, 2015
### conduct MLM GWAS

library("data.table", lib="~/bin/Rlib/")
library("GenABEL.data", lib="~/bin/Rlib/")
library("GenABEL", lib="~/bin/Rlib/")




### load the data
gm1 <- load.gwaa.data(phe="largedata/pheno_282_dong.txt", 
                     gen="largedata/10.Dong/ZeaGBSv27_278_agpv4.raw", force=T)
head(gm@phdata)
#gm@gtdata


### load the data
gm2 <- load.gwaa.data(phe="largedata/pheno_282_dong_lc.txt", 
                     gen="/home/jolyang/dbcenter/HapMap/HapMap3/agpv4_chr3_agpv3_151M_282set.raw", force=T)
head(gm@phdata)
gm2@phdata$id <- toupper(gm2@phdata$id)
gm2@gtdata@idnames <- toupper(gm2@gtdata@idnames)

gm <- merge.gwaa.data(gm1, gm2)


################# QC1 #########################
qc1 <- check.marker(gm, p.level=0, callrate=0.5, ibs.mrk=0, maf=0.05, perid.call=0.5)
# 0 people excluded because too high autosomal heterozygosity (FDR <1%)
# In total, 310943 (100%) markers passed all criteria
# In total, 278 (100%) people passed all criteria



### using the simple model GWAS
tem <- data.frame(snpid=qc1$call$name, chr=qc1$call$chromosome, pos=qc1$call$map)
tb <- subset(tem, snpid %in% qc1$snpok)
# GRMZM2G039867 tru1
# AGPv4 annotation
# (Chr3: 151329451..151332389)
BP <- 500000
chr3snp <- as.character(subset(tb, chr == 3 & pos > 151329451 - BP & pos < 151332389 + BP)$snpid)
# length(chr3snp)

#207
gm.qc <- gm[qc1$idok, qc1$snpok]
gm.qc.chr3 <- gm[qc1$idok, chr3snp]



############### kinship ##########################
gkin <- ibs(gm.qc, weight="freq")


############### MLM #################
head(gm@phdata)
# MainSpikeLength
t1 <- polygenic(MainSpikeLength, data=gm.qc, kin=gkin)
res1.mm <- mmscore(t1, data=gm.qc, snpsubset= chr3snp)

# NumberofTilleringPlants
t2 <- polygenic(NumberofTilleringPlants, data=gm.qc, kin=gkin)
res2.mm <- mmscore(t2, data=gm.qc, snpsubset= chr3snp)

# SecondaryBranchNumber
t3 <- polygenic(SecondaryBranchNumber, data=gm.qc, kin=gkin)
res3.mm <- mmscore(t3, data=gm.qc, snpsubset= chr3snp)

# Spikelets.MainSpike
t4 <- polygenic(Spikelets.MainSpike, data=gm.qc, kin=gkin)
res4.mm <- mmscore(t4, data=gm.qc, snpsubset= chr3snp)

# Spikelets.PrimaryBranch
t5 <- polygenic(Spikelets.PrimaryBranch, data=gm.qc, kin=gkin)
res5.mm <- mmscore(t5, data=gm.qc, snpsubset= chr3snp)

# TasselBranchLength
t6 <- polygenic(TasselBranchLength, data=gm.qc, kin=gkin)
res6.mm <- mmscore(t6, data=gm.qc, snpsubset= chr3snp)

# TasselLength
t7 <- polygenic(TasselLength, data=gm.qc, kin=gkin)
res7.mm <- mmscore(t7, data=gm.qc, snpsubset= chr3snp)

# TasselPrimaryBranches
t8 <- polygenic(TasselPrimaryBranches, data=gm.qc, kin=gkin)
res8.mm <- mmscore(t8, data=gm.qc, snpsubset= chr3snp)


save(file="largedata/dong_gwas_agpv4.RData", 
     list=c("gm.qc","res1.mm", "res2.mm", "res3.mm", "res4.mm", "res5.mm", "res6.mm", "res7.mm", "res8.mm"))


res1 <- mmscore(t1, data=gm.qc)
res2 <- mmscore(t2, data=gm.qc)
res3 <- mmscore(t3, data=gm.qc)
res4 <- mmscore(t4, data=gm.qc)
res5 <- mmscore(t5, data=gm.qc)
res6 <- mmscore(t6, data=gm.qc)
res7 <- mmscore(t7, data=gm.qc)
res8 <- mmscore(t8, data=gm.qc)

save(file="largedata/dong_gwas_agpv4_all.RData", 
     list=c("gm.qc","res1", "res2", "res3", "res4", "res5", "res6", "res7", "res8"))


########################################
load("largedata/dong_gwas_agpv4.RData")

#plot(res1.mm, main="MainSpikeLength", pch=16, col="cadetblue")
plot(res2.mm, main="NumberofTilleringPlants", pch=16, col="cadetblue")
#plot(res2.mm, df="Pc1df", main="NumberofTilleringPlants", pch=16, col="cadetblue")
abline(v=151332389)

res2 <- results(res2.mm)
#plot(res3.mm, main="10 kernel weight", pch=16, col="cadetblue")

plot(res4.mm, main="Spikelets.MainSpike", pch=16, col="cadetblue")
abline(v=151332389)

#plot(res5.mm, main="Spikelets.PrimaryBranch", pch=16, col="cadetblue")
#plot(res6.mm, main="TasselBranchLength", pch=16, col="cadetblue")
#plot(res7.mm, main="10 kernel weight", pch=16, col="cadetblue")
#plot(res8.mm, main="10 kernel weight", pch=16, col="cadetblue")
#abline(v=151332389)


#151329451..151332389


########################################
load("largedata/dong_gwas_agpv4_all.RData")

plot(res1, ystart=2, main="MainSpikeLength", pch=16, cex=0.5, col="cadetblue")
plot(res2, ystart=2, cex=0.5, main="NumberofTilleringPlants", pch=16, col="cadetblue")
plot(res3, ystart=2, cex=0.5, main="10 kernel weight", pch=16, col="cadetblue")
plot(res4, ystart=2, cex=0.5,main="Spikelets.MainSpike", pch=16, col="cadetblue")
plot(res5, ystart=2, cex=0.5, main="Spikelets.PrimaryBranch", pch=16, col="cadetblue")
plot(res6,ystart=2, cex=0.5, main="TasselBranchLength", pch=16, col="cadetblue")
plot(res7, ystart=2, cex=0.5,main="10 kernel weight", pch=16, col="cadetblue")
plot(res8, ystart=2, cex=0.5, main="10 kernel weight", pch=16, col="cadetblue")
abline(v=151332389)


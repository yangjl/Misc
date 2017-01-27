### Jinliang Yang
### May 12th, 2015
### conduct MLM GWAS

library("data.table", lib="~/bin/Rlib/")
library("GenABEL.data", lib="~/bin/Rlib/")
library("GenABEL", lib="~/bin/Rlib/")




### load the data
gm <- load.gwaa.data(phe="largedata/pheno_282_dong.txt", 
                     gen="largedata/10.Dong/ZeaGBSv27_278_agpv3.raw", force=T)
head(gm@phdata)
#gm@gtdata

################# QC1 #########################
qc1 <- check.marker(gm, p.level=0, callrate=0.5, ibs.mrk=0, maf=0.05, perid.call=0.5)
# 0 people excluded because too high autosomal heterozygosity (FDR <1%)
# In total, 310943 (100%) markers passed all criteria
# In total, 278 (100%) people passed all criteria

### using the simple model GWAS
tem <- data.frame(snpid=qc1$call$name, chr=qc1$call$chromosome, pos=qc1$call$map)
tb <- subset(tem, snpid %in% qc1$snpok)
# GRMZM2G039867 tru1
BP <- 1000000
chr3snp <- as.character(subset(tb, chr == 3 & pos > 150088574 - BP & pos < 150091550 + BP)$snpid)
# length(chr3snp)

#207
gm.qc <- gm[qc1$idok, qc1$snpok]
gm.qc.chr3 <- gm[qc1$idok, chr3snp]

############### simple model ##########################
par(mfrow=c(1,2))
res1 <- qtscore(CobDiameter, data=gm.qc.chr3, trait = "gaussian" )
res2 <- qtscore(TKW, data=gm.qc1, trait = "gaussian")

plot(res1)
abline(v=127466000)
plot(res2)
abline(v=127466000)

############### kinship ##########################
gkin <- ibs(gm.qc, weight="freq")
res1.eg <- egscore(X10KW, data=gm.qc1, kinship.matrix=gkin)
res2.eg <- egscore(TKW, data=gm.qc1, kinship.matrix=gkin)

plot(res1.eg)
abline(v=127466000)
plot(res2.eg)
abline(v=127466000)

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








save(file="largedata/dong_gwas_res.RData", 
     list=c("gm.qc","res1.mm", "res2.mm", "res3.mm", "res4.mm", "res5.mm", "res6.mm", "res7.mm", "res8.mm"))


load("largedata/dong_gwas_res.RData")

#plot(res1.mm, main="10 kernel weight", pch=16, col="cadetblue")
plot(res2.mm, main="NumberofTilleringPlants", pch=16, col="cadetblue")
#plot(res3.mm, main="10 kernel weight", pch=16, col="cadetblue")
plot(res4.mm, main="Spikelets.MainSpike", pch=16, col="cadetblue")
plot(res5.mm, main="Spikelets.PrimaryBranch", pch=16, col="cadetblue")
#plot(res6.mm, main="TasselBranchLength", pch=16, col="cadetblue")
#plot(res7.mm, main="10 kernel weight", pch=16, col="cadetblue")
#plot(res8.mm, main="10 kernel weight", pch=16, col="cadetblue")





abline(v=127466000)
plot(res2.mm, main="total kernel weight", pch=16, col="cadetblue")
abline(v=127466000)


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

######### heritability and effect size #################
geth2 <- function(snpid = "S5_127255673", gm=gm, res=res1.mm,
                  traits=c("TasselPrimaryBranches", "SecondaryBranchNumber")){
    snp <- as.numeric(gtdata(gm[, snpid]))
    #table(snp,exclude=NULL)
    # introduce missing data indicator
    nomissInd <- (!is.na(snp))
    #table(nomissInd)
    gkin <- ibs(gm[nomissInd,], w="freq")
    
    pol1 <- polygenic(TasselPrimaryBranches, gm[nomissInd,], kin=gkin, quiet=TRUE)
    pol2 <- polygenic(SecondaryBranchNumber, gm[nomissInd,], kin=gkin, quiet=TRUE)
    
    pol <- polygenic(NumberofTilleringPlants, gm[nomissInd,], kin=gkin, quiet=TRUE, llfun='polylik')
    #pol1$esth2
    #pol0 <- polygenic(TKW, gm[nomissInd,], kin=gkin, quiet=TRUE)
    #pol0$esth2
    
    ### ===> heritability
    propVarExBySNP1 <- res[snpid,"chi2.1df"]/res[snpid,"N"]
    propHerExBySNP1 <- propVarExBySNP1/pol1$esth2
    propVarExBySNP2 <- res[snpid,"chi2.1df"]/res[snpid,"N"]
    propHerExBySNP2 <- propVarExBySNP2/pol2$esth2
    res <- data.frame(snp=snpid, trait=traits, 
                      varp=c(propVarExBySNP1, propVarExBySNP2), 
                      varh2=c(propHerExBySNP1, propHerExBySNP2),
                      h2=c(pol1$esth2, pol2$esth2))
    return(res)
}

################# get heritability #########################
load("largedata/dong_gwas_agpv4.RData")
res2 <- subset(results(res2.mm), !is.na(P1df))
idx <- which.max(-log10(res2$Pc1df))
res2[idx,]

geth2(snpid = "3_150058036_G_A", gm=gm, res=res2.mm, traits="NumberofTilleringPlants")
geth2(snpid = "3_150053905_T_A", gm=gm, res=res3.mm)
geth2(snpid = "3_150095715_C_T", gm=gm, res=res8.mm)


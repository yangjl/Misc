### Jinliang Yang
### May 12th, 2015
###

#### Germplasm GWAS using GenABEL


commonGenoPheno <- function(){
    #### ==>
    pheno <- read.table("largedata/4.sweet/pheno_ames282.txt", header=TRUE)
    
    #### ==>
    geno <- fread("largedata/4.sweet/ZeaGBSv27_Ames282.illumina", sep="\t", header=TRUE)
    geno <- as.data.frame(geno)
    
    names(geno) <- toupper(names(geno))
    idx <- names(geno)[names(geno) %in% pheno$id]
    message(sprintf("###>>> common ids [ %s ]", length(idx)))
    geno <- geno[, c("NAME", "CHR", "POS", idx)]
    names(geno)[1:3] <- c("Name", "Chr", "Pos")
    write.table(geno, "largedata/4.sweet/ZeaGBSv27_Ames264.illumina", sep="\t", row.names=FALSE, quote=FALSE)
    
}

library("data.table", lib="~/bin/Rlib/")
library("GenABEL.data", lib="~/bin/Rlib/")
library("GenABEL", lib="~/bin/Rlib/")

commonGenoPheno()
###>>> common ids [ 264 ]
convert.snp.illumina(infile="largedata/4.sweet/ZeaGBSv27_Ames264.illumina", 
                     outfile="largedata/4.sweet/geno_ames264.raw", strand = "+", bcast = 10000000)


gm <- load.gwaa.data(phe="largedata/4.sweet/pheno_ames282.txt", gen="largedata/4.sweet/geno_ames264.raw", force=T)
head(gm@phdata)
#gm@gtdata


frq <- fread("largedata/4.sweet/ZeaGBSv27_Ames264.illumina", header=TRUE)
frq <- as.data.frame(frq)
frq <- frq[, 1:3]
chr5snp <- subset(frq, Chr == 5 &Pos > 122000000 & Pos < 132000000)$Name
################# QC1 #########################
qc1 <- check.marker(gm, p.level=0, callrate=0.5, ibs.mrk=0, maf=0.05, perid.call=0.5)
gm.qc1 <- gm[qc1$idok, chr5snp]


res1 <- qtscore(X10KW, data=gm.qc1, trait = "gaussian" )
plot(res1)

res2 <- qtscore(KC, data=gm.qc1)
plot(res2)

############### kinship ##########################
gkin <- ibs(gm, weight="freq")
res2.eg <- egscore(X10KW, data=gm.qc1, kinship.matrix=gkin)
plot(res2.eg)

res2.eg <- egscore(TKW, data=gm.qc1, kinship.matrix=gkin)
plot(res2.eg)

res2.eg <- egscore(KC, data=gm.qc1, kinship.matrix=gkin)
plot(res2.eg)

############### MLM #################
h2a <- polygenic(TKW, data=gm, kin=gkin)
res3.mm <- mmscore(h2a, data=gm, snpsubset=chr5snp)
plot(res3.mm)

gw.mm <- mmscore(h2a, data=gm)
plot(gw.mm)

summary(glm(TKW ~ S5_127255673, data=gm.qc1))

######### heritability

snp <- as.numeric(gtdata(gm[,"S5_127255673"]))
table(snp,exclude=NULL)
# introduce missing data indicator
nomissInd <- (!is.na(snp))
table(nomissInd)

gkin <- ibs(gm[nomissInd,], w="freq")
pol1 <- polygenic(TKW~snp[nomissInd], gm[nomissInd,], kin=gkin, quiet=TRUE)
pol1$esth2
pol0 <- polygenic(TKW, gm[nomissInd,], kin=gkin, quiet=TRUE)
pol0$esth2

### ===> heritability
propVarExBySNP <- res3.mm["S5_127255673","chi2.1df"]/res3.mm["S5_127255673","N"]
propHerExBySNP <- propVarExBySNP/h2a$esth2

pheno <- phdata(gm)
df <- merge(pheno, snp, by="row.names")
df <- subset(df, !is.na(S5_127255673))
boxplot(TKW ~ as.factor(S5_127255673), data=df)

t.test(subset(df, S5_127255673==0)$TKW, subset(df, S5_127255673==2)$TKW)

### Jinliang Yang
### May 14th, 2015
### plot the simple results and compute the size of effect and heritability

library("data.table", lib="~/bin/Rlib/")
library("GenABEL.data", lib="~/bin/Rlib/")
library("GenABEL", lib="~/bin/Rlib/")

load("cache/gwas_res.RData")
res1 <- results(res1)
res1$qval <- p.adjust(res1$P1df, method = "fdr")
res2 <- results(res2)
res2$qval <- p.adjust(res2$P1df, method = "fdr")
res1.eg <- results(res1.eg)
res2.eg <- results(res2.eg)
res1.mm <- results(res1.mm)
res2.mm <- results(res2.mm)

write.table(res1.mm[, 1:11], "data/Table_gwas_10kw.txt", sep="\t", quote=FALSE)
write.table(res2.mm[, 1:11], "data/Table_gwas_tkw.txt", sep="\t", quote=FALSE)

par(mfrow=c(1,2))
plot(x=res1.mm$Position, y=-log10(res1.mm$Pc1df), xlab="Chromosome 5", main="10 Kernel Weight",
     ylab="-log10(q-value)", pch=19, col="cadetblue")
abline(v=127466000, col="blue", lty=2, lwd=2)
abline(h=-log10(0.01), col="red", lty=2, lwd=2)
#S5_128108485

plot(x=res2.mm$Position, y=-log10(res2.mm$Pc1df), xlab="Chromosome 5", main="Total Kernel Weight",
     ylab="-log10(q-value)", pch=19, col="cadetblue")
abline(v=127466000, col="blue", lty=2, lwd=2)
abline(h=-log10(0.01), col="red", lty=2, lwd=2)
#S5_127255673

######### heritability and effect size #################
geth2 <- function(snpid = "S5_127255673", gm=gm, res1=res1.mm, model=as.formula("TKW")){
    snp <- as.numeric(gtdata(gm[, snpid]))
    #table(snp,exclude=NULL)
    # introduce missing data indicator
    nomissInd <- (!is.na(snp))
    #table(nomissInd)
    gkin <- ibs(gm[nomissInd,], w="freq")
    
    pol1 <- polygenic(TKW, gm[nomissInd,], kin=gkin, quiet=TRUE)
    pol2 <- polygenic(X10KW, gm[nomissInd,], kin=gkin, quiet=TRUE)
    #pol1$esth2
    #pol0 <- polygenic(TKW, gm[nomissInd,], kin=gkin, quiet=TRUE)
    #pol0$esth2
    
    ### ===> heritability
    propVarExBySNP1 <- res1[snpid,"chi2.1df"]/res1[snpid,"N"]
    propHerExBySNP1 <- propVarExBySNP1/pol1$esth2
    propVarExBySNP2 <- res1[snpid,"chi2.1df"]/res1[snpid,"N"]
    propHerExBySNP2 <- propVarExBySNP2/pol1$esth2
    res <- data.frame(snp=snpid, trait=c("TKW", "X10KW"), 
                      varp=c(propVarExBySNP1, propVarExBySNP1), 
                      varh2=c(propHerExBySNP1, propHerExBySNP1) )
    return(res)
}

geth2(snpid = "S5_127255673", gm=gm, res1=res1.mm, model=as.formula("TKW"))





propVarExBySNP <- res2.mm["S5_127255673","chi2.1df"]/res2.mm["S5_127255673","N"]
propHerExBySNP <- propVarExBySNP/h2a1$esth2

pheno <- phdata(gm)
df <- merge(pheno, snp, by="row.names")
df <- subset(df, !is.na(S5_127255673))
boxplot(TKW ~ as.factor(S5_127255673), data=df)

t.test(subset(df, S5_127255673==0)$TKW, subset(df, S5_127255673==2)$TKW)
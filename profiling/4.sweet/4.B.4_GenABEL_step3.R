### Jinliang Yang
### May 14th, 2015
### plot the simple results and compute the size of effect and heritability

library("data.table", lib="~/bin/Rlib/")
library("GenABEL.data", lib="~/bin/Rlib/")
library("GenABEL", lib="~/bin/Rlib/")

load("largedata/dong_gwas_res.RData")

#plot(res1.mm, main="10 kernel weight", pch=16, col="cadetblue")
plot(res2.mm, main="NumberofTilleringPlants", pch=16, col="cadetblue")
#plot(res3.mm, main="10 kernel weight", pch=16, col="cadetblue")
plot(res4.mm, main="Spikelets.MainSpike", pch=16, col="cadetblue")
plot(res5.mm, main="Spikelets.PrimaryBranch", pch=16, col="cadetblue")
#plot(res6.mm, main="TasselBranchLength", pch=16, col="cadetblue")
#plot(res7.mm, main="10 kernel weight", pch=16, col="cadetblue")
#plot(res8.mm, main="10 kernel weight", pch=16, col="cadetblue")

res2 <- results(res2.mm)
res2$qval <- p.adjust(res2$P1df, method = "fdr")

res4 <- results(res4.mm)
res5 <- results(res5.mm)

write.table(res1.mm[, 1:11], "data/Table_gwas_10kw.txt", sep="\t", quote=FALSE)
write.table(res2.mm[, 1:11], "data/Table_gwas_tkw.txt", sep="\t", quote=FALSE)

pdf("graphs/Figure_results.pdf", width=10, height=5)
par(mfrow=c(1,3))
plot(x=res2$Position, y=-log10(res2$Pc1df), xlab="Chromosome 3", main="Number of Tillering Plants",
     ylab="-log10(q-value)", pch=19, cex=1.3, col="cadetblue")
abline(v=150088574, col="grey", lty=2, lwd=2)
abline(v=150091550, col="grey", lty=2, lwd=2)
abline(h=-log10(0.01), col="red", lty=2, lwd=1)
#S5_128108485 B73=T

plot(x=res4$Position, y=-log10(res4$Pc1df), xlab="Chromosome 3", main="Spikelets Main Spike",
     ylab="-log10(q-value)", pch=19, cex=1.3, col="cadetblue")
abline(v=150088574, col="grey", lty=2, lwd=2)
abline(v=150091550, col="grey", lty=2, lwd=2)
abline(h=-log10(0.01), col="red", lty=2, lwd=1)
#S5_127255673 B73=A

plot(x=res5$Position, y=-log10(res5$Pc1df), xlab="Chromosome 3", main="Spikelets Primary Branch",
     ylab="-log10(q-value)", pch=19, cex=1.3, col="cadetblue")
abline(v=150088574, col="grey", lty=2, lwd=2)
abline(v=150091550, col="grey", lty=2, lwd=2)
abline(h=-log10(0.01), col="red", lty=2, lwd=1)
dev.off()

######### heritability and effect size #################
geth2 <- function(snpid = "S5_127255673", gm=gm, res1=res1.mm){
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
    propHerExBySNP2 <- propVarExBySNP2/pol2$esth2
    res <- data.frame(snp=snpid, trait=c("TKW", "X10KW"), 
                      varp=c(propVarExBySNP1, propVarExBySNP2), 
                      varh2=c(propHerExBySNP1, propHerExBySNP2),
                      h2=c(pol1$esth2, pol2$esth2))
    return(res)
}

h21 <- geth2(snpid = "S5_128108485", gm=gm, res1=res1.mm)
h21 <- h21[-1,]

h22 <- geth2(snpid = "S5_127255673", gm=gm, res1=res2.mm)
h22 <- h22[-2,]
plotbox_1 <- function(){
    snp <- as.numeric(gtdata(gm[, "S5_127255673"]))
    pheno <- phdata(gm)
    df <- merge(pheno, snp, by="row.names")
    df <- subset(df, !is.na(S5_127255673))
    boxplot(TKW ~ as.factor(S5_127255673), data=df, names=c("B73-like","non-B73"))
    print(t.test(subset(df, S5_127255673==0)$TKW, subset(df, S5_127255673==2)$TKW))
}
plotbox_1()

plotbox_2 <- function(){
    snp <- as.numeric(gtdata(gm[, "S5_128108485"]))
    pheno <- phdata(gm)
    df <- merge(pheno, snp, by="row.names")
    df <- subset(df, !is.na(S5_128108485))
    boxplot(X10KW ~ as.factor(S5_128108485), data=df, names=c("B73-like","non-B73"))
}
plotbox_2()
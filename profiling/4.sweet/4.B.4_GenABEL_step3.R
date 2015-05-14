### Jinliang Yang
### May 14th, 2015
### plot the simple results

plot(res1.mm, main="10 kernel weight", pch=16, col="cadetblue")
abline(v=127466000)
plot(res2.mm, main="total kernel weight", pch=16, col="cadetblue")
abline(v=127466000)

######### heritability and effect size #################
summary(glm(TKW ~ S5_127255673, data=gm.qc1))

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
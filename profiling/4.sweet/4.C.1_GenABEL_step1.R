### Jinliang Yang
### May 12th, 2015
###

#### Germplasm GWAS using GenABEL

setwd("/Users/yangjl/Documents/workingSpace/NAMGWAS/Verification/GWAS_germplasm")
load("Germplasm_formatting.RData")
ls()

library("GenABEL")
convert.snp.illumina(inf="geno_transformed.txt", out="germ.geno", strand="file")

gm <- load.gwaa.data(phe="pheno_germplasm.txt", gen="germ.geno", force=T)
gm@phdata
gm@gtdata

################# QC1 #########################
qc1 <- check.marker(gm, p.level=0, callrate=0.5, ibs.mrk=0, maf=0.05, perid.call=0.5)
summary(gm.qc1)

qc1$snpok
SNP67
mysnpok <- SNP67[SNP67 %in% qc1$snpok == TRUE]

gm.qc1 <- gm[qc1$idok, qc1$mysnpok]
res1 <- qtscore(krn ~ GrinID, data=gm.qc1, trait = "binomial" )
plot(res1)

res2 <- qtscore(mean, data=gm.qc1)
plot(res2)


pval1 <- results(res1)
pval1$Pbonf <- p.adjust(pval1$P1df, "bonferroni")
pval1 <- pval1[,c("Chromosome","Position","P1df", "Pbonf")]

pval1$snpid <- paste(pval1$Chromosome, pval1$Position, sep="_")
pval1 <- merge(pval1, subsnp82, by="snpid")
pval1[order(pval1$P1df),]

plot(qc2.pca, main="GWAS by controlling structure")
abline(h=1.2, lty=2)


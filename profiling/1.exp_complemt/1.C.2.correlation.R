### Jinliang Yang
### Jan. 31th, 2015
### test the correlationship between exp and phenotype

source("profiling/1.exp_complemt/1.C.1.count_PAE.R")

#### RNAseq expression data for FGSv2, with GERPavg and 220 accessions
expdata <- getexpdata()

#### count number of complementation weighted by mean GERP score  
res <- countexp(exp=expdata, present=30, absent=0.01, method="sum")

#### get the phenotypic data with GERP weighted # of complementary genes
pheno <- getpheno(res=res)

### calculate the correlations of the two
out <- getcor(pheno=pheno)







###
myres

### remove outlier:
pheno <- subset(pheno, comp < 600)
par(mfrow=c(2,3))
for(i in 3:7){
    x <- pheno$comp
    y <- as.numeric(as.character(pheno[,i]))
    df <- data.frame(x=x, y=y)
    df <- subset(df, !is.na(x) & !is.na(y))
    print(cor.test(df$x, df$y))
    plot(df$x, df$y, xlab="Number of complementary genes",
         ylab="Phenotypic value", main=names(pheno)[i], col="blue")
    reg1 <- lm(df$y ~ df$x)
    abline(reg1, col="red", lwd=2)
}







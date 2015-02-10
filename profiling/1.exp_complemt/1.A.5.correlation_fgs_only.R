### Jinliang Yang
### Feb. 1st, 2015
### correlation using the FGS only

#### Ames Hybrids
hyb <- read.table("data/Ames_hybrids_BLUEs.txt", header=TRUE)
names(hyb)[1] <- "F1"
hyb$P1 <- gsub("x.*", "", hyb$F1)
hyb$P2 <- gsub(".*x", "", hyb$F1)

hyb$P2 <- toupper(gsub(" ", "", hyb$P2))




###
res <- read.csv("cache/exp_res_FGS.csv")
res$comp <- res$comp1 + res$comp2

tem1 <- merge(hyb, res[, c("accenumb", "comp")], by.x="P2", by.y="accenumb")
tem2 <- merge(hyb, res[, c("geno", "comp")], by.x="P2", by.y="geno")
pheno <- rbind(tem1, tem2)
pheno <- subset(pheno, P1 == "PI601322")

pheno[pheno[,3]<0, ][,3] <- "NA"

####### original
for(i in 3:7){
    print(cor.test(pheno$comp, as.numeric(as.character(pheno[,i]))))
    plot(pheno$comp, as.numeric(as.character(pheno[,i]))) 
}


### remove outlier:
pheno <- subset(pheno, comp < 100 & comp > 50)

par(mfrow=c(2,3))
for(i in 3:7){
    print(cor.test(pheno$comp, as.numeric(as.character(pheno[,i]))))
    plot(pheno$comp, as.numeric(as.character(pheno[,i])), xlab="Number of complementary genes",
         ylab="Phenotypic value", main=names(pheno)[i], col="blue") 
}

### remove outlier:

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










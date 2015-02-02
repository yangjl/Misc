### Jinliang Yang
### Jan. 31th, 2015
### test the correlationship between exp and phenotype


#### Ames Hybrids
hyb <- read.table("data/Ames_hybrids_BLUEs.txt", header=TRUE)
names(hyb)[1] <- "F1"
hyb$P1 <- gsub("x.*", "", hyb$F1)
hyb$P2 <- gsub(".*x", "", hyb$F1)

hyb$P2 <- toupper(gsub(" ", "", hyb$P2))




###
res <- read.csv("cache/exp_res_p50a1.csv")

res <- read.csv("cache/exp_res_p5a0.1.csv")
res <- read.csv("cache/exp_res_p10a1.csv")
res$comp <- res$comp1 + res$comp2



###### loop to selet the parameters
for(k in 5:100){
    res <- countexp(exp=exp, present=k, absent=0.1)
    res <- merge(xpid, res, by.x="geno", by.y="plantid")
    
    
    tem1 <- merge(hyb, res[, c("accenumb", "comp")], by.x="P2", by.y="accenumb")
    tem2 <- merge(hyb, res[, c("geno", "comp")], by.x="P2", by.y="geno")
    pheno <- rbind(tem1, tem2)
    pheno <- subset(pheno, P1 == "PI601322")
    
    pheno[pheno[,3]<0, ][,3] <- "NA"
    x <- pheno$comp
    y <- as.numeric(as.character(pheno[,i]))
    if(cor.test(df$x, df$y)$p.value < 0.05){
        print(k)
        stop("I find it!")
    }
}


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







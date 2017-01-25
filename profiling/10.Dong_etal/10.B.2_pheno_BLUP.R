### Jinliang Yang
### 01-22-2017
### study NAM and 282 phenotypes

calt <- function(mypheno, outdir="largedata/BLUPs/", t="CobDiameter", plot=TRUE){
    
    myt <- subset(mypheno, trait %in% t)
    loc <- data.frame(table(myt$location))
    loc <- subset(loc, Freq > 0)
    myt$genotype <- as.factor(myt$genotype)
    myt$location <- as.factor(myt$location)
    myt$value <- as.numeric(myt$value)
    
    if(plot){
        hist(myt$value, main=t)
    }
    
    fit <- get_BLUP(data = myt, model = value ~ (1|genotype) + (1|location), 
                    which.factor = "genotype", outfile = paste0(outdir, t, ".csv"))
    
    h2 <- get_H2(fit, numerator  ="genotype", 
                 denominator = data.frame(f = c("genotype", "Residual"), df = c(1, nrow(loc))))
    return(h2)
}

#####
library(g3tools)
library(lme4)

pheno <- read.csv("largedata/pheno_long.csv")
tb <- data.frame(table(pheno$trait))
idx <- grep("Z0..E0...", pheno$genotype)
mypheno <- pheno[-idx,]

toi <- tb[c(26, 37, 41, 46, 47, 49, 50, 51), ]

toi$H2 <- 0
for(k in 1:nrow(toi)){
    toi$H2[k] <- calt(mypheno, outdir="largedata/BLUPs/", t=toi$Var1[k], plot=TRUE)
}



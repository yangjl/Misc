### Jinliang Yang
### May 12th, 2015
### 


nathan <- read.csv("data/Nathan_p1_subset.csv")

genotype <- unique(nathan$INBRED)



nathan <- read.csv("data/Nathan_p1_subset.csv")
nathan <- nathan[nathan$X10KW_Inbred !=".",]
nathan$X10KW_Inbred <- as.numeric(as.character(nathan$X10KW_Inbred))/10
names(nathan)[7] <- "AKW_Inbred" 

akwblue <- getblue(amesdata=ear, panzeadata=nathan, trait="AKW", model=AKW~Genotype, nathan=TRUE)

library(nlme)
BLUE <- function(data=krn, model=KRN~Genotype, random=~1|Datasource/Farm/Year, trait="KRN", intercept="B73"){
    lmeout1 <- lme(model, data=data, random=random);
    ped.hat1 <- lmeout1$coef$fixed;
    ped.hat1[-1] <- ped.hat1[-1]+ped.hat1[1];
    names(ped.hat1)[1]=intercept;
    names(ped.hat1) <- gsub("Genotype", "", names(ped.hat1));
    tped <- data.frame(Genotype=names(ped.hat1), trait=ped.hat1)
    names(tped)[2] <- trait;
    return(tped)
}


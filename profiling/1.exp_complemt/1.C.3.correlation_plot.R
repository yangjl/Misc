### Jinliang Yang
### Feb. 1st, 2015
### correlation using the FGS only


parplot1 <- function(out1=out1){
    par(mfrow=c(1, 2))
    traits <- unique(out1$trait)
    cols <- c("blue", "black", "wheat", "darkgreen")
    tem <- subset(out1, trait == "BLUE_EH")
    plot(x=tem$present, y=tem$r, type="l", ylim=c(-0.3, 0.25), lwd=3, col="red",
         xlab="Cutoff for Present (FPKM)", ylab="Correlation (r)", main="Sum of GERP (Absent: FPKM =0.01)")
    for(i in 1:4){
        tem <- subset(out1, trait == traits[i])
        lines(x=tem$present, y=tem$r, type="l", lwd=3, col=cols[i])
    }
    legend("bottomright", legend=gsub("BLUE_", "", traits), lwd=3, col=c(cols, "red"))
    
    tem <- subset(out1, trait == "BLUE_EH")
    plot(x=tem$present, y=-log(tem$pval), type="l", lwd=3, col="red", ylim=c(0, 5),
         xlab="Cutoff for Present (FPKM)", ylab="-log(pval)", main="Sum of GERP (Absent: FPKM =0.01)")
    for(i in 1:4){
        tem <- subset(out1, trait == traits[i])
        lines(x=tem$present, y=-log(tem$pval), type="l", lwd=3, col=cols[i])
    }
    abline(h=-log(0.05), lty=2, lwd=2, col="grey")
    abline(v=69, lty=2, lwd=2, col="grey")
}

####################
out1 <- read.csv("largedata/1.gc/corout1_ab01_sum.csv")
parplot1(out1=out1)

parplot2 <- function(out1=out2){
    par(mfrow=c(1, 2))
    traits <- unique(out1$trait)
    cols <- c("blue", "black", "wheat", "darkgreen")
    tem <- subset(out1, trait == "BLUE_EH")
    plot(x=tem$present, y=tem$r, type="l", ylim=c(-0.3, 0.25), lwd=3, col="red",
         xlab="Cutoff for Present (FPKM)", ylab="Correlation (r)", main="Mean of GERP (Absent: FPKM =0.01)")
    for(i in 1:4){
        tem <- subset(out1, trait == traits[i])
        lines(x=tem$present, y=tem$r, type="l", lwd=3, col=cols[i])
    }
    legend("bottomright", legend=gsub("BLUE_", "", traits), lwd=3, col=c(cols, "red"))
    
    tem <- subset(out1, trait == "BLUE_EH")
    plot(x=tem$present, y=-log(tem$pval), type="l", lwd=3, col="red", ylim=c(0, 5),
         xlab="Cutoff for Present (FPKM)", ylab="-log(pval)", main="Mean of GERP (Absent: FPKM =0.01)")
    for(i in 1:4){
        tem <- subset(out1, trait == traits[i])
        lines(x=tem$present, y=-log(tem$pval), type="l", lwd=3, col=cols[i])
    }
    abline(h=-log(0.05), lty=2, lwd=2, col="grey")
    abline(v=20, lty=2, lwd=2, col="grey")
}

####################
out2 <- read.csv("largedata/1.gc/corout2_ab01_avg.csv")
parplot2(out1=out2)


parplot3 <- function(out=out3){
    par(mfrow=c(1, 2))
    traits <- unique(out$trait)
    cols <- c("blue", "black", "wheat", "darkgreen")
    tem <- subset(out, trait == "BLUE_EH")
    plot(x= -log2(tem$absent), y= tem$r, type="l", lwd=3, col="red", ylim=c(-0.1, 0.25),
         xlab="Cutoff for Absent (-log2(FPKM))", ylab="Correlation (r)", main="Sum of GERP (Present: FPKM =69)")
    for(i in 1:4){
        tem <- subset(out, trait == traits[i])
        lines(x=-log2(tem$absent), y=tem$r, type="l", lwd=3, col=cols[i])
    }
    abline(v=-log2(0.01), lty=2, col="red")
    legend("bottomright", legend=gsub("BLUE_", "", traits), lwd=3, col=c(cols, "red"))
    
    tem <- subset(out, trait == "BLUE_EH")
    plot(x= -log2(tem$absent), y= -log(tem$pval), type="l", lwd=3, col="red", ylim=c(0, 6),
         xlab="Cutoff for Present (-log2(FPKM))", ylab="-log(pval)", main="Sum of GERP (Present: FPKM =69)")
    for(i in 1:4){
        tem <- subset(out, trait == traits[i])
        lines(x=-log2(tem$absent), y=-log(tem$pval), type="l", lwd=3, col=cols[i])
    }
    abline(h=-log(0.05), lty=2, lwd=2, col="red")
    
}


out3 <- read.csv("largedata/1.gc/corout3_pre69_sum.csv")
parplot34(out=out3)

parplot4 <- function(out=out3){
    par(mfrow=c(1, 2))
    traits <- unique(out$trait)
    cols <- c("blue", "black", "wheat", "darkgreen")
    tem <- subset(out, trait == "BLUE_EH")
    plot(x= -log2(tem$absent), y= tem$r, type="l", lwd=3, col="red", ylim=c(-0.3, 0.2),
         xlab="Cutoff for Absent (-log2(FPKM))", ylab="Correlation (r)", main="Mean of GERP (Present: FPKM=20)")
    for(i in 1:4){
        tem <- subset(out, trait == traits[i])
        lines(x=-log2(tem$absent), y=tem$r, type="l", lwd=3, col=cols[i])
    }
    abline(v=-log2(0.01), lty=2, col="red")
    legend("bottomright", legend=gsub("BLUE_", "", traits), lwd=3, col=c(cols, "red"))
    
    tem <- subset(out, trait == "BLUE_EH")
    plot(x= -log2(tem$absent), y= -log(tem$pval), type="l", lwd=3, col="red", ylim=c(0, 3.5),
         xlab="Cutoff for Present (-log2(FPKM))", ylab="-log(pval)", main="Mean of GERP (Present: FPKM=20)")
    for(i in 1:4){
        tem <- subset(out, trait == traits[i])
        lines(x=-log2(tem$absent), y=-log(tem$pval), type="l", lwd=3, col=cols[i])
    }
    abline(h=-log(0.05), lty=2, lwd=2, col="red")
    
}


out4 <- read.csv("largedata/1.gc/corout4_pre20_avg.csv")
parplot4(out=out4)




#### Ames Hybrids
hyb <- read.table("data/Ames_hybrids_BLUEs.txt", header=TRUE)
names(hyb)[1] <- "F1"
hyb$P1 <- gsub("x.*", "", hyb$F1)
hyb$P2 <- gsub(".*x", "", hyb$F1)

hyb$P2 <- toupper(gsub(" ", "", hyb$P2))




###
res <- read.csv("cache/exp_res_p30a01.csv")
res$comp <- res$comp1gerp + res$comp2gerp

tem1 <- merge(hyb, res[, c("accenumb", "comp")], by.x="P2", by.y="accenumb")
tem2 <- merge(hyb, res[, c("geno", "comp")], by.x="P2", by.y="geno")
pheno <- rbind(tem1, tem2)
pheno <- subset(pheno, P1 == "PI601322")

pheno[pheno[,3]<0, ][,3] <- "NA"

####### original
for(i in 3:7){
    print(cor.test(log2(pheno$comp), as.numeric(as.character(pheno[,i]))))
    plot(pheno$comp, as.numeric(as.character(pheno[,i]))) 
}


### remove outlier:
pheno <- subset(pheno, comp < 100 & comp > 50)

par(mfrow=c(2,3))
for(i in 3:7){
    print(cor.test(log2(pheno$comp, as.numeric(as.character(pheno[,i]))))
    plot(pheno$comp, as.numeric(as.character(pheno[,i])), xlab="Number of complementary genes",
         ylab="Phenotypic value", main=names(pheno)[i], col="blue") 
}

0.2







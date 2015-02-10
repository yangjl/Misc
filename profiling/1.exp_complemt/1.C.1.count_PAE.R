### Jinliang Yang
### Jan. 31th, 2015
### PHZ51


#####################
countexp <- function(exp=exp0, xpid=xpid, present=10, absent=0.1){
    
    ### absent only in the PHZ51
    ### exp <- subset(exp, PHZ51 < 0.1)
    #### 1. absent in PHZ51, but present in others
    res1 <- data.frame()
    exp1 <- subset(exp, PHZ51 < absent)
    for(i in 3:ncol(exp1)){
        temexp <- subset(exp1, exp1[,i] > present)
        pid <- names(exp1)[i]
        temres <- data.frame(plantid=pid, comp1num=nrow(temexp),
                             comp1gerp=sum(temexp$gerpsum))
        res1 <- rbind(res1, temres)
    }
    
    #### 2. present in PHZ51, but absent in others
    res2 <- data.frame()
    exp2 <- subset(exp, PHZ51 > present)
    for(i in 3:ncol(exp2)){
        temexp <- subset(exp2, exp2[,i] < absent)
        pid <- names(exp2)[i]
        temres <- data.frame(plantid=pid, comp2num=nrow(temexp), 
                             comp2gerp=sum(temexp$gerpsum))
        res2 <- rbind(res2, temres)
    }
    
    res <- merge(res1, res2, by="plantid")
    res <- merge(xpid, res, by.x="geno", by.y="plantid")
    return(res)
    
}


### get the expression data
xpid <- read.csv("data/rnaseq_genoid.csv")
exp0 <- read.csv("largedata/1.gc/maize_gene_503lines.csv")
exp <- exp0[, as.character(xpid$geno)]
##########
gerp <- read.table("largedata/1.gc/sumgerp_in_gene.txt", header=FALSE)
names(gerp) <- c("chr", "start", "end", "geneid","gerpsum")

exp2 <- merge(gerp[, c("geneid", "gerpsum")],exp, by.x="geneid", by.y="row.names")



#### for PAV genes and FGS 
res <- countexp(exp=exp2, xpid=xpid, present=10, absent=0.1)
write.table(res, "cache/exp_res.csv", row.names=FALSE, quote=FALSE, sep=",")





res <- countexp(exp=exp, present=5, absent=0.1)
res <- merge(xpid, res, by.x="geno", by.y="plantid")
write.table(res, "cache/exp_res_p5a0.1.csv", row.names=FALSE, quote=FALSE, sep=",")

res <- countexp(exp=exp, present=50, absent=0.1)
res <- merge(xpid, res, by.x="geno", by.y="plantid")
write.table(res, "cache/exp_res_p50a0.1.csv", row.names=FALSE, quote=FALSE, sep=",")

res <- countexp(exp=exp, present=50, absent=1)
res <- merge(xpid, res, by.x="geno", by.y="plantid")
write.table(res, "cache/exp_res_p50a1.csv", row.names=FALSE, quote=FALSE, sep=",")

res <- countexp(exp=exp, present=10, absent=1)
res <- merge(xpid, res, by.x="geno", by.y="plantid")
write.table(res, "cache/exp_res_p10a1.csv", row.names=FALSE, quote=FALSE, sep=",")


#### FGS only 
idx <- grep("^joint", row.names(exp))
exp2 <- exp[-idx, ]
res2 <- countexp(exp=exp2, present=10, absent=0.1)
res2 <- merge(xpid, res2, by.x="geno", by.y="plantid")
write.table(res2, "cache/exp_res_FGS.csv", row.names=FALSE, quote=FALSE, sep=",")






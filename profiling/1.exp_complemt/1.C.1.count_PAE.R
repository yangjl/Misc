### Jinliang Yang
### Jan. 31th, 2015
### PHZ51

### get the expression data
getexpdata <- function(){
    xpid <- read.csv("data/rnaseq_genoid.csv")
    exp0 <- read.csv("largedata/1.gc/maize_gene_503lines.csv")
    exp <- exp0[, as.character(xpid$geno)]
    ##########
    gerp <- read.table("largedata/1.gc/sumgerp_in_gene.txt", header=FALSE)
    names(gerp) <- c("chr", "start", "end", "geneid","gerpsum")
    gerp$gerpavg <- gerp$gerpsum/(gerp$end - gerp$start +1)
    
    exp2 <- merge(gerp[, c("geneid","gerpsum", "gerpavg")], exp, by.x="geneid", by.y="row.names")
    return(exp2)
}


#####################
countexp <- function(exp=exp0, present=10, absent=0.1, method="sum"){
    
    ### absent only in the PHZ51
    ### exp <- subset(exp, PHZ51 < 0.1)
    #### 1. absent in PHZ51, but present in others
    res1 <- data.frame()
    exp1 <- subset(exp, PHZ51 < absent)
    for(i in 4:ncol(exp1)){
        temexp <- subset(exp1, exp1[,i] > present)
        pid <- names(exp1)[i]
        if(method == "sum"){
            temres <- data.frame(plantid=pid, comp1num=nrow(temexp),
                                 comp1gerp=sum(temexp$gerpsum))
        }
        if(method == "avg"){
            temres <- data.frame(plantid=pid, comp1num=nrow(temexp),
                                 comp1gerp=sum(temexp$gerpavg))
        }
        
        res1 <- rbind(res1, temres)
    }
    
    #### 2. present in PHZ51, but absent in others
    res2 <- data.frame()
    exp2 <- subset(exp, PHZ51 > present)
    for(i in 4:ncol(exp2)){
        temexp <- subset(exp2, exp2[,i] < absent)
        pid <- names(exp2)[i]
        if(method == "sum"){
            temres <- data.frame(plantid=pid, comp2num=nrow(temexp), 
                                 comp2gerp=sum(temexp$gerpsum))
        }
        if(method == "avg"){
            temres <- data.frame(plantid=pid, comp2num=nrow(temexp), 
                                 comp2gerp=sum(temexp$gerpavg))
        }
        res2 <- rbind(res2, temres)
    }
    
    ##########
    xpid <- read.csv("data/rnaseq_genoid.csv")
    res <- merge(res1, res2, by="plantid")
    res <- merge(xpid, res, by.x="geno", by.y="plantid")
    return(res)
    
}



### get the phenotypic data
getpheno <- function(res=res){
    #### Ames Hybrids
    hyb <- read.table("data/Ames_hybrids_BLUEs.txt", header=TRUE)
    names(hyb)[1] <- "F1"
    hyb$P1 <- gsub("x.*", "", hyb$F1)
    hyb$P2 <- gsub(".*x", "", hyb$F1)
    
    hyb$P2 <- toupper(gsub(" ", "", hyb$P2))
    hyb[hyb == -999] <- NA
    
    ##### single parameter checking
    #res <- read.csv("cache/exp_res.csv")
    res$comp <- res$comp1gerp + res$comp2gerp
    tem1 <- merge(hyb, res[, c("accenumb", "comp")], by.x="P2", by.y="accenumb")
    tem2 <- merge(hyb, res[, c("geno", "comp")], by.x="P2", by.y="geno")
    pheno <- rbind(tem1, tem2)
    pheno <- subset(pheno, P1 == "PI601322")
    #pheno$comp <- log2(pheno$comp)
    return(pheno)
}

##### Get correlation for 7 traits with one set of parameters
getcor <- function(pheno=pheno, method = "sum"){
    if(method == "sum"){
        x <- log2(pheno$comp) 
    }else if(method == "avg"){
        x <- pheno$comp
    }else{
        stop("No such method, => 'sum' or 'avg' ")
    }
    
    myres2 <- data.frame()
    for(i in 3:7){
        y <- as.numeric(as.character(pheno[,i]))
        df <- data.frame(x=x, y=y)
        df <- subset(df, !is.na(x) & !is.na(y))
        pval <- cor.test(df$x, df$y)$p.value
        
        print(sprintf("###>>> trait= [ %s ],  pval = [ %s ]", names(pheno)[i], pval))
        tem <- data.frame(trait=names(pheno)[i], pval=pval, r=cor(df$x, df$y))
        myres2 <- rbind(myres2, tem)    
        
    }
    return(myres2) 
}


### Jinliang Yang
### 01-22-2017
### study NAM and 282 phenotypes

library("tidyr")
gatherdata <- function(){
    
    ### trait mx from Panzea
    pheno <- read.delim("data/traitMatrix_maize282NAM_v15-130212.txt", header=TRUE)
    names(pheno)[1] <- "genotype"
    pheno[2,1] <- "genotype"
    
    nm <- gsub("\\..$|\\...$", "", names(pheno))
    nmtb <- data.frame(id=names(pheno), uid=nm)
    tid <- unique(nm[-1])
    #length(tid) #58 trait id
    
    message(sprintf("###>>> loaded [ %s ] traits from matrix ...", length(tid)))
    
    df <- data.frame()
    for(i in 1:length(tid)){
        idx <- subset(nmtb, uid %in% tid[i])
        
        tpheno <- pheno[, c("genotype", as.character(idx$id))]
        newid <- as.character(t(tpheno[1, ]))
        tpheno <- tpheno[-1, ]
        names(tpheno)[2:ncol(tpheno)] <- newid[-1]
        
        tem <- gather(tpheno, 2:ncol(tpheno), key="location",  value="value", na.rm=TRUE)
        tem$trait <- tid[i]
        tem$value <- as.numeric(as.character(tem$value))
        tem <- tem[!is.na(tem$value), ]
        df <- rbind(df, tem)
        print(i)
    }
    
    return(df)
}

########
df <- gatherdata()
df <- df[!is.na(df$genotype),]
write.table(df, "largedata/pheno_long.csv", sep=",", row.names=FALSE, quote=FALSE)


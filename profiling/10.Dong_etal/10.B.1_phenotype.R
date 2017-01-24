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


cw <- subset(df, trait %in% "CobWeight")
outlier <- tail(sort(cw$value), n=30)
cw <- cw[cw$value < 100, ]


hist(cw$value)


source("~/Documents/Github/zmSNPtools/Rcodes/mixed_model.R")

    
getpheno <- function(traits="X10KW_Inbred"){
    nathan <- read.csv("data/Nathan_p1_subset.csv")
    a1 <- nrow(nathan)
    nathan$INBRED <- toupper(nathan$INBRED)
    #genotype <- unique(nathan$INBRED)
    #nathan <- read.csv("data/Nathan_p1_subset.csv")
    
    toi <- traits[1]
    nathan <- nathan[nathan[, toi] !=".",]
    nathan[, toi] <- as.numeric(as.character(nathan[, toi]))
    message(sprintf("###>>> total [ %s ], remaining [ %s ] for trait [ %s ]", a1, nrow(nathan), toi))
    out <- mixed_model(data = nathan, model = as.formula(paste(toi, "~ INBRED")), random = ~1 | Env, 
                          trait = toi)
    out$sex <-1
    out <- out[, c(1,3,2)]
    for(i in 2:length(traits)){
        toi <- traits[i]
        nathan <- nathan[nathan[, toi] !=".",]
        nathan[, toi] <- as.numeric(as.character(nathan[, toi]))
        message(sprintf("total [ %s ], remaining [ %s ] for trait [ %s ]", a1, nrow(nathan), toi))
        myblue <- mixed_model(data = nathan, model = as.formula(paste(toi, "~ INBRED")), random = ~1 | Env, 
                              trait = toi)
        out <- merge(out, myblue, by="Genotype", all=TRUE)
    }
    
    nm <- names(out)
    nm <- gsub("_.*", "", nm)
    names(out) <- nm
    return(out)
}

##########
mytraits <- c("X10KW_Inbred","KC_Inbred","CD_Inbred","CL_Inbred","CW_Inbred","TKW_Inbred")
pheno <- getpheno(traits=mytraits)
names(pheno)[1] <- "id"

write.table(pheno, "data/pheno_ames282.txt", row.names=FALSE, quote=FALSE, sep="\t")

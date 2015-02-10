### Jinliang Yang
### Feb. 10th, 2015
### loop through different parameters

source("profiling/1.exp_complemt/1.C.1.count_PAE.R")

###### loop to selet the parameters
loopP <- function(absent_cutoff = 0.01, method="sum"){
    
    myres <- data.frame()
    abrange <- seq(1, 100, by=1)
    expdata <- getexpdata()
    for(k in abrange){
        res <- countexp(exp=expdata, present= k, absent= absent_cutoff, method=method)
        pheno <- getpheno(res=res)
        out <- getcor(pheno=pheno, method = method)
        out$present <- k
        out$absent <- absent_cutoff
        myres <- rbind(myres, out)
    }
    return(myres)
}

###### loop to selet the parameters
loopA <- function(present_cutoff = 69, method="sum"){
    
    myres <- data.frame()
    abrange <- seq(0.001, 0.1, by=0.001)
    expdata <- getexpdata()
    for(k in abrange){
        res <- countexp(exp=expdata, present= present_cutoff, absent= k, method=method)
        pheno <- getpheno(res=res)
        out <- getcor(pheno=pheno, method = method)
        out$present <- present_cutoff
        out$absent <- k
        myres <- rbind(myres, out)
    }
    return(myres)
}
##################
out1 <- loopP(absent_cutoff = 0.01, method="sum")
write.table(out1, "largedata/1.gc/corout1_ab01_sum.csv", sep=",", row.names=FALSE, quote=FALSE)

out2 <- loopP(absent_cutoff = 0.01, method="avg")
write.table(out2, "largedata/1.gc/corout2_ab01_avg.csv", sep=",", row.names=FALSE, quote=FALSE)

out3 <- loopA(present_cutoff = 69, method="sum")
write.table(out3, "largedata/1.gc/corout3_pre69_sum.csv", sep=",", row.names=FALSE, quote=FALSE)

out4 <- loopA(present_cutoff = 20, method="avg")
write.table(out4, "largedata/1.gc/corout4_pre20_avg.csv", sep=",", row.names=FALSE, quote=FALSE)



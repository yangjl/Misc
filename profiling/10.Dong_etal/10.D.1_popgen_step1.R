### Jinliang Yang
### 02-01-2017
### purpose: pop-gen analysis




taxa <- read.csv("data/hmp3_taxa.csv")
table(taxa$Dataset)

line <- subset(taxa, !(Dataset %in% "German") & !(Dataset %in% "CIMMYT/BGI"))
length(unique(line$Prefixed_Taxon))
#1098

imp <- subset(line, Category %in% " Improved") #
land <- subset(line, Category %in% "Landrace" | Category %in% "andrace") #23
teo <- subset(line, Category %in% "Parviglumis") #27

# GRMZM2G039867 tru1
#BP <- 1000000
#chr3snp <- as.character(subset(tb, chr == 3 & pos > 150088574 - BP & pos < 150091550 + BP)$snpid)
# length(chr3snp)

library(PopGenome)

vcf_handle <- .Call("VCF_open",filename="/home/jolyang/dbcenter/HapMap/HapMap3/merged_flt_ad_c1.vcf.gz")
ind <- .Call("VCF_getSampleNames",vcf_handle)
#sid <- ind[1:1000]
d282 <- ind[grep("282set", ind)]

BP <- 1000000
gdat <- readVCF("/home/jolyang/dbcenter/HapMap/HapMap3/merged_flt_ad_c3.vcf.gz", numcols=10000, tid="3",
                samplenames=c(d282, as.character(land$Prefixed_Taxon), as.character(teo$Prefixed_Taxon)),  
                include.unknown=TRUE,
                frompos = 150088574 - BP, topos = 150091550 + BP, approx=FALSE, out="", 
                parallel=FALSE, gffpath=FALSE)


gdat@n.sites
gdat@region.names
gdat@region.data

get.sum.data(gdat)


### assign populations
gdat <- set.populations(gdat, list(d282, 
                                   as.character(land$Prefixed_Taxon), 
                                   as.character(teo$Prefixed_Taxon)), 
                        diploid=TRUE)


# Sliding window analyses
slide <- sliding.window.transform(gdat, 10000, 10000, type=2)

# total number of windows
length(slide@region.names) #20020

# Statistics
slide <- diversity.stats(slide)
nucdiv <- slide@nuc.diversity.within
# the values have to be normalized by the number of nucleotides in each window
nucdiv <- nucdiv/200
head(nucdiv)

write.csv(nucdiv, "largedata/pai_tru1.csv", row.names=FALSE)

nucdiv <- read.csv("cache/pai_tru1.csv")
nucdiv <- nucdiv[, -1]
# Generate output
# Smoothing lines via spline interpolation
ids <- 1:20020
s = 0.3
loess.nucdiv1 <- loess(nucdiv[,1] ~ ids, span=s)
loess.nucdiv2 <- loess(nucdiv[,2] ~ ids, span=s)
loess.nucdiv3 <- loess(nucdiv[,3] ~ ids, span=s)
plot(predict(loess.nucdiv1), type = "l", xaxt="n", xlab="position (Mb)",
     ylab="nucleotide diversity", main = "Chromosome 2L (10kb windows)", ylim=c(0.0001, 0.0004))
lines(predict(loess.nucdiv2), col="blue")
lines(predict(loess.nucdiv3), col="red")

abline(v=10010)

axis(1,c(1,1000,2000,3000,4000,5000),
     c("0","10","20","30","40","50"))
# create the legend
legend("topright",c("M","S","X"),col=c("black","blue","red"), lty=c(1,1,1))


nucdiv <- read.csv("largedata/pai_tru1.csv")
# Generate output
# Smoothing lines via spline interpolation
ids <- 1:200
s = 0.3
loess.nucdiv1 <- loess(nucdiv[,1] ~ ids, span=s)
loess.nucdiv2 <- loess(nucdiv[,2] ~ ids, span=s)
loess.nucdiv3 <- loess(nucdiv[,3] ~ ids, span=s)
plot(predict(loess.nucdiv1), type = "l", xaxt="n", xlab="position (Mb)",
     ylab="nucleotide diversity", main = "Chromosome 2L (10kb windows)", ylim=c(0.1, 0.3))
lines(predict(loess.nucdiv2), col="blue")
lines(predict(loess.nucdiv3), col="red")

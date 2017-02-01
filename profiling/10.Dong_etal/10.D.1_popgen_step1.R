### Jinliang Yang




line <- read.table("~/dbcenter/HapMap/HapMap3/vcf.header")
newline <- toupper(as.vector(t(line)))
d282 <- grep("282SET", newline)

taxa <- read.csv("data/hmp3_taxa.csv")
table(taxa$Dataset)

line <- subset(taxa, !(Dataset %in% "German") & !(Dataset %in% "CIMMYT/BGI"))
length(unique(line$Prefixed_Taxon))
#1098

imp <- subset(line, Category %in% " Improved") #
land <- subset(line, Category %in% "Landrace" | Category %in% "andrace") #23
teo <- subset(line, Category %in% "Parviglumis") #27

#d282 <- subset(line, )

library(PopGenome)

vcf_handle <- .Call("VCF_open",filename="/home/jolyang/dbcenter/HapMap/HapMap3/merged_flt_ad_c1.vcf.gz")
ind <- .Call("VCF_getSampleNames",vcf_handle)
#sid <- ind[1:1000]
d282 <- ind[grep("282set", ind)]

GENOME2 <- readVCF("/home/jolyang/dbcenter/HapMap/HapMap3/merged_flt_ad_c10.vcf.gz", numcols=10000, tid="10",
                   samplenames=d282,  include.unknown=TRUE,
                   frompos = 1, topos = 1000000, approx=FALSE, out="", parallel=FALSE, gffpath=FALSE)

GENOME.class@n.sites
GENOME.class@region.names
GENOME.class@region.data

get.sum.data(GENOME2)








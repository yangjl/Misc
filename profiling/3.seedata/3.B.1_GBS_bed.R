### Jinliang Yang
### 4/7/2015
### transform GBS format to BED+ format with haplotype call

seeds <- read.delim("data/seeds_09.02.2015_22.38.10.txt")
#[1] 22022    51
subseed <- subset(seeds, general_identifier %in% idtab$GID)

idtab <- read.csv("data/SeeD_SID_to_GID.csv")
length(unique(idtab$GID))
length(unique(idtab$SampleID))
# Note: some of the accessions were genotyped multiple times


##### transform GBS to BED format
library(data.table)
source("lib/gbs2bed.R")


gbs2plink(gbsfile="largedata/genotypes/SB.imputed.hmp",
        outfile="largedata/genotypes/SB.imputed")


system("plink --tfile SB.imputed --make-bed --out SB.imputed")
system("plink --bfile SB.imputed --freq --missing --out SB.imputed")
system("plink --bfile SB.imputed --het --ibc --out SB.imputed")

system("plink --bfile SB.imputed --indep-pairwise 100 10 0.05 --r2 --out SB.imputed")






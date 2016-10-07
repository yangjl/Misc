### Jinliang Yang
### Feb. 9th, 2015

####### using bedtools to compute the GERP score in genes
cd largedata/1.gc/
bedtools intersect -a FGSv2.bed4 -b gerpv2.bed4 -wa -wb > gerpToFGS_v2.bed
bedtools groupby -i gerpToFGS_v2.bed -g 1,2,3,4 -c 8 -o sum > sumgerp_in_gene.txt


##########
gerp <- read.table("largedata/1.gc/sumgerp_in_gene.txt", header=FALSE)
names(gerp) <- c("chr", "start", "end", "geneid","gerpsum")
gerp$gerpavg <- gerp$gerpsum/(gerp$end - gerp$start +1)

par(mfrow=c(1,2))
hist(log2(gerp$gerpsum), breaks=30, xlab="Log2 (sum of GERP in genes)", col="wheat",
     main="Sum of GERP score in FGSv2")
hist(gerp$gerpavg, breaks=30, xlab="mean GERP in genes", col="wheat",
     main="Mean GERP score in FGSv2")

###########
exp0 <- read.csv("largedata/1.gc/maize_gene_503lines.csv")

res <- apply(exp0[, 4:506], 1, mean)
xp <- data.frame(geneid=names(res), exp=res)

res2 <- merge(gerp, xp, by="geneid")
res2 <- subset(res2, exp > 1/500 & exp < 200)

res2$logexp <- log2(res2$exp)
#res2 <- subset(res2, abs(logexp) > 10)
plot(res2$logexp, res2$gerpavg)


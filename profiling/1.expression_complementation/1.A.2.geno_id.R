### Jinliang Yang
### Jan. 31th, 2015
### purpose: find the useful genotypes

#### GRIN lookuptable
grin <- read.csv("largedata/1.gc/GRIN.csv", header=TRUE)
grintab <- grin[, c("accenumb", "accename")]

### Note: plantid: uppercase, no space, special characters([ ._]()) => ""
grintab$plantid <- grintab$accename
grintab$accenumb <- toupper(gsub(" ", "", grintab$accenumb))
grintab$plantid <- toupper(gsub(" ", "", grintab$plantid))
grintab[grintab$plantid=="", ]$plantid <- grintab[grintab$plantid=="", ]$accenumb

#test <- gsub("[\\.\\(\\)\\_]", "", "a._b-)c)(d.")
grintab$plantid <- gsub("[\\.\\(\\)\\_]", "", grintab$plantid)


#############################
geno <- read.csv("largedata/1.gc/maize_gene_503lines.csv", nrows=5)
genodf <- data.frame(geno=names(geno)[-1:-3], type="xp")
genodf$geno2 <- toupper(gsub(" ", "", genodf$geno))
genodf$geno2 <- gsub("[\\.\\(\\)\\_]", "", genodf$geno2)

geno <- merge(genodf, grintab, by.x="geno2", by.y="plantid")
length(unique(geno$geno)) #426




#### Ames Hybrids
hyb <- read.table("data/Ames_hybrids_BLUEs.txt", header=TRUE)
names(hyb)[1] <- "F1"
hyb$P1 <- gsub("x.*", "", hyb$F1)
hyb$P2 <- gsub(".*x", "", hyb$F1)

hyb$P2 <- toupper(gsub(" ", "", hyb$P2))

length(unique(c(hyb$P1, hyb$P2)))

####PHZ51 in the transcriptomic dataset
subset(geno, accenumb %in% unique(hyb$P1))

sub0 <- subset(geno, accenumb %in% unique(hyb$P1)) #1 PHZ51
sub1 <- subset(geno, geno %in% hyb$P2) #39
sub2 <- subset(geno, accenumb %in% hyb$P2) #191

xpgeno <- rbind(sub0, sub1, sub2)
xpgeno <- xpgeno[!duplicated(xpgeno$geno),]
write.table(xpgeno[, -3], "data/rnaseq_genoid.csv", sep=",", row.names=FALSE, quote=FALSE)



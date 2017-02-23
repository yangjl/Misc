### Jinliang Yang
### Feb 22th, 2017


library(data.table)
d <- fread("largedata/Hufford2012/757736/Hufford_et_al._2012_genestats_MZ.txt")
t <- fread("largedata/Hufford2012/757736/Hufford_et_al._2012_genestats_TEO.txt")
tot <- nrow(d)

d <- subset(d, !is.na(ThetaPi))
mean(d$ThetaPi/d$seqbp, na.rm=T) # 0.007215115
pai <- 0.00128233
sum(d$ThetaPi/d$seqbp < pai)/tot

mean(d$TajD, na.rm=T)


t <- subset(t, !is.na(ThetaPi))
mean(t$ThetaPi/t$seqbp, na.rm=T) # 0.007801881
mean(t$TajD, na.rm=T)


hist(log(d$ThetaPi/d$seqbp), breaks=50)
idx <- grep("GRMZM2G039867", d$locus)
tru1 <- d[idx,]
abline(v=log(tru1$ThetaPi/tru1$seqbp))
abline(v=log(0.001318671))





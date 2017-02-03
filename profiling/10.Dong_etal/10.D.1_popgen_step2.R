### Jinliang Yang
### 02-01-2017
### purpose: pop-gen analysis

slide <- F_ST.stats(slide, mode="nucleotide")
pairwise.FST <- t(slide@nuc.F_ST.pairwise)
head(slide@nuc.F_ST.vs.all)


write.csv(slide@nuc.F_ST.vs.all, "largedata/Fst_tru1_win100.csv", row.names=FALSE)
write.csv(pairwise.FST, "largedata/pw_Fst_tru1_win1k.csv", row.names=FALSE)

f <- read.csv("largedata/Fst_tru1_win100.csv")
# Generate output
# Smoothing lines via spline interpolation
ids <- 1:nrow(f)
s = 0.1
loess.nucdiv1 <- loess(f[,1] ~ ids, span=s)
loess.nucdiv2 <- loess(f[,2] ~ ids, span=s)
loess.nucdiv3 <- loess(f[,3] ~ ids, span=s)
plot(predict(loess.nucdiv1), type = "l", xaxt="n", xlab="position (Mb)",
     ylab="nucleotide diversity", main = "Chromosome 2L (10kb windows)")
lines(predict(loess.nucdiv2), col="blue")
lines(predict(loess.nucdiv3), col="red")

abline(v=1000)


#########################


slide <- neutrality.stats(slide, FAST=TRUE)
res <- get.neutrality(slide)
write.csv(res[[1]], "largedata/neu_d282.csv", row.names=FALSE)
write.csv(res[[3]], "largedata/neu_teo.csv", row.names=FALSE)


neu <- read.csv("largedata/neu_d282.csv")
teo <- read.csv("largedata/neu_teo.csv")
#plot(1:nrow(neu), neu[,1])
# Generate output
# Smoothing lines via spline interpolation
ids <- 1:nrow(neu)
s = 0.2

loess.nucdiv1 <- loess(neu[,1] ~ ids, span=s)
loess.nucdiv2 <- loess(teo[,1] ~ ids, span=s)
plot(predict(loess.nucdiv1), type = "l", xaxt="n", xlab="position (Mb)",
     ylab="Tajima's D", main = "Chromosome 2L (10kb windows)", lwd=2)
lines(predict(loess.nucdiv2), col="red", lwd=2, lty=2)
legend("topright",c("M","S","X"),col=c("black","blue","red"), lty=c(1,1,1))
abline(v=10000)

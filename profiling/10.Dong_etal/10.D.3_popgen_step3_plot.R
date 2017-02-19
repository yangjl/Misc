## Jinliang Yang
## 02-02-2017
## plot popgen analysis results

BP <- 5000
start <- 151329451 - BP
end <- 151332389 + BP

neu <- read.csv("largedata/neu_d282.csv")
teo <- read.csv("largedata/neu_teo.csv")
ids <- 1:nrow(neu)
s = 0.2



pdf("graphs/Fig3.pdf",  width=16, height=10)
par(mfrow=c(2,1))
lo1 <- loess(neu[,1] ~ ids, span=s)
lo2 <- loess(teo[,1] ~ ids, span=s)
plot(predict(lo1), type = "l", xaxt="n", xlab="position (kb)",
     ylab="Tajima's D", main = "Chr 3 (win=100bp, step=25bp)", lwd=3, col="darkgreen")
lines(predict(lo2), col="red", lwd=3, lty=2)
abline(h=0, col="grey")
legend("topright", c("Maize","Teosinte"), col=c("darkgreen","red"), lty=c(1,2), lwd=c(3,3))
#axis(1,c(1,1000,2000,3000,4000,5000), c("0","10","20","30","40","50"))
abline(v=204, col="grey")
abline(v=514-200, col="grey")
axis(1,c(1,100,200,300,400,500), 
     round(c(start,start+25*100,start+25*200,start+25*300,start+25*400,start+25*500)/1000))


###############################
nucdiv <- read.csv("largedata/agpv4_pai_tru1_win100_25.csv")
# Generate output
# Smoothing lines via spline interpolation
ids <- 1:nrow(nucdiv)
s = 0.2
loess.nucdiv1 <- loess(nucdiv[,1] ~ ids, span=s)
#loess.nucdiv2 <- loess(nucdiv[,2] ~ ids, span=s)
loess.nucdiv3 <- loess(nucdiv[,3] ~ ids, span=s)
plot(predict(loess.nucdiv1), type = "l", xaxt="n", xlab="position (kb)", col="darkgreen",
     ylab="Nucleotide Diversity", main = "Chr 3 (win=100bp, step=25bp)", ylim=c(0, 0.06), lwd=3)
#lines(predict(loess.nucdiv2), col="blue", lwd=2)
lines(predict(loess.nucdiv3), col="red", lwd=3)
legend("topright", c("Maize","Teosinte"), col=c("darkgreen","red"), lty=c(1,2), lwd=c(3,3))
#axis(1,c(1,1000,2000,3000,4000,5000), c("0","10","20","30","40","50"))
abline(v=204, col="grey")
abline(v=514-200, col="grey")
axis(1,c(1,100,200,300,400,500), 
     round(c(start,start+25*100,start+25*200,start+25*300,start+25*400,start+25*500)/1000))

dev.off()

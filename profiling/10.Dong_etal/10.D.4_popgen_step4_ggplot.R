## Jinliang Yang
## 02-02-2017
## plot popgen analysis results

library(ggplot2)

BP <- 5000
start <- 151329451 - BP
end <- 151332389 + BP

neu <- read.csv("largedata/neu_d282.csv")
teo <- read.csv("largedata/neu_teo.csv")
ids <- 1:nrow(neu)
s = 0.2

neu$pop <- "Maize"
neu$x <- seq(from=start, to=end, length.out=nrow(neu))
teo$pop <- "Teosinte"
teo$x <- x <- seq(from=start, to=end, length.out=nrow(neu))
out <- rbind(neu, teo)



fsize=14

p1 <- ggplot(out, aes(x=x, y=Tajima.D, colour=factor(pop, levels=c("Maize", "Teosinte")) )) +
    labs(colour="") +
    theme_bw() +
    xlab("Chr3 (bp)") +
    ylab("Tajima'D") +
    scale_color_manual(values=c("#458b74", "#8b2323")) +
    #scale_linetype_manual(values=lty1) +
    guides(colour = guide_legend()) +
    geom_smooth(method="loess", span=0.5) +
    geom_hline(yintercept = 0, col="grey") +
    geom_vline(xintercept = c(151329451, 151332389), col="grey") +
    theme(axis.text.y = element_text(angle = 90, hjust = 0.5),
          axis.text=element_text(size=fsize),
          axis.title=element_text(size=fsize, face="bold"),
          legend.title = element_text(size=fsize, face="bold"),
          legend.text = element_text(size=fsize))
p1


library(tidyr)
nucdiv <- read.csv("largedata/agpv4_pai_tru1_win100_25.csv")
nucdiv$x <- seq(from=start, to=end, length.out=nrow(neu))
nucdiv <- nucdiv[, -2]
names(nucdiv)[1:2] <- c("Maize", "Teosinte")
out2 <- gather(nucdiv, pop, pai, 1:2)


p2 <- ggplot(out2, aes(x=x, y=pai, colour=factor(pop, levels=c("Maize", "Teosinte")) )) +
    labs(colour="") +
    theme_bw() +
    xlab("") +
    ylab("Nucleotide Diversity") +
    scale_color_manual(values=c("#458b74", "#8b2323")) +
    #scale_linetype_manual(values=lty1) +
    guides(colour = guide_legend()) +
    geom_smooth(method="loess", span=0.5) +
    geom_hline(yintercept = 0, col="grey") +
    geom_vline(xintercept = c(151329451, 151332389), col="grey") +
    theme(axis.text.y = element_text(angle = 90, hjust = 0.5),
          axis.text=element_text(size=fsize),
          axis.title=element_text(size=fsize, face="bold"),
          legend.title = element_text(size=fsize, face="bold"),
          legend.text = element_text(size=fsize))
p2

library(cowplot)

pdf("graphs/Fig3.pdf", height=6, width=8)
plot_grid(p2,p1, cols=1, labels=c("A", "B"))
dev.off()

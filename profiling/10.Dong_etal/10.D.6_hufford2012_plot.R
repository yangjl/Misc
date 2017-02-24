
library(data.table)
library(tidyr)
#AGPv2:GRMZM2G039867, chr3:150,049,274..150,052,250

d <- fread("largedata/Hufford2012/757736/Hufford_et_al._2012_10kb_statistics.txt", data.table=FALSE)

BP = 50000
subd <- subset(d, chr==3 & winstart > 150049274 -BP & winend < 150052250 + BP)
tem <- subd[, c("chr", "winstart", "winend", "ThetaPiMZ", "ThetaPiteo")]
out <- gather(tem, pop, pi, 4:5)


#LocusID L S n D startingtheta inheritancescalar
df <- subd[, c("chr", "winstart", "winend", "seqbp", "SMZ", "div_sites", "ThetaPiMZ")]
df <- subset(df, !is.na(div_sites))
df$lid <- paste0("Chr", df$chr, "g", df$winstart)
df$n <- 37
df$inheritancescalar <- 1

df <- df[, c("lid", "seqbp", "SMZ", "n", "div_sites", "ThetaPiMZ", "inheritancescalar")]

cat(c(nrow(df), 1, "a1", 23.2, "\n"),file="largedata/infile.txt", sep="\t")
write.table(df, "largedata/infile.txt", append = TRUE, row.names=FALSE, 
            quote=FALSE, sep="\t", col.names=FALSE)


mz <- fread("largedata/Hufford2012/757736/Hufford_et_al._2012_genestats_MZ.txt", data.table=FALSE)

p2 <- ggplot(out, aes(x=winstart, y=pi, colour=factor(pop, levels=c("ThetaPiMZ", "ThetaPiteo")) )) +
    labs(colour="") +
    theme_bw() +
    xlab("") +
    ylab("Nucleotide Diversity") +
    scale_color_manual(values=c("#458b74", "#8b2323")) +
    #scale_linetype_manual(values=lty1) +
    guides(colour = guide_legend()) +
    geom_smooth(method="loess", span=0.6) +
    geom_hline(yintercept = mean(d$ThetaPiMZ, na.rm=T), col="#458b74") +
    geom_hline(yintercept = mean(d$ThetaPiteo, na.rm=T), col="#8b2323") +
    geom_vline(xintercept = c(150049274, 150052250), col="grey") +
    theme(axis.text.y = element_text(angle = 90, hjust = 0.5),
          axis.text=element_text(size=fsize),
          axis.title=element_text(size=fsize, face="bold"),
          legend.title = element_text(size=fsize, face="bold"),
          legend.text = element_text(size=fsize))
p2

p1 <- ggplot(out, aes(x=x, y=Tajima.D, colour=factor(pop, levels=c("Maize", "Teosinte")) )) +
    labs(colour="") +
    theme_bw() +
    xlab("Chr3 (bp)") +
    ylab("Tajima'D") +
    scale_color_manual(values=c("#458b74", "#8b2323")) +
    #scale_linetype_manual(values=lty1) +
    guides(colour = guide_legend()) +
    geom_smooth(method="loess", span=0.5) +
    geom_hline(yintercept = 0.5334655, col="#458b74") +
    geom_hline(yintercept = 0.3942796, col="#8b2323") +
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



library(cowplot)

pdf("graphs/Fig3.pdf", height=6, width=8)
plot_grid(p2,p1, cols=1, labels=c("A", "B"))
dev.off()


2*(log(lh1) - log(lh2))

hist(rchisq(n=1000, df=2))

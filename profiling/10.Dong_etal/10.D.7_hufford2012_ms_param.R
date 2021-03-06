### Jinliang Yang
### 03-20-2017

library(data.table)
library(tidyr)
library(ggplot2)
#AGPv2:GRMZM2G039867, chr3:150,049,274..150,052,250
# AGPv4: 151329451 151332389 + BP
d <- fread("largedata/Hufford2012/757736/Hufford_et_al._2012_10kb_statistics.txt", data.table=FALSE)
fsize=12

BP = 50000
subd <- subset(d, chr==3 & winstart > 150049274 -BP & winend < 150052250 + BP)

#bed <- subd[, c(1, 9, 10)]
#bed$uid <- paste(bed[,1], bed[,2], sep="_")
#write.table(bed, "data/chr3_v2.bed", sep="\t", row.names=FALSE, quote=FALSE)

#nbed <- read.table("data/output_chr3_v2.bed")
#nbed <- nbed[!duplicated(nbed$V4),]

v24 <- 151329451 - 150049274

tem <- subd[, c("chr", "winstart", "winend", "ThetaPiMZ", "ThetaPiLR", "ThetaPiteo")]
out <- gather(tem, pop, pi, 4:6)
out <- subset(out, !is.na(pi))

out[out$pop == "ThetaPiMZ", ]$pop <- "Maize"
out[out$pop == "ThetaPiteo", ]$pop <- "Teosinte"
out[out$pop == "ThetaPiLR", ]$pop <- "Landrace"

out$winstart <- out$winstart + v24

p2 <- ggplot(out, aes(x=winstart, y=pi, colour=factor(pop, levels=c("Maize", "Landrace", "Teosinte")) )) +
    labs(colour="") +
    theme_bw() +
    xlab("") +
    ylab("Nucleotide Diversity") +
    scale_color_manual(values=c("#458b74", "#6495ed", "#8b2323")) +
    #scale_linetype_manual(values=lty1) +
    guides(colour = guide_legend()) +
    geom_smooth(method="loess", span=0.6) +
    geom_hline(yintercept = mean(d$ThetaPiMZ, na.rm=T), col="#458b74") +
    geom_hline(yintercept = mean(d$ThetaPiteo, na.rm=T), col="#8b2323") +
    geom_hline(yintercept = mean(d$ThetaPiLR, na.rm=T), col="#8b2323") +
    geom_vline(xintercept = c(150049274, 150052250) +v24, col="grey") +
    theme(axis.text.y = element_text(angle = 90, hjust = 0.5),
          axis.text=element_text(size=fsize),
          axis.title=element_text(size=fsize, face="bold"),
          legend.title = element_text(size=fsize, face="bold"),
          legend.text = element_text(size=fsize) )
          #axis.text.x= element_blank() ) +
          #axis.ticks.x=element_blank())
    

pdf("graphs/Fig4.pdf", height=4, width=6)
p2
dev.off()




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

#############
library(data.table)
library(tidyr)
library(ggplot2)
d <- fread("largedata/Hufford2012/757736/Hufford_et_al._2012_10kb_statistics.txt", data.table=FALSE)


BP = 50000
subd <- subset(d, chr==3 & winstart > 150049274 -BP & winend < 150052250 + BP)
tem <- subd[, c("chr", "winstart", "winend", "TajDMZ", "TajDLR", "TajDteo")]
out <- gather(tem, pop, tajd, 4:6)
out <- subset(out, !is.na(tajd))

out[out$pop == "TajDMZ", ]$pop <- "Maize"
out[out$pop == "TajDLR", ]$pop <- "Landrace"
out[out$pop == "TajDteo", ]$pop <- "Teosinte"

out$winstart <- out$winstart + v24

p1 <- ggplot(out, aes(x=winstart, y=tajd, colour=factor(pop, levels=c("Maize", "Landrace", "Teosinte")) )) +
    labs(colour="") +
    theme_bw() +
    xlab("") +
    ylab("Tajima'D") +
    scale_color_manual(values=c("#458b74", "#6495ed", "#8b2323")) +
    #scale_linetype_manual(values=lty1) +
    guides(colour = guide_legend()) +
    geom_smooth(method="loess", span= 1.5) +
    geom_hline(yintercept = mean(d$TajDMZ, na.rm=T), col="#458b74") +
    geom_hline(yintercept = mean(d$TajDLR, na.rm=T), col="#6495ed") +
    geom_hline(yintercept = mean(d$TajDteo, na.rm=T), col="#8b2323") +
    geom_vline(xintercept = c(150049274, 150052250) +v24, col="grey") +
    theme(axis.text.y = element_text(angle = 90, hjust = 0.5),
          axis.text=element_text(size=fsize),
          axis.title=element_text(size=fsize, face="bold"),
          legend.title = element_text(size=fsize, face="bold"),
          legend.text = element_text(size=fsize))
p1


library(cowplot)

pdf("graphs/Fig3.pdf", height=5, width=13, onefile=FALSE)
plot_grid(p2, p1, ncol=2, labels=c("A", "B"))
dev.off()


a0 <- read.table("largedata/outfile0.txt", header=T)
a2 <- read.table("largedata/outfile2.txt", header=T, skip = 1)
    
t1 <- 2*(log(a0[,9]) - log(a2[,9]))
t2 <- 2*abs(log(a0[,11]) - log(a2[,11]))

dchisq(t1, df=1)
dchisq(t2, df=1)

hist(rchisq(n=1000, df=2))

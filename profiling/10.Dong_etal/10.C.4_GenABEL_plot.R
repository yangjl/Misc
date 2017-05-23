### Jinliang Yang
### Feb, 17th, 2017
### purpose: plot GWAS results

library("GenomeGraphs")



#151329451..151332389
library("biomaRt")

listMarts(host="plants.ensembl.org")
ensembl <- useMart(biomart="plants_mart", host="plants.ensembl.org")
listDatasets(ensembl)

zm <- useMart(biomart="plants_mart", host="plants.ensembl.org", dataset="zmays_eg_gene")

# check attributes are available to select. More information on ensembl data base
out <- listAttributes(zm)

# check which filters are available
listFilters(zm)

gene <- makeGene(id = "ENSG00000095203", type="ensembl_gene_id", biomart = mart)

# fuction to get  gene id's and gene name from data base
gid <- getBM(attributes=c("ensembl_gene_id", "external_gene_id", "chromosome_name", "start_position", "end_position"), 
                       values="gene_names", mart= zm) 
#151329451..151332389
"Zm00001d042111"
myg <- subset(gid, chromosome_name==3 & start_position >= 151320000 & end_position <= 151333389)

#seq <- getSequence(chromosome=3, start = 151320000 -BP, end = 151333389 + BP, mart = zm)
#seq <- getSequence(id = "Zm00001d042111", mart = zm)



# bcftools view agpv4_chr3_agpv3_140-160M.vcf.gz -r 3:151300000-15140000 -Ov -o agpv4_chr3_151300000-151400000.vcf
# scp -P 2022 farm:/home/jolyang/dbcenter/HapMap/HapMap3/agpv4_chr3_agpv3_140-160M.vcf.gz ~/Desktop
########################################
load("largedata/dong_gwas_agpv4.RData")

plot(res1.mm, main="MainSpikeLength", pch=16, col="cadetblue")
plot(res2.mm, main="NumberofTilleringPlants", pch=16, col="cadetblue")
#plot(res2.mm, df="Pc1df", main="NumberofTilleringPlants", pch=16, col="cadetblue")
abline(v=151332389)

res2 <- results(res2.mm)
plot(res3.mm, main="10 kernel weight", pch=16, col="cadetblue")

plot(res4.mm, main="Spikelets.MainSpike", pch=16, col="cadetblue")
abline(v=151332389)

plot(res3.mm, main="Spikelets.MainSpike", pch=16, col="cadetblue")
plot(res4.mm, main="NumberofTilleringPlants", pch=16, col="cadetblue")
plot(res5.mm, main="Spikelets.PrimaryBranch", pch=16, col="cadetblue")
#plot(res6.mm, main="TasselBranchLength", pch=16, col="cadetblue")
#plot(res7.mm, main="10 kernel weight", pch=16, col="cadetblue")
plot(res8.mm, main="10 kernel weight", pch=16, col="cadetblue")
#abline(v=151332389)

BP <- 20000

res2 <- subset(results(res2.mm), !is.na(P1df))
pdf("graphs/fig3.pdf", width=8, height=4)
BP <- 50000
res2 <- subset(results(res2.mm), !is.na(P1df))
get_plot(res=res2, mart=zm, BP, pos=c(151328862 - BP, 151332856 + BP), 
         cutoff=10^-2.8, trait="Number of Tillering Plants", ylim=c(0, 5))
dev.off()

get_plot(res=res2, mart=zm, BP, pos=c(151328862 - BP, 151332856 + BP), 
         trait="Number of Tillering Plants", ylim=c(0, 5))

pdf("graphs/fig1.pdf", width=8, height=4)
BP <- 50000
res3 <- subset(results(res3.mm), !is.na(P1df))
get_plot(res=res3, mart=zm, BP, pos=c(151328862 - BP, 151332856 + BP), 
         trait="Secondary Branch Number", ylim=c(0, 7))
dev.off()

subset(res3, P1df == min(res3$P1df))


res4 <- subset(results(res4.mm), !is.na(P1df))
get_plot(res=res4, mart=zm, BP, pos=c(151328862 - BP, 151332856 + BP), 
         trait="Spikelets MainSpike", ylim=c(0, 3.5))

BP <- 50000
res5 <- subset(results(res5.mm), !is.na(P1df))
get_plot(res=res5, mart=zm, BP, pos=c(151328862 - BP, 151332856 + BP), 
         trait="Spikelets Primary Branch", ylim=c(0, 4))

res6 <- subset(results(res6.mm), !is.na(P1df))
get_plot(res=res6, mart=zm, BP, pos=c(151328862 - BP, 151332856 + BP), 
         trait="Tassel Branch Length", ylim=c(0, 3.5))

pdf("graphs/fig2.pdf", width=8, height=4)
BP <- 50000
res8 <- subset(results(res8.mm), !is.na(P1df))
get_plot(res=res8, mart=zm, BP, pos=c(151328862 - BP, 151332856 + BP), 
         trait="Tassel Primary Branches", ylim=c(0, 7))
dev.off()
subset(res8, P1df == min(res8$P1df))
# 3_150095715_C_T

bed <- read.table("largedata/output_peak_chr3.bed")
subset(bed, V2 > 150328862 & V3 < 153469526)

##############    
get_plot <- function(res, mart=zm, BP, pos=c(151328862 - BP, 151332856 + BP), cutoff=NULL, trait, ylim=c(0, 7)){
    
    res <- subset(res, Position > pos[1] & Position < pos[2])
    res$qval <- p.adjust(res$P1df, method = "fdr")
    res <- res[order(res$qval), ]
    out <- subset(res, qval < 0.05)
    
    if(nrow(out) == 0){
        message("### NO SNP detected! smallest qval =[ %s ]", min(res$qval))
    }else{
        res$l <- abs(res$qval - 0.05)
        mycutoff <- res[res$l == min(res$l),]$P1df
    }
    if(!is.null(cutoff)){
        mycutoff <- cutoff
    }
    
    seg <- makeSegmentation(start=pos[1], end= pos[2], value= -log10(mycutoff)[1],
                            dp = DisplayPars(color = "darkred", lwd=1, lty = 2))
    
    pList = list(title <- makeTitle(text = trait, color = "darkred"),
                 
                 "-log10(P value)" = makeBaseTrack(base = res$Position, value = -log10(res$P1df), 
                                                   trackOverlay=seg,
                                                   dp = DisplayPars(lwd = .6, color = "darkblue", ylim = ylim)),
                 
                 "Gene" = makeGeneRegion(chromosome = 3, start = pos[1], end = pos[2], 
                                         strand = "+", biomart = mart, 
                                         dp = DisplayPars(plotId = FALSE, idRotation = 0, cex = 1.2)),
                 
                 makeGenomeAxis(dp = DisplayPars(size = 1)))
    
    gdPlot(pList, minBase = pos[1], maxBase = pos[2])
}




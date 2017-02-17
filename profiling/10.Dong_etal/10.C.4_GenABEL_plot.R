

library("GenomeGraphs")
########################################
load("largedata/dong_gwas_agpv4.RData")

#plot(res1.mm, main="MainSpikeLength", pch=16, col="cadetblue")
plot(res2.mm, main="NumberofTilleringPlants", pch=16, col="cadetblue")
#plot(res2.mm, df="Pc1df", main="NumberofTilleringPlants", pch=16, col="cadetblue")
abline(v=151332389)

res2 <- results(res2.mm)
#plot(res3.mm, main="10 kernel weight", pch=16, col="cadetblue")

plot(res4.mm, main="Spikelets.MainSpike", pch=16, col="cadetblue")
abline(v=151332389)

#plot(res5.mm, main="Spikelets.PrimaryBranch", pch=16, col="cadetblue")
#plot(res6.mm, main="TasselBranchLength", pch=16, col="cadetblue")
#plot(res7.mm, main="10 kernel weight", pch=16, col="cadetblue")
#plot(res8.mm, main="10 kernel weight", pch=16, col="cadetblue")
#abline(v=151332389)


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

gene <- makeGene(id = "Zm00001d042111", type="ensembl_gene_id", biomart = zm)
transcript <- makeTranscript(id = "Zm00001d042111", type="ensembl_gene_id", biomart = zm)
gdPlot(list(gene, transcript))

seq <- getSequence(chromosome=3, start = 151320000 -BP, end = 151333389 + BP, mart = zm)
seq <- getSequence(id = "Zm00001d042111", mart = zm)


BP <- 100000
pList = list("-log10(P value)" = makeBaseTrack(base = res2$Position, value = -log10(res2$P1df), 
                                           dp = DisplayPars(lwd = .3, cex=2, color = "darkblue", ylim = c(0,3.5))),
             
             "Gene" = makeGeneRegion(chromosome = 3, start = 151320000 -BP, end = 151333389 + BP, 
                                  strand = "+", biomart = zm, 
                                  dp = DisplayPars(plotId = TRUE, idRotation = 0, cex = 1.2)),
             "Pos" = makeGeneRegion(chromosome = 3, start = 151320000 -BP, end = 151333389 + BP, 
                                         strand = "-", biomart = zm, 
                                         dp = DisplayPars(plotId = FALSE, idRotation = 90, cex = 1.2)),
             makeGenomeAxis(dp = DisplayPars(size = 1)))

gdPlot(pList, minBase = 151320000 -BP, maxBase = 151333389 + BP, 
       overlay = makeRectangleOverlay(start = 151320000 -BP, end = 151333389 + BP, region = c(4,8), 
                                      dp = DisplayPars(alpha = .5)))

# bcftools view agpv4_chr3_agpv3_140-160M.vcf.gz -r 3:151300000-15140000 -Ov -o agpv4_chr3_151300000-151400000.vcf
# scp -P 2022 farm:/home/jolyang/dbcenter/HapMap/HapMap3/agpv4_chr3_agpv3_140-160M.vcf.gz ~/Desktop


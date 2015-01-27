# Jinliang Yang
# 1.26.2015
# purpose: compute the gene expression complementation for the bybrids

### reformat the data
gene <- read.delim("largedata/1.gc/maize_503genotypes_raw_expression_values_FPKM.txt", header=TRUE)
names(gene) <- c(names(gene)[-1], "N")
gene <- gene[, 1:506]


idx <- grepl("joint", row.names(gene))


gene[is.na(gene)] <- 0.00001
gene[gene==0] <- 0.0001
hist(-log10(t(gene[6, 4:506])), main="GRMZM2G125201", xlab="-log10(FPKM)")

range(gene$Mo17)




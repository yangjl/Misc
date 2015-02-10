### Jinliang Yang
### Feb 9th, 2015
### check whether the coordinates are AGPv2


######################## Filtering the duplicated annotation #################
fgsv2 <- read.table("~/dbcenter/AGP/AGPv2/ZmB73_5b_FGS.gff", header=FALSE)
names(fgsv2) <- c("seqname", "source", "feature", "start", "end", "score",
                  "strand", "frame", "attribute")

gene <- subset(fgsv2, feature=="gene")
#39656
gene$attribute <- gsub("ID=", "", gene$attribute)
gene$attribute <- gsub(";.*", "", gene$attribute)
gene <- subset(gene, seqname %in% 1:10)

bed4gene <- gene[, c("seqname", "start", "end", "attribute")]
bed4gene$start <- bed4gene$start - 1

write.table(bed4gene, "largedata/1.gc/FGSv2.bed4", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")


######################## Pan-transcriptome data #################
exp <- read.csv("largedata/1.gc/maize_gene_503lines.csv")
expos <- exp[, 1:3]

######################## merge the two #################

two <- merge(gene, expos, by.x="attribute", by.y="row.names")
#[1] 39455    12

# note: Pan-transcriptome data is 0-based coding

sum(two$start-1 != two$position_left)
which(two$end != two$position_right)




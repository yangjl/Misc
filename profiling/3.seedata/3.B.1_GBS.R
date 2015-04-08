### Jinliang Yang
### 4/7/2015
### transform GBS format to BED+ format with haplotype call

seeds <- read.delim("data/seeds_09.02.2015_22.38.10.txt")
#[1] 22022    51
subseed <- subset(seeds, general_identifier %in% idtab$GID)




idtab <- read.csv("data/SeeD_SID_to_GID.csv")
length(unique(idtab$GID))
length(unique(idtab$SampleID))
# Note: some of the accessions were genotyped multiple times


#####
gbs10 <- read.delim("~/dbcenter/seeds_data/All_SeeD_2.7_chr10_no_filter.unimputed.hmp.txt", header=TRUE)

gbs2bed(gbsfile="~/dbcenter/seeds_data/All_SeeD_2.7_chr10_no_filter.unimputed.hmp.txt",
        outfile="~/dbcenter/seeds_data/chr10_filetered_unimputed.hmp")
    

gbs2bed <- function(gbsfile="~/dbcenter/seeds_data/All_SeeD_2.7_chr10_no_filter.unimputed.hmp.txt",
                    outfile="~/dbcenter/seeds_data/chr10_filetered_unimputed.hmp"){
    
    ### read in GBS file
    tab5 <- read.delim(gbsfile, nrows=5, header=TRUE)
    classes <- sapply(tab5, class)
    gbs <- read.delim(gbsfile, header=TRUE, colClasses=classes)
    message(sprintf("Loaded [ %s ] SNPs and [ %s ] cols for file [%s]!", nrow(gbs), ncol(gbs), gbsfile))
    
    ### change to BED5+ format
    gbs <- gbs[, c(3,4,4,1,2,5, 12:ncol(gbs))]
    names(gbs)[1:6] <- c("chr", "start", "end", "snpid", "alleles", "nchar")
    nms <- names(gbs)
    nms <- gsub("\\..*$", "", nms)
    names(gbs) <- nms
    gbs$start <- gbs$start -1
    message(sprintf("Changed to BED5+ format and start filtering ..."))
    
    ### filter SNPs contain multiple alleles
    gbs$nchar <- nchar(as.character(gbs$alleles))
    subg <- subset(gbs, nchar == 3)
    subg <- subg[, -6]
    #idx <- grep("-", subg$alleles)
    #subg <- subg[-idx,]
    message(sprintf("Remaining [ %s ] sites with two variations!", nrow(subg)))
    
    message(sprintf("Start to transforming, recoding and writing ..."))
    ### change to two identical haplotypes
    for(i in 6:ncol(subg)){
        subg[, i] <- paste(subg[, i], subg[, i])
        #print(i)
    }
    
    ###change IUPAC Ambiguity Codes
    #M    A or C	K
    #R	A or G	Y
    #W	A or T	W
    #S	C or G	S
    #Y	C or T	R
    #K	G or T	M
    subg[subg=="M M"] <- "A C"
    subg[subg=="R R"] <- "A G"
    subg[subg=="W W"] <- "A T"
    subg[subg=="S S"] <- "C G"
    subg[subg=="Y Y"] <- "C T"
    subg[subg=="K K"] <- "G T"
    subg[subg=="0 0"] <- "N N"
    
    write.table(subg, outfile, sep="\t", row.names=TRUE, col.names=FALSE, quote=FALSE)
    message(sprintf("DONE!"))
}





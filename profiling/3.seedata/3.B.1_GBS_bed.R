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


##### transform GBS to BED format
source("lib/gbs2bed.R")
for(i in 5:9){
    gbs2bed(gbsfile= paste0("~/dbcenter/seeds_data/All_SeeD_2.7_chr", i, "_no_filter.unimputed.hmp.txt"),
            outfile= paste0("~/dbcenter/seeds_data/chr", i, "_filetered_unimputed.hmp"))
    
}

##### run the following python code to get the SNP frq and missing rate
snpfrq -p ~/dbcenter/seeds_data/ -i chr10_filetered_unimputed.hmp -s 6 -m "0N" -a 0 -b 1 -c 2 -o chr10_filetered_unimputed.frq
snpfrq -p ~/dbcenter/seeds_data/ -i chr9_filetered_unimputed.hmp -s 6 -m "0N" -a 0 -b 1 -c 2 -o chr9_filetered_unimputed.frq
snpfrq -p ~/dbcenter/seeds_data/ -i chr8_filetered_unimputed.hmp -s 6 -m "0N" -a 0 -b 1 -c 2 -o chr8_filetered_unimputed.frq
snpfrq -p ~/dbcenter/seeds_data/ -i chr7_filetered_unimputed.hmp -s 6 -m "0N" -a 0 -b 1 -c 2 -o chr7_filetered_unimputed.frq
snpfrq -p ~/dbcenter/seeds_data/ -i chr6_filetered_unimputed.hmp -s 6 -m "0N" -a 0 -b 1 -c 2 -o chr6_filetered_unimputed.frq
snpfrq -p ~/dbcenter/seeds_data/ -i chr5_filetered_unimputed.hmp -s 6 -m "0N" -a 0 -b 1 -c 2 -o chr5_filetered_unimputed.frq
snpfrq -p ~/dbcenter/seeds_data/ -i chr4_filetered_unimputed.hmp -s 6 -m "0N" -a 0 -b 1 -c 2 -o chr4_filetered_unimputed.frq
snpfrq -p ~/dbcenter/seeds_data/ -i chr3_filetered_unimputed.hmp -s 6 -m "0N" -a 0 -b 1 -c 2 -o chr3_filetered_unimputed.frq
snpfrq -p ~/dbcenter/seeds_data/ -i chr2_filetered_unimputed.hmp -s 6 -m "0N" -a 0 -b 1 -c 2 -o chr2_filetered_unimputed.frq
snpfrq -p ~/dbcenter/seeds_data/ -i chr1_filetered_unimputed.hmp -s 6 -m "0N" -a 0 -b 1 -c 2 -o chr1_filetered_unimputed.frq








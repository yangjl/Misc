### Jinliang Yang
### 12/4/2015


library("tasselr")
library("imputeR")

# Note: at least 100G memory was needed to load the hdf5 file
# load h5file
allzea <- initTasselHDF5("~/dbcenter/AllZeaGBS/ZeaGBSv27_publicSamples_imputedV5_AGPv2-150114.h5", version="5")
allzea <- loadBiallelicGenotypes(allzea, verbose = TRUE)


#ncoding genotypes... Error in loadBiallelicGenotypes(allzea, verbose = TRUE) : 
#    genotype includes allele that's not ref or alt


allzea@ranges         allzea@ref            allzea@alt            
  allzea@genotypes      allzea@samples

geno <- allzea@genotypes

##### copy the following codes to shell
srun.x11 -p bigmemh --ntasks=8 --nodelist=bigmem4
module load gcc jdk/1.8 tassel/5

# to export - whatever the hdf5 file is - if you want plink -'Plink' or 'Hapmap' for hmp
run_pipeline.pl -Xmx64g -fork1 -h5 ZeaGBSv27_publicSamples_imputedV5_AGPv2-150114.h5 \\
-includeTaxaInfile /home/jolyang/Documents/Github/SeeDs/data/Taxa_ames_3324.txt -export ZeaGBSv27_Ames_3324.h5 -exportType HDF5 -runfork1

#first you'll need to sort the GBS 2.7 - this takes a while (1 hour) - so just copy my version from here:/group/jrigrp4/Justin_Kate/GBS2.7
# 'sortedGBS.hmp.txt' - if you want to do it yourself
run_pipeline.pl -Xmx64g -SortGenotypeFilePlugin -inputFile ZeaGBSv27_publicSamples_imputedV5_AGPv2-150114.hmp.txt -outputFile ZeaGBSv27_sorted -fileType Hapmap

# for GBS 2.7, the "keep_list_NAM_children.txt" is just a one column list of GBS runs I want to keep
run_pipeline.pl -Xmx64g -fork1 -h ZeaGBSv27_sorted.hmp.txt -includeTaxaInfile /home/jolyang/Documents/Github/SeeDs/data/Taxa_ames_3324.txt -export ZeaGBSv27_Ames -exportType Hapmap -runfork1







# reformat to imputeR object
ob1 <- imputeRob(h5=teo, missingcode=3)
save(file="largedata/teo.RData", list="ob1")

#calculating missing rates ... done.


allzea <- read.csv("data//AllZeaGBSv2.7_publicSamples_metadata20140411.csv")
table(allzea$Project)
table(allzea$GermplasmSet)
ames <- subset(allzea, GermplasmSet %in% "Ames" & Project != "Imputation Test")
ames$taxa <- paste(ames$DNASample, ames$LibraryPrepID, sep=":")
#2010 Ames Lines        AMES Inbreds             Ames282      ApeKI 384-plex 
#2628                 696
write.table(ames$taxa, "data/Taxa_ames_3324.txt", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
# note this ames lines including 2010 Ames lines and 696 Ames inbred!
# Note: space delimted => %s/\n/ /g in vim


##### copy the following codes to shell
srun.x11 -p bigmemh --ntasks=8 --nodelist=bigmem4
module load gcc jdk/1.8 tassel/5

# to export - whatever the hdf5 file is - if you want plink -'Plink' or 'Hapmap' for hmp
run_pipeline.pl -Xmx64g -fork1 -h5 ZeaGBSv27_publicSamples_imputedV5_AGPv2-150114.h5 -export -exportType Hapmap -runfork1

#first you'll need to sort the GBS 2.7 - this takes a while (1 hour) - so just copy my version from here:/group/jrigrp4/Justin_Kate/GBS2.7
# 'sortedGBS.hmp.txt' - if you want to do it yourself
run_pipeline.pl -Xmx64g -SortGenotypeFilePlugin -inputFile ZeaGBSv27_publicSamples_imputedV5_AGPv2-150114.hmp.txt -outputFile ZeaGBSv27_sorted -fileType Hapmap

# for GBS 2.7, the "keep_list_NAM_children.txt" is just a one column list of GBS runs I want to keep
run_pipeline.pl -Xmx64g -fork1 -h ZeaGBSv27_sorted.hmp.txt -includeTaxaInfile /home/jolyang/Documents/Github/SeeDs/data/Taxa_ames_3324.txt -export ZeaGBSv27_Ames -exportType Hapmap -runfork1


# Note: at least 100G memory was needed to load the hdf5 file
# load h5file
test <- initTasselHDF5("~/dbcenter/AllZeaGBS/ZeaGBSv27_publicSamples_imputedV5_AGPv2-150114.h5", version="5")
gbs <- loadBiallelicGenotypes(test, verbose = TRUE)

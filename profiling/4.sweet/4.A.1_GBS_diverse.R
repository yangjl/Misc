### Jinliang Yang
### May 8th, 2015
### find association for the sweet gene


allzea <- read.csv("~/Documents/Github/SeeDs/data/AllZeaGBSv2.7_publicSamples_metadata20140411.csv")
table(allzea$Project)
table(allzea$GermplasmSet)
ames282 <- subset(allzea, Project %in% "Ames282") 
length(unique(ames282$DNASample)) #283                

ames282$taxa <- paste(ames282$DNASample, ames282$LibraryPrepID, sep=":")
#2010 Ames Lines        AMES Inbreds             Ames282      ApeKI 384-plex 
#2628                 696
write.table(ames282$taxa, "data/Taxa_ames282_288.txt", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
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




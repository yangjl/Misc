### Jinliang Yang

allzea <- read.csv("data/AllZeaGBSv2.7_publicSamples_metadata20140411.csv")
dim(allzea)
#[1] 17280    18

table(allzea$Project)
table(allzea$GermplasmSet)

nam <- subset(allzea, GermplasmSet %in% "NAM")

nam$taxa <- paste(nam$DNASample, nam$LibraryPrepID, sep=":")

write.table(nam$taxa, "~/dbcenter/AllZeaGBS/Taxa_nam_5873.txt", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
# Note: space delimted => %s/\n/ /g in vim




##### copy the following codes to shell
srun.x11 -p bigmemh --ntasks=8 --nodelist=bigmem4
module load gcc jdk/1.8 tassel/5

# to export - whatever the hdf5 file is - if you want plink -'Plink' or 'Hapmap' for hmp
run_pipeline.pl -Xmx64g -fork1 -h5 ZeaGBSv27_publicSamples_imputedV5_AGPv2-150114.h5 \\
-includeTaxaInfile /home/jolyang/Documents/Github/SeeDs/data/Taxa_ames_3324.txt -export ZeaGBSv27_Ames_3324.h5 -exportType HDF5 -runfork1


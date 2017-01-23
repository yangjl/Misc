### Jinliang Yang
### purpose: transform ZeaGBSv2.7 Ames inbred => bed+ format

######=======>
library(data.table, lib = "~/bin/Rlib")

source("lib/gbs2bed_ames282.R")
gbs2bed_ames(gbsfile="~/dbcenter/AllZeaGBS/ZeaGBSv27_Ames282_agpv3.hmp.txt",
             outfile="~/dbcenter/AllZeaGBS/ZeaGBSv27_Ames282_agpv3.bed5")

### log
# Read 954794 rows and 299 (of 299) columns from 0.554 GB file in 00:00:35
# Changed to BED5+ format and start filtering ...
# Remaining [ 515770 ] sites with two variations!
#    Start to IUPAC=>N transforming, recoding and writing ...

# not working
# snpfrq -i ZeaGBSv27_Ames282_agpv3.bed5 -s 6 -m "N" -a 0 -b 1 -o ZeaGBSv27_Ames282_agpv3.frq
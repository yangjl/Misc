
# load modules
module load gcc jdk/1.8 tassel/5
# to export hdf5 file to Hapmap format
run_pipeline.pl -Xmx64g -fork1 -h5 ZeaGBSv27_publicSamples_imputedV5_AGPv2-150114.h5 \\
-export -exportType Hapmap -runfork1

# sort the GBS2.7 - this takes a while (1 hour)
run_pipeline.pl -Xmx64g -SortGenotypeFilePlugin -inputFile AllZeaGBSv2.7_publicSamples_imputedV3b_agpv3.hmp.gz -outputFile ZeaGBSv27_agpv3_sorted -fileType Hapmap

# for GBS2.7, the "keep_list_NAM_children.txt" is just a one column list to keep
run_pipeline.pl -Xmx64g -fork1 -h ZeaGBSv27_agpv3_sorted.hmp.txt -includeTaxaInfile ~/Documents/Github/Misc/data/FullName_ames282_288.txt -export ZeaGBSv27_Ames282_agpv3 -exportType Hapmap -runfork1



# load modules
module load gcc jdk/1.8 tassel/5
# to export hdf5 file to Hapmap format
run_pipeline.pl -Xmx64g -fork1 -vcf agpv4_chr3_agpv3_140-160M.vcf.gz -export -exportType Hapmap -runfork1

# sort the GBS2.7 - this takes a while (1 hour)
run_pipeline.pl -Xmx64g -SortGenotypeFilePlugin -inputFile AllZeaGBSv2.7_publicSamples_imputedV3b_agpv3.hmp.gz -outputFile ZeaGBSv27_agpv3_sorted -fileType Hapmap

# for GBS2.7, the "keep_list_NAM_children.txt" is just a one column list to keep
run_pipeline.pl -Xmx64g -fork1 -h ZeaGBSv27_agpv3_sorted.hmp.txt -includeTaxaInfile ~/Documents/Github/Misc/data/FullName_ames282_288.txt -export ZeaGBSv27_Ames282_agpv3 -exportType Hapmap -runfork1



year <- as.character(seeds$colldate)

year <- gsub("-.*", "", year)
year <- as.numeric(as.character(year))
hist(year)


######
install.packages("raster")
library(raster)

#Get bio variables at 10 minutes of a degree resolution
biovar = getData('worldclim', var='bio', res=10)

biovar5 = getData('worldclim', var="bio", res=5)

w = getData('worldclim', var='tmin', res=0.5, lon=5, lat=45)
###Get elevation data with SRTM - go to website to determine relevant tiles, 
#each of these tiles are large - may be easier to do on website

#elev.europe <- getData('SRTM', lon = 5, lat = 45)
#elev.europe1 <- getData('SRTM', lon = 0, lat = 40)
#elev.europe2 <- getData('SRTM', lon = -5, lat = 35)



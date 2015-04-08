### Jinliang Yang
### experiments on map plotting

install.packages("dismo")

library(dismo)
library(maps)       # Provides functions that let us plot the maps
library(mapdata)    # Contains the hi-resolution points that mark out the countries.
library(ggmap)
map('worldHires')

seeds <- read.delim("data/seeds_09.02.2015_22.38.10.txt")

for(i in 1:nrow(seeds)){
    points(seeds$longitude[i], seeds$latitude[i],  col=2, pch=18)
}

points(-1.615672,54.977768,col=2,pch=18)



### 
library(rworldmap)

par(mai=c(0,0,0.2,0), xaxs="i", yaxs="i")
mapBubbles(dF=getMap(), oceanCol='lightblue', landCol='wheat')



library(rworldmap)
newmap <- getMap(resolution = "coarse")
plot(newmap, oceanCol='lightblue', landCol='wheat')

plot(newmap, xlim = c(-20, 59), ylim = c(35, 71), asp = 1)


> library(ggmap)
> map <- get_map(location = 'Europe', zoom = 4)
> mapPoints <- ggmap(map) +
    +   geom_point(aes(x = lon, y = lat, size = sqrt(flights)), data = airportD, alpha = .5)

library(ggmap)
mapImageData1 <- get_map(location = c(lon = -0.016179, lat = 51.538525),
                         color = "color",
                         source = "google",
                         maptype = "satellite",
                         zoom = 17)
ggmap(mapImageData1,
      extent = "device",
      ylab = "Latitude",
      xlab = "Longitude")

geom_point(aes(x = Longitude, y = Latitude), data = data,
           alpha = .5, color="darkred", size = 3)

##lowerleftlon, lowerleftlat, upperrightlon, upperrightlat
myloc <- c(-130, -56, -34, 55)
mymap <- get_map(location=myloc, source="google", crop=TRUE, color="bw")
ggmap(mymap) + 
geom_point(aes(x = longitude, y = latitude), data = seeds,
           alpha = .5, color="darkred", size = 1)

### location of Totontepec
#17.216667, -95.983333


#http://gis.stackexchange.com/questions/119736/ggmap-create-circle-symbol-where-radius-represents-distance-miles-or-km
d <- data.frame(lat = c(33.79245), lon = c(-84.32130))

coordinates(d) <- ~ lon + lat
projection(d) <- "+init=epsg:4326"

emory <- gmap(extent(-130, -56, -34, 55), zoom = 10, scale = 2)
d_mrc <- spTransform(d, CRS = CRS(projection(mymap)))
# Miles to meters conversion
mile2meter <- function(x) {
    x * 1609.344
}

# Buffer creation
d_mrc_bff <- gBuffer(d_mrc, width = mile2meter(20))
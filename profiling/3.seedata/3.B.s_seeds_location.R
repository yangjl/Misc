### Jinliang Yang
### 4.12.2015
### find SeeDs distance relative to Totontepec

getDistance <- function(df=idinfo, lon=-95.983333,lan=17.216667){
    
    df$dist <- sqrt((df$longitude - lon)^2 + (df$latitude - lan)^2)
    return(df)
}
#########
DefineDist <- function(){
    locinfo <- read.delim("data/seeds_09.02.2015_22.38.10.txt", header=TRUE)
    
    ids <- read.csv("data/SeeD_SID_to_GID.csv")
    # 4710    2
    # Note: some of about 600 with no collection information
    length(unique(ids$GID))
    #4020
    # Note: 
    idinfo <- merge(ids, locinfo, by.x="GID", by.y="general_identifier")
    length(unique(idinfo$GID))
    #3493
    
    ### location of Totontepec
    ### -95.983333
    ### 17.216667
    
    idinfo <- getDistance(df=idinfo, lon=-95.983333,lan=17.216667)
    return(idinfo)    
}

idinfo <- DefineDist()
write.table(idinfo, "data/SeeD_idinfo.csv", sep=",", row.names=FALSE, quote=FALSE)



hist(idinfo$elevation, breaks=50, main="SeeDs Data", xlab="Elevation", col="grey")
######################################################################################
library(ggmap)
##lowerleftlon, lowerleftlat, upperrightlon, upperrightlat
myloc <- c(-130, -56, -34, 55)

x=-95.983333
y=17.216667
r=1
angles <- seq(0, 2*pi, length.out=360)
mymap <- get_map(location=myloc, source="google", crop=TRUE, color="bw")
ggmap(mymap) + 
    geom_point(aes(x = longitude, y = latitude), data = idinfo,
               alpha = .5, color="darkred", size = 1) +
    geom_line(aes(r*cos(angles)+x,r*sin(angles)+y), alpha=.6, color="blue")

### location of Totontepec
####################################################################
##lowerleftlon, lowerleftlat, upperrightlon, upperrightlat
myloc <- c(-100, 10, -90, 25)

x=-95.983333
y=17.216667
r=1
angles <- seq(0, 2*pi, length.out=360)
mymap <- get_map(location=myloc, source="google", crop=TRUE, color="bw")
ggmap(mymap) + 
    geom_point(aes(x = longitude, y = latitude), data = idinfo,
               alpha = .5, color="darkred", size = 2) +
    geom_line(aes(r*cos(angles)+x,r*sin(angles)+y), alpha=.3, color="blue")


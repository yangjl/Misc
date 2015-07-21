### Jinliang Yang
### July 21th, 2015

#"Variable   Description                   Units"
#"   1       Database name                 "
#"   2       Date: year,month,day          yyyymmdd"
#"   3       Observation time              hhmm"
#"   4       Precipitation, amount         Millimeters"
#"   5       Precipitation, type           (coded)"
#"   6       Air temperature, maximum      Celsius"
#"   7       Air temperature, minimum      Celsius"
#"   8       Air temperature, observed     Celsius"
#"   9       Weather conditions            (coded)"
#"  10       Wind, direction               N,NE,E,SE,S,SW,W,NW, 0=calm"
#"  11       Wind, speed                   Meters per second"
#"  12       Bulb temperature, wet         Celsius"
#"  13       Bulb temperature, dry         Celsius"
#"  14       Soil temperature, maximum     Celsius"
#"  15       Soil temperature, minimum     Celsius"
#"  16       Pan evaporation               Millimeters"
#"  17       Solar radiation               Watts per sq. meters"
#"  18       Reference evapotranspiration  Millimeters"
#"  19       Relative humidity, minimum    Percent"
#"  20       Relative humidity, maximum    Percent"


davis <- read.csv("data/davis_weather.csv")

davis$year <- gsub("....$", "", davis$Date)
davis <- subset(davis, !is.na(Date))
davis$md <- gsub("^....", "", davis$Date)
davis$month <- gsub("..$", "", davis$md)
davis$day <- gsub("^..", "", davis$md)

davis$month <- as.numeric(as.character(davis$month))
davis$day <- as.numeric(as.character(davis$day))


######################
plot(0, 0, xlim=c(1, 400), ylim=range(davis$Air.max, na.rm=TRUE), type="n", xlab="", ylab="Temperature (oC)")
for(i in unique(davis$year)){
    d <- subset(davis, year == i)
    for(j in 1:12){
        dm <- subset(d, month == j)
        abline(v=31*(dm$month-1), col="grey")
        lines(dm$day + 31*(dm$month-1), dm$Air.max, col=sample(colours(), 5))
    }  
}
    





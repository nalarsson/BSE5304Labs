# Natalie Larsson
# BSE 5304G
# Lab 01
# 1/25/2023

print("Hello World")

## SETUP -----
if (!require("pacman")) install.packages("pacman")
p_load(rgdal,parallel,ggplot2,dplyr,patchwork,rnoaa,hrbrthemes)
system("git config --global user.email 'nlarsson@vt.edu' ") 
system("git config --global user.name 'Natalie Larsson' ")


## PROBLEM 1 CODE -----

# hometown: Frederick, Maryland
ht_lat <- 39.45
ht_lon <- -77.50

# find stations near Frederick, Maryland
stns=meteo_distance(
  station_data=ghcnd_stations(),
  lat=ht_lat,
  long=ht_lon,
  units = "deg",
  radius = 20,
  limit = NULL
)

# using Sharpsburg, MD station
# station data only goes through March of 2022 -- OKed by Fuka
WXData=meteo_pull_monitors(
  monitors=stns[277,1],   
  keep_flags = FALSE,
  date_min = "2016-01-01",
  date_max = NULL,
  var = c("TMAX","TMIN","PRCP") 
)

#build data frame to preserve original data
weather <- data.frame(
  day <- WXData$date,
  tmin <- WXData$tmin*.1, #covert to celcius from tenth of degree celcius
  tmax <- WXData$tmax*.1,
  pr <- WXData$prcp*.01 #convert to cm from tenth of mm
)

# colors
pr.color <- "#167d7e"
tmin.color <- "#f57e00"
tmax.color <-"#FF0000"
  
# create ggplot
weather.chart <- ggplot(weather, aes(x=day)) +
  geom_line(aes(y=tmin), color=tmin.color, size=1.25) +
  geom_line(aes(y=tmax), color=tmax.color, size=1.25) +
  geom_line(aes(y=pr), color=pr.color, size=1.25) +
  ggtitle("Daily Precipitation, Minimum Temperature, \nand Maximum Temperature near Frederick, MD") +
  xlab("Year") + 
  scale_y_continuous(name="Temperature (°C)",
                     sec.axis=sec_axis(~., name="Precipitation (cm)")) +
  theme_ipsum() +
  theme(
    axis.title.y = element_text(color = tmax.color, size=13),
    axis.title.y.right = element_text(color = pr.color, size=13)
  )

weather.chart
pdf()
dev.off()


## PROBLEM 2 CODE -----
# call function from source
source("https://goo.gl/Cb8zGn")

myflowgage_id="01643000"
myflowgage=get_usgs_gage(myflowgage_id,
            begin_date="2016-02-01",end_date="2023-02-01")
class(myflowgage)

plot(myflowgage$flowdata$mdate,myflowgage$flowdata$flow,
     main=myflowgage$gagename,xlab = "Date",
     ylab="Flow m^3/day",type="l")

conversion <- 24*60*60
myflowgage$flowdata$flowcms <- myflowgage$flowdata$flow/conversion
myflowgage$flowdata$flowdepth <- myflowgage$flowdata$flowcms/myflowgage$area*1000

weatherandflow.chart <- weather.chart +
  geom_line(myflowgage$flowdata, mapping=aes(x=mdate, y=flowdepth)) +
  ggtitle("Daily Precipitation, Minimum Temperature, \n Maximum Temperature and Streamflow near Frederick, MD") +
  scale_y_continuous(name="Temperature (°C)",
                     sec.axis=sec_axis(~., name="Precipitation and Streamflow (cm)"))
                                                                                      
weatherandflow.chart

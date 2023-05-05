if(!require("pacman")) install.packages("pacman")
pacman::p_load(meteoForecast, EcoHydRology)
pacman::p_load(ggplot2,dplyr,patchwork,rnoaa)
pacman::p_load(operators,topmodel,DEoptim,soilDB,sp,curl,httr,
               rnoaa,raster,shapefiles,rgdal,elevatr,terra,progress,lubridate,zoo)


mygitdir=rstudioapi::getActiveProject()
mypdfdir=paste0(mygitdir,"/pdfs",LabNo)
dir.create(mypdfdir)
#
setwd(mygitdir)
system("git config --global user.email 'nlarsson@vt.edu' ")
system("git config --global user.name 'Natalie Larsson' ")
system("git config pull.rebase false")


Tempvars = grepVar("temp", service="gfs", complete=TRUE)
Precipvars = grepVar("precip", service="gfs", complete=TRUE)

testDay <- Sys.Date() - 1

install.packages("lattice")
library(lattice)
## Beware: the x-axis labels display time using your local timezone.

today=Sys.Date()

myflowgage_id="01654500" #Accotink Creek near Annandale VA # "0205551460"
myflowgage=get_usgs_gage(myflowgage_id, begin_date="2015-01-01",
                         end_date=today)

myflowgage$flowdata$Qmm = myflowgage$flowdata$flow/myflowgage$area/10^3

source("https://raw.githubusercontent.com/Rojakaveh/FillMissWX/main/FillMissWX.R")
WXData=FillMissWX(declat=myflowgage$declat, declon=myflowgage$declon,
                  StnRadius=30,minstns=10,date_min="2010-01-01",
                  date_max=today,targElev=myflowgage$elev,
                  method = "IDW",alfa=2)



## Multiple variables
Precipvars <- getPoint(c(myflowgage$declon, myflowgage$declat), 
                 var = c("Total_precipitation_surface_Mixed_intervals_Accumulation"), 
                 service = "gfs",
                 day = today)

GFSprecip = aggregate(Precipvars, as.Date(time(Precipvars)), sum)

Tempvars <- getPoint(c(myflowgage$declon, myflowgage$declat), 
                       var = c("Temperature_surface"), 
                       service = "gfs",
                       day = today)

GFSMaxTemp = aggregate(Tempvars, as.Date(time(Tempvars)), max)
GFSMinTemp = aggregate(Tempvars, as.Date(time(Tempvars)), min)

plot(GFSMinTemp, ylim=c(276,300))
lines(GFSMaxTemp)
title(main="Maximum and Minimum Temperature at Accotink Creek, VA", 
      xlab="Date", ylab="Temperature (K)")

# covert zoo class to data.frame, stack onto max temp, min temp precip
# add forecast data to WXData

BasinData=merge(WXData,myflowgage$flowdata,by.x="date",by.y="mdate")

#dir.create("~/pngs")
setwd("~/pngs")
graphdir="~/pngs"

# Temp graph
png("TempForecast.png")
plot(GFSMinTemp, ylim=c(276,310), 
     xlab="Date", ylab="Temperature (K)")
lines(GFSMaxTemp)
title(main="Forecasted Maximum and Minimum Temperature at Accotink Creek, VA")
dev.off()

# Precip graph
png("PrecipForecast.png")
plot(GFSprecip, 
     xlab="Date", ylab="Precipitation (kg/m^2)")
title(main="Forecasted Precipitation at Accotink Creek, VA")
dev.off()


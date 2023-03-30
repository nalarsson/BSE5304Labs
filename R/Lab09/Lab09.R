#Original Gage -----
LabNo="/Lab09"
myflowgage_id="0205551460"  # Old Friendly Gage
#
# What needs to be loaded
#
if (!require("pacman")) install.packages("pacman")
myhomedir=Sys.getenv("HOME")
datadir=paste0(myhomedir,"/data",LabNo)
dir.create(datadir,recursive = T)
srcdir=paste0(myhomedir,"/src")
dir.create(srcdir,recursive = T)

mygitdir=rstudioapi::getActiveProject()
mypdfdir=paste0(mygitdir,"/pdfs",LabNo)
dir.create(mypdfdir)
#
setwd(mygitdir)
system("git config --global user.email 'nlarsson@vt.edu' ")
system("git config --global user.name 'Natalie Larsson' ")
system("git config pull.rebase false")

pacman::p_load(httr,EcoHydRology,curl,elevatr,raster,rgdal,
                 data.table,foreign,maptools,dataRetrieval,gdistance)
setwd(datadir)
# Or to see what might be attached
siteNo = "0205551460"
search()
myflowgage=get_usgs_gage(siteNo, begin_date ="2022-01-01", end_date="2023-01-01")
intersect(search(), objects())
objects()  # This will list the objects you have.
rm(list=objects()) # Removes ALL the objects… so be careful here.

# https://owi.usgs.gov/R/dataRetrieval.html
# https://owi.usgs.gov/R/training-curriculum/usgs-packages/dataRetrieval-readNWIS/

# ?dataRetrieval  # Review the man page for this package
# ?readNWISuv
# ??readNWISdv
# ?readNWISdata
#
# Need to figure out which data to download. 
# https://nwis.waterdata.usgs.gov/nwis/pmcodes?radio_pm_search=param_group&pm_group=All+--+include+all+parameter+groups&pm_search=&casrn_search=&srsname_search=&format=html_table&show=parameter_group_nm&show=parameter_nm&show=casrn&show=srsname&show=parameter_units
# 
url="https://nwis.waterdata.usgs.gov/nwis/pmcodes?radio_pm_search=param_group&pm_group=All+--+include+all+parameter+groups&pm_search=&casrn_search=&srsname_search=&format=html_table&show=parameter_group_nm&show=parameter_nm&show=casrn&show=srsname&show=parameter_units"
#browseURL(url)
#
# Yeah, these databases are complex to get to know, remember our 
# SSURGO?
#
# Before you begin your modeling project, confirm what your model outputs 
# has a value to calibrate against, i.e. match parameter and units. For 
# this lab we are looking for Gage Height, while historically, we have been 
# looking at Discharge. NOT ALL PARAMETERS ARE AVAILABLE!
#
url="https://help.waterdata.usgs.gov/parameter_cd?group_cd=%"
#browseURL(url)

#View(parameterCdFile)
##############################################
# 0205551460 LICK RUN ABOVE PATTON AVENUE AT ROANOKE, VA
##############################################
make_usgs_gage_list = function(siteNo = "0205551460", 
                               parameterCd = c("00060","00065"),
                               start.date = "2017-05-01",  # Not frozen to not frozen
                               end.date = "2017-11-01"  )  {# to still not frozen
  #
  # For each gage location, let's keep the data organized as a 
  # list.
  USGS=list()   # Organize the data in a nice list as in previous labs
  USGS[["flowdata"]]<- get_usgs_gage(siteNo, start.date,end.date)
  USGS[["flowdata"]]<- readNWISuv(siteNumbers = siteNo,parameterCd = parameterCd,startDate = start.date,endDate = end.date)
  #head(USGS$flowdata)  # Note that we have 00060 and 00065...
  #  agency_cd	site_no        	dateTime X_00060_00000 X_00060_00000_cd
  #1  	USGS 0205551460 2017-05-01 04:00:00      	6.38            	A
  #2  	USGS 0205551460 2017-05-01 04:05:00      	6.38            	A
  #  X_00065_00000 X_00065_00000_cd tz_cd
  #1      	2.74            	A   UTC
  #2      	2.74            	A   UTC
  #
  # And of course we want to work in SI units so:
  USGS$flowdata$depth_m=USGS$flowdata$X_00065_00000*0.3048
  # m/ft depth
  USGS$flowdata$cms=USGS$flowdata$X_00060_00000*.02832
  # m3/ft3 flow
  #
  # Let's add in the USGS gage site information to the list and inspect
  USGS[["site"]]=readNWISsite(siteNo)
  #head(USGS$site)
  class(USGS$site$dec_lat_va)
  #
  # Set the Manning Coefficient in the USGS Gage's Site Table
  #
  url="https://www.google.com/search?q=manning%27s+n+for+stream"
  #browseURL(url)
  url="https://www.fsl.orst.edu/geowater/FX3/help/8_Hydraulic_Reference/Mannings_n_Tables.htm"
  #browseURL(url)
  USGS$site$man_n=.035/1.49
  #
  # Create a SpatialPointsDataFrame out of the site dataframe in the USGS list
  coordinates(USGS$site)=~dec_long_va+dec_lat_va
  return(USGS)
}

# pacman::p_load(useful)
# compare.list(testlist,USGS0205551460)
#Might be nice to have a function that gets these data, no?

# Of course since somebodies said they hate making functions we
# will go through this here to show how fun it is!
#

# If done correctly, your function should be able to populate
# the lists for the remaining gages!
#

#gages from homework!!
USGS02056000=make_usgs_gage_list(siteNo = "02056000")
USGS0205551460=make_usgs_gage_list(siteNo ="0205551460" )
USGS02055100=make_usgs_gage_list(siteNo ="02055100" )
USGS02055000=make_usgs_gage_list(siteNo ="02055000" )
USGS02054530=make_usgs_gage_list(siteNo ="02054530" )


#Data from DEM:
ab_ll=rbind(USGS02056000$site,
                USGS0205551460$site,
                USGS02055100$site,
                USGS02055000$site,
                USGS02054530$site)
class(ab_ll)
ab_ll@proj4string
proj4_utm = paste0("+proj=utm +zone=",
                     trunc((180+coordinates(USGS02055000$site)[1])/6+1), 
                     " +datum=WGS84 +units=m +no_defs")
print(proj4_utm)
# Lat/Lon (_ll) is much easier!
proj4_ll = "+proj=longlat"
crs_ll=CRS(proj4_ll)
crs_utm=CRS(proj4_utm)
proj4string(ab_ll)=proj4_ll
ab_utm=spTransform(ab_ll,crs_utm)
ab_utm@coords

#what changes here to change zoom extent?
mydem=get_aws_terrain(locations=ab_utm@coords, 
                        z = 12, prj = proj4_utm,expand=1)
#
# Lets plot the DEM and the gage locations so we can guess 
# what gages connect with what gages
#


pdf()
plot(mydem)
plot(ab_utm,add=T)
text(ab_utm, labels=ab_utm@data$site_no, cex=0.6, font=2,pos=1)
# Streams as USGS sees them, I know I can get an overview of streams with the 
# USGS H
#???
url="https://prd-tnm.s3.amazonaws.com/StagedProducts/Hydrography/NHD/HU8/Shape/NHD_H_03010101_HU8_Shape.zip"
curl_download(url,"NHD_H_03010101_HU8_Shape.zip")
unzip("NHD_H_03010101_HU8_Shape.zip",exdir="03010101")
streams=readOGR("03010101/Shape/NHDFlowline.dbf")
streams_utm=spTransform(streams,crs_utm)
plot(streams_utm,col="blue",add=T)
dev.off()

# #
# # Breaking things for educational purposes
# View(USGS02056000$flowdata)
# USGS02056000$flowdata=USGS02056000$flowdata[,c(1,2,3,4,5,8,10)]
# View(USGS02056000$flowdata)
# # Oh Noooooo!!!! This gage for some reason doesn't have "Gage height"
# # 00065! What can we do!?!? OK, no worries, we do have "Discharge" 00060	
###########################################
# 02056000 ROANOKE RIVER AT NIAGARA, VA
###########################################

# Assume in the inventory link that for this gage, our Gage height is missing. 
head(USGS02056000$flowdata,2)  # Note that we have 00060 but missing 00065...
#  agency_cd  site_no            dateTime X_00060_00000 X_00060_00000_cd tz_cd
#1      USGS 02056000 2017-05-01 04:00:00           876                A   UTC
#2      USGS 02056000 2017-05-01 04:15:00           876                A   UTC
# Hrm? No depth associated with the flow? BUT USGS maintains rating curves
# explain what a rating curve is: https://en.wikipedia.org/wiki/Rating_curve
# and use the readNWISrating() function to grab it for this gage

#no other gages need this method

USGS02056000[["rating"]]=readNWISrating(USGS02056000$site$site_no)
plot(USGS02056000$rating$DEP,USGS02056000$rating$INDEP,xlab="DEP",ylab="INDEP")
#

# Note that this is very similar to what we saw in the previous gage's results
# and as it turns out, we can use it to estimate a 00065 measurement as 
# we did for the previous gage.

#what is this doing to get zero?
USGS02056000$flowdata$X_00065_00000=approx(USGS02056000$rating$DEP,
                                             USGS02056000$rating$INDEP, xout = USGS02056000$flowdata$X_00060_00000, ties = mean)$y
points(USGS02056000$flowdata$X_00060_00000,USGS02056000$flowdata$X_00065_00000,
         col="red")
#
USGS02056000$flowdata$depth_m=USGS02056000$flowdata$X_00065_00000*0.3048
# m/ft depth
# depth created for USGS02056000

# A quick readthrough of the Example 1: Hiking around Maunga Whau
# in the package vignette. 
# vignette("Overview", package = "gdistance")
# Set the starting and ending locations
# determine the river reach length and slope using the gdistance package.
#


A=SpatialPoints(USGS0205551460$site)# Up gradient site Lick Run
B=SpatialPoints(USGS02056000$site) # Down gradient site ROA River atNiagara
proj4string(A)=proj4_ll
proj4string(B)=proj4_ll
A_utm=spTransform(A,crs_utm)
B_utm=spTransform(B,crs_utm)
# Cut the DEM down to a more manageable size
cropmydem=crop(mydem,extend(extent(ab_utm),600))
cropmydem=trim(cropmydem)
cropmydem=cropmydem*1000.0
plot(cropmydem)
plot(ab_utm,add=T)
# Set up the weighting functions
altDiff <- function(x){x[2] - x[1]}
hd <- transition(cropmydem, altDiff, 8, symm=FALSE)
slope <- geoCorrection(hd)
adj <- adjacent(cropmydem, cells=1:ncell(cropmydem), pairs=TRUE, directions=8)
speed <- slope
speed[adj] <- 6 * exp(-3.5 * abs(slope[adj] + 0.05))
Conductance <- geoCorrection(speed)
# Find and plot the flow path
AtoB <- shortestPath(Conductance, A_utm, B_utm, output="SpatialLines")
plot(AtoB,add=T)
plot(streams_utm,col="blue",add=T)
plot(AtoB,add=T)
SpatialLinesLengths(AtoB)


USGS0205551460$site$L=SpatialLinesLengths(AtoB) # km to m
USGS0205551460$site$L # reach length in m
#
#
# Getting slope, we will extract the slope for points A and B from the DEM and # divide the difference by the length in m, this gives us a much better 
# estimate of slope than taking the point slopes at the gage site
#
USGS0205551460$site$slope=(extract(mydem,A_utm)-
                               extract(mydem,B_utm))/USGS0205551460$site$L
USGS0205551460$site$slope

# So now we have flow depth (y "$depth_m"), manning's n ("$man_n"), Q ("$cms"), and slope ("$slope") rearrange to solve for B
# B=(n*Q)/(y^(5/3)*sqrt(So))
USGS0205551460$flowdata$B=(USGS0205551460$site$man_n*
                               USGS0205551460$flowdata$cms)/(USGS0205551460$flowdata$depth_m^(5/3)*
                                                               sqrt(USGS0205551460$site$slope))
head(USGS0205551460$flowdata)
#  agency_cd	site_no        	dateTime X_00060_00000 X_00060_00000_cd
#1  	USGS 05267000 2017-05-01 04:00:00      	6.38            	A
#2  	USGS 05267000 2017-05-01 04:05:00      	6.38            	A
#  X_00065_00000 X_00065_00000_cd tz_cd   	cms  depth_m    	B
#1      	2.74            	A   UTC 0.1806816 0.835152 0.103032
#2      	2.74            	A   UTC 0.1806816 0.835152 0.103032
#
# Lets look at how B changes with flow.    
plot(USGS0205551460$flowdata$dateTime,USGS0205551460$flowdata$B, main="LICK RUN TO ROANOKE RIVER AT NIAGARA, VA")
# Does this seem reasonable (...like order of magnitude reasonable)? You can 
# perform a quick and dirty check using google earth and measuring the channel 
# width in a few places.
#

#PLOT FIX IN HW Q1
# plot(USGS0205551460$flowdata$cms,USGS0205551460$flowdata$depth_m, main="LICK RUN TO ROANOKE RIVER AT NIAGARA, VA")
# 
# #adjust plot
# basedepth = min(USGS0205551460$flowdata$depth_m, na.rm=TRUE)
# USGS0205551460$flowdata$depth_m = USGS0205551460$flowdata$depth_m - basedepth
# baseflow = min(USGS0205551460$flowdata$cms, na.rm=TRUE)
# USGS0205551460$flowdata$cms = USGS0205551460$flowdata$cms - baseflow
# 
# pdf("~/2023/BSE5304Labs09/R/Lab09/adjustedgraph.pdf")
# plot(USGS0205551460$flowdata$cms,USGS0205551460$flowdata$depth_m, main="LICK RUN TO ROANOKE RIVER AT NIAGARA, VA")
# dev.off()
# ck
# USGS0205551460$flowdata$ck = 5/3*(sqrt(USGS0205551460$site$slope)/USGS0205551460$site$man_n)*USGS0205551460$flowdata$depth_m^(2/3)
#   # ANS
#   mean(USGS0205551460$flowdata$ck,na.rm=T)
# # [1] 2.547238 for this example, confirm this result
# USGS0205551460$flowdata$dt = USGS0205551460$site$L/USGS0205551460$flowdata$ck
#   mean(USGS0205551460$flowdata$dt,na.rm=T)
# # [1] 6328.655  for this example, confirm this result
# plot(USGS0205551460$flowdata$dateTime,USGS0205551460$flowdata$dt)
# 
# USGS0205551460$flowdata$outTime=USGS0205551460$flowdata$dateTime+
#   USGS0205551460$flowdata$dt

# Find the beginning of  Waves assuming a new wave starts at 110% of prior 
# flow. This might need to change for your homework



# Gage  USGS 0205551460 ------
# A=SpatialPoints(USGS0205551460$site)# Up gradient site Lick Run
# B=SpatialPoints(USGS02056000$site) # Down gradient site ROA River atNiagara
# proj4string(A)=proj4_ll
# proj4string(B)=proj4_ll
# A_utm=spTransform(A,crs_utm)
# B_utm=spTransform(B,crs_utm)
# # Cut the DEM down to a more manageable size
# cropmydem=crop(mydem,extend(extent(ab_utm),600))
# cropmydem=trim(cropmydem)
# cropmydem=cropmydem*1000.0
# plot(cropmydem)
# plot(ab_utm,add=T)
# # Set up the weighting functions
# altDiff <- function(x){x[2] - x[1]}
# hd <- transition(cropmydem, altDiff, 8, symm=FALSE)
# slope <- geoCorrection(hd)
# adj <- adjacent(cropmydem, cells=1:ncell(cropmydem), pairs=TRUE, directions=8)
# speed <- slope
# speed[adj] <- 6 * exp(-3.5 * abs(slope[adj] + 0.05))
# Conductance <- geoCorrection(speed)
# # Find and plot the flow path
# AtoB <- shortestPath(Conductance, A_utm, B_utm, output="SpatialLines")
# plot(AtoB,add=T)
# plot(streams_utm,col="blue",add=T)
# plot(AtoB,add=T)
# SpatialLinesLengths(AtoB)



  USGS0205551460$flowdata$ck = 5/3*(sqrt(USGS0205551460$site$slope)/USGS0205551460$site$man_n)*USGS0205551460$flowdata$depth_m^(2/3)
  USGS0205551460$flowdata$dt = USGS0205551460$site$L/USGS0205551460$flowdata$ck
  USGS0205551460$flowdata$outTime = USGS0205551460$flowdata$dateTime + USGS0205551460$flowdata$dt
  
  #how change this??
  WaveStartDecPercent=1.50 #where does the wave start?
  
  USGS0205551460$flowdata$newwave=
    USGS0205551460$flowdata$cms *WaveStartDecPercent <
    data.table::shift(USGS0205551460$flowdata$cms)
  summary(USGS0205551460$flowdata$newwave)

  # Add plot of the point found
  len=length(USGS0205551460$flowdata$newwave)
  USGS0205551460$flowdata$newwave[is.na(USGS0205551460$flowdata$newwave)]=F
  # Removes repeated finds by going through loop backwords
  for (i in seq(len,2)){
    print(i)
    if(USGS0205551460$flowdata$newwave[i]==T &
      USGS0205551460$flowdata$newwave[i-1]==T){
      USGS0205551460$flowdata$newwave[i]=F
    }
  }
  plot(USGS0205551460$flowdata$dateTime,USGS0205551460$flowdata$cms,type="l")
  points(USGS0205551460$flowdata$dateTime[USGS0205551460$flowdata$newwave],
         USGS0205551460$flowdata$cms[USGS0205551460$flowdata$newwave],col=2)

  # Find the time locations where waves begin
  which(USGS0205551460$flowdata$newwave == TRUE)
  plot(USGS0205551460$flowdata$dateTime,USGS0205551460$flowdata$cms,
       type="l",xlim=c(USGS0205551460$flowdata$dateTime[1109],
                       USGS0205551460$flowdata$dateTime[1109+200]))
  
  #out time not working??
  lines(USGS0205551460$flowdata$outTime,USGS0205551460$flowdata$cms,col=2)
  

  USGS0205551460$flowdata$ck = 5/3*(sqrt(USGS0205551460$site$slope)/USGS0205551460$site$man_n)*USGS0205551460$flowdata$depth_m^(2/3)
  USGS0205551460$flowdata$dt = USGS0205551460$site$L/USGS0205551460$flowdata$ck
  USGS0205551460$flowdata$outTime = USGS0205551460$flowdata$dateTime + USGS0205551460$flowdata$dt
  
  #how change this??
  WaveStartDecPercent=1.50 #where does the wave start?
  
  USGS0205551460$flowdata$newwave=
    USGS0205551460$flowdata$cms *WaveStartDecPercent <
    data.table::shift(USGS0205551460$flowdata$cms)
  summary(USGS0205551460$flowdata$newwave)

  # Add plot of the point found
  len=length(USGS0205551460$flowdata$newwave)
  USGS0205551460$flowdata$newwave[is.na(USGS0205551460$flowdata$newwave)]=F
  # Removes repeated finds by going through loop backwords
  for (i in seq(len,2)){
    print(i)
    if(USGS0205551460$flowdata$newwave[i]==T &
      USGS0205551460$flowdata$newwave[i-1]==T){
      USGS0205551460$flowdata$newwave[i]=F
    }
  }
  plot(USGS0205551460$flowdata$dateTime,USGS0205551460$flowdata$cms,type="l")
  points(USGS0205551460$flowdata$dateTime[USGS0205551460$flowdata$newwave],
         USGS0205551460$flowdata$cms[USGS0205551460$flowdata$newwave],col=2)

  # Find the time locations where waves begin
  which(USGS0205551460$flowdata$newwave == TRUE)
  pdf("~/2023/BSE5304Labs09/R/Lab09/gage0205551460.pdf")
  plot(USGS0205551460$flowdata$dateTime,USGS0205551460$flowdata$cms,
       type="l", xlim=c(USGS0205551460$flowdata$dateTime[1109],
                       USGS0205551460$flowdata$dateTime[1109+200]))
  
  lines(USGS0205551460$flowdata$outTime,USGS0205551460$flowdata$cms,col=2)
  dev.off()

# Gage 02055100 -----
  A=SpatialPoints(USGS02055100$site)# Up gradient site Lick Run
  B=SpatialPoints(USGS02056000$site) # Down gradient site ROA River atNiagara
  proj4string(A)=proj4_ll
  proj4string(B)=proj4_ll
  A_utm=spTransform(A,crs_utm)
  B_utm=spTransform(B,crs_utm)
  # Cut the DEM down to a more manageable size
  cropmydem=crop(mydem,extend(extent(ab_utm),600))
  cropmydem=trim(cropmydem)
  cropmydem=cropmydem*1000.0
  plot(cropmydem)
  plot(ab_utm,add=T)
  # Set up the weighting functions
  altDiff <- function(x){x[2] - x[1]}
  hd <- transition(cropmydem, altDiff, 8, symm=FALSE)
  slope <- geoCorrection(hd)
  adj <- adjacent(cropmydem, cells=1:ncell(cropmydem), pairs=TRUE, directions=8)
  speed <- slope
  speed[adj] <- 6 * exp(-3.5 * abs(slope[adj] + 0.05))
  Conductance <- geoCorrection(speed)
  # Find and plot the flow path
  AtoB <- shortestPath(Conductance, A_utm, B_utm, output="SpatialLines")
  plot(AtoB,add=T)
  plot(streams_utm,col="blue",add=T)
  plot(AtoB,add=T)
  SpatialLinesLengths(AtoB)
  
  
  
  
  USGS02055100$site$L=SpatialLinesLengths(AtoB) # km to m
  USGS02055100$site$L # reach length in m
  #
  #
  # Getting slope, we will extract the slope for points A and B from the DEM and # divide the difference by the length in m, this gives us a much better 
  # estimate of slope than taking the point slopes at the gage site
  #
  USGS02055100$site$slope=(extract(mydem,A_utm)-
                               extract(mydem,B_utm))/USGS02055100$site$L
  USGS02055100$site$slope
  
  # So now we have flow depth (y "$depth_m"), manning's n ("$man_n"), Q ("$cms"), and slope ("$slope") rearrange to solve for B
  # B=(n*Q)/(y^(5/3)*sqrt(So))
  USGS02055100$flowdata$B=(USGS02055100$site$man_n*
                               USGS02055100$flowdata$cms)/(USGS02055100$flowdata$depth_m^(5/3)*
                                                               sqrt(USGS02055100$site$slope))
  
  USGS02055100$flowdata$ck = 5/3*(sqrt(USGS02055100$site$slope)/USGS02055100$site$man_n)*USGS02055100$flowdata$depth_m^(2/3)
  USGS02055100$flowdata$dt = USGS02055100$site$L/USGS02055100$flowdata$ck
  USGS02055100$flowdata$outTime = USGS02055100$flowdata$dateTime + USGS02055100$flowdata$dt
  
  #how change this??
  WaveStartDecPercent=1 #where does the wave start?
  
  USGS02055100$flowdata$newwave=
    USGS02055100$flowdata$cms *WaveStartDecPercent <
    data.table::shift(USGS02055100$flowdata$cms)
  summary(USGS02055100$flowdata$newwave)
  
  # Add plot of the point found
  len=length(USGS02055100$flowdata$newwave)
  USGS02055100$flowdata$newwave[is.na(USGS02055100$flowdata$newwave)]=F
  # Removes repeated finds by going through loop backwords
  for (i in seq(len,2)){
    print(i)
    if(USGS02055100$flowdata$newwave[i]==T &
       USGS02055100$flowdata$newwave[i-1]==T){
      USGS02055100$flowdata$newwave[i]=F
    }
  }
  plot(USGS02055100$flowdata$dateTime,USGS02055100$flowdata$cms,type="l")
  points(USGS02055100$flowdata$dateTime[USGS02055100$flowdata$newwave],
         USGS02055100$flowdata$cms[USGS02055100$flowdata$newwave],col=2)
  
  # Find the time locations where waves begin
  which(USGS02055100$flowdata$newwave == TRUE)
  pdf("~/2023/BSE5304Labs09/R/Lab09/gage02055100.pdf")
  plot(USGS02055100$flowdata$dateTime,USGS02055100$flowdata$cms,
       type="l",xlim=c(USGS02055100$flowdata$dateTime[1109],
                       USGS02055100$flowdata$dateTime[1109+200]), ylim=c(0,3))
  
  #out time not working??
  lines(USGS02055100$flowdata$outTime,USGS02055100$flowdata$cms,col=2)
  dev.off()
  
  
  
# Gage 02054530  -----
  A=SpatialPoints(USGS02054530$site)# Up gradient site Lick Run
  B=SpatialPoints(USGS02056000$site) # Down gradient site ROA River atNiagara
  proj4string(A)=proj4_ll
  proj4string(B)=proj4_ll
  A_utm=spTransform(A,crs_utm)
  B_utm=spTransform(B,crs_utm)
  # Cut the DEM down to a more manageable size
  cropmydem=crop(mydem,extend(extent(ab_utm),600))
  cropmydem=trim(cropmydem)
  cropmydem=cropmydem*1000.0
  plot(cropmydem)
  plot(ab_utm,add=T)
  # Set up the weighting functions
  altDiff <- function(x){x[2] - x[1]}
  hd <- transition(cropmydem, altDiff, 8, symm=FALSE)
  slope <- geoCorrection(hd)
  adj <- adjacent(cropmydem, cells=1:ncell(cropmydem), pairs=TRUE, directions=8)
  speed <- slope
  speed[adj] <- 6 * exp(-3.5 * abs(slope[adj] + 0.05))
  Conductance <- geoCorrection(speed)
  # Find and plot the flow path
  AtoB <- shortestPath(Conductance, A_utm, B_utm, output="SpatialLines")
  plot(AtoB,add=T)
  plot(streams_utm,col="blue",add=T)
  plot(AtoB,add=T)
  SpatialLinesLengths(AtoB)
  
  
  
  
  USGS02054530$site$L=SpatialLinesLengths(AtoB) # km to m
  USGS02054530$site$L # reach length in m
  #
  #
  # Getting slope, we will extract the slope for points A and B from the DEM and # divide the difference by the length in m, this gives us a much better 
  # estimate of slope than taking the point slopes at the gage site
  #
  USGS02054530$site$slope=(extract(mydem,A_utm)-
                               extract(mydem,B_utm))/USGS02054530$site$L
  USGS02054530$site$slope
  
  # So now we have flow depth (y "$depth_m"), manning's n ("$man_n"), Q ("$cms"), and slope ("$slope") rearrange to solve for B
  # B=(n*Q)/(y^(5/3)*sqrt(So))
  USGS02054530$flowdata$B=(USGS02054530$site$man_n*
                               USGS02054530$flowdata$cms)/(USGS02054530$flowdata$depth_m^(5/3)*
                                                               sqrt(USGS02054530$site$slope))
  
  USGS02054530$flowdata$ck = 5/3*(sqrt(USGS02054530$site$slope)/USGS02054530$site$man_n)*USGS02054530$flowdata$depth_m^(2/3)
  USGS02054530$flowdata$dt = USGS02054530$site$L/USGS02054530$flowdata$ck
  USGS02054530$flowdata$outTime = USGS02054530$flowdata$dateTime + USGS02054530$flowdata$dt
  
  #how change this??
  WaveStartDecPercent=.8 #where does the wave start?
  
  USGS02054530$flowdata$newwave=
    USGS02054530$flowdata$cms *WaveStartDecPercent <
    data.table::shift(USGS02054530$flowdata$cms)
  summary(USGS02054530$flowdata$newwave)
  
  # Add plot of the point found
  len=length(USGS02054530$flowdata$newwave)
  USGS02054530$flowdata$newwave[is.na(USGS02054530$flowdata$newwave)]=F
  # Removes repeated finds by going through loop backwords
  for (i in seq(len,2)){
    print(i)
    if(USGS02054530$flowdata$newwave[i]==T &
       USGS02054530$flowdata$newwave[i-1]==T){
      USGS02054530$flowdata$newwave[i]=F
    }
  }
  plot(USGS02054530$flowdata$dateTime,USGS02054530$flowdata$cms,type="l")
  points(USGS02054530$flowdata$dateTime[USGS02054530$flowdata$newwave],
         USGS02054530$flowdata$cms[USGS02054530$flowdata$newwave],col=2)
  
  # Find the time locations where waves begin
  which(USGS02054530$flowdata$newwave == TRUE)
  pdf("~/2023/BSE5304Labs09/R/Lab09/gage02054530.pdf")
  plot(USGS02054530$flowdata$dateTime,USGS02054530$flowdata$cms,
       type="l",xlim=c(USGS02054530$flowdata$dateTime[1109],
                       USGS02054530$flowdata$dateTime[1109+200]))
  
  #out time not working??
  lines(USGS02054530$flowdata$outTime,USGS02054530$flowdata$cms,col=2)
  dev.off()
  
# Gage USGS 02055000 -----
  
  A=SpatialPoints(USGS02055000$site)# Up gradient site Lick Run
  B=SpatialPoints(USGS02056000$site) # Down gradient site ROA River atNiagara
  proj4string(A)=proj4_ll
  proj4string(B)=proj4_ll
  A_utm=spTransform(A,crs_utm)
  B_utm=spTransform(B,crs_utm)
  # Cut the DEM down to a more manageable size
  cropmydem=crop(mydem,extend(extent(ab_utm),600))
  cropmydem=trim(cropmydem)
  cropmydem=cropmydem*1000.0
  plot(cropmydem)
  plot(ab_utm,add=T)
  # Set up the weighting functions
  altDiff <- function(x){x[2] - x[1]}
  hd <- transition(cropmydem, altDiff, 8, symm=FALSE)
  slope <- geoCorrection(hd)
  adj <- adjacent(cropmydem, cells=1:ncell(cropmydem), pairs=TRUE, directions=8)
  speed <- slope
  speed[adj] <- 6 * exp(-3.5 * abs(slope[adj] + 0.05))
  Conductance <- geoCorrection(speed)
  # Find and plot the flow path
  AtoB <- shortestPath(Conductance, A_utm, B_utm, output="SpatialLines")
  plot(AtoB,add=T)
  plot(streams_utm,col="blue",add=T)
  plot(AtoB,add=T)
  SpatialLinesLengths(AtoB)
  
  
  
  
  USGS02055000$site$L=SpatialLinesLengths(AtoB) # km to m
  USGS02055000$site$L # reach length in m
  #
  #
  # Getting slope, we will extract the slope for points A and B from the DEM and # divide the difference by the length in m, this gives us a much better 
  # estimate of slope than taking the point slopes at the gage site
  #
  USGS02055000$site$slope=(extract(mydem,A_utm)-
                             extract(mydem,B_utm))/USGS02055000$site$L
  USGS02055000$site$slope
  
  # So now we have flow depth (y "$depth_m"), manning's n ("$man_n"), Q ("$cms"), and slope ("$slope") rearrange to solve for B
  # B=(n*Q)/(y^(5/3)*sqrt(So))
  USGS02055000$flowdata$B=(USGS02055000$site$man_n*
                             USGS02055000$flowdata$cms)/(USGS02055000$flowdata$depth_m^(5/3)*
                                                           sqrt(USGS02055000$site$slope))
  
  USGS02055000$flowdata$ck = 5/3*(sqrt(USGS02055000$site$slope)/USGS02055000$site$man_n)*USGS02055000$flowdata$depth_m^(2/3)
  USGS02055000$flowdata$dt = USGS02055000$site$L/USGS02055000$flowdata$ck
  USGS02055000$flowdata$outTime = USGS02055000$flowdata$dateTime + USGS02055000$flowdata$dt
  
  #how change this??
  WaveStartDecPercent=1.10 #where does the wave start?
  
  USGS02055000$flowdata$newwave=
    USGS02055000$flowdata$cms *WaveStartDecPercent <
    data.table::shift(USGS02055000$flowdata$cms)
  summary(USGS02055000$flowdata$newwave)
  
  # Add plot of the point found
  len=length(USGS02055000$flowdata$newwave)
  USGS02055000$flowdata$newwave[is.na(USGS02055000$flowdata$newwave)]=F
  # Removes repeated finds by going through loop backwords
  for (i in seq(len,2)){
    print(i)
    if(USGS02055000$flowdata$newwave[i]==T &
       USGS02055000$flowdata$newwave[i-1]==T){
      USGS02055000$flowdata$newwave[i]=F
    }
  }
  plot(USGS02055000$flowdata$dateTime,USGS02055000$flowdata$cms,type="l")
  points(USGS02055000$flowdata$dateTime[USGS02055000$flowdata$newwave],
         USGS02055000$flowdata$cms[USGS02055000$flowdata$newwave],col=2)
  
  # Find the time locations where waves begin
  which(USGS02055000$flowdata$newwave == TRUE)
  pdf("~/2023/BSE5304Labs09/R/Lab09/gage02055000.pdf")
  plot(USGS02055000$flowdata$dateTime,USGS02055000$flowdata$cms,
       type="l",xlim=c(USGS02055000$flowdata$dateTime[1109],
                       USGS02055000$flowdata$dateTime[1109+200]))
  
  #out time not working??
  lines(USGS02055000$flowdata$outTime,USGS02055000$flowdata$cms,col=2)
  dev.off()

  
# MY OWN RIVER SYSTEM - MONONGAHELA ------------

#upstream
  
# 03072655 Monongahela River near Masontown, PA
# 03075070 Monongahela River at Elizabeth, PA
# 03086000 Ohio River at Sewickley, PA

#downstream

USGS03072655=make_usgs_gage_list(siteNo = "03072655")
USGS03075070=make_usgs_gage_list(siteNo = "03075070")
USGS03086000=make_usgs_gage_list(siteNo ="03086000")


ab_ll=rbind(USGS03072655$site,
            USGS03075070$site,
            USGS03086000$site)
class(ab_ll)
ab_ll@proj4string
proj4_utm = paste0("+proj=utm +zone=",
                   trunc((180+coordinates(USGS03075070$site)[1])/6+1),
                   " +datum=WGS84 +units=m +no_defs")
print(proj4_utm)
# Lat/Lon (_ll) is much easier!
proj4_ll = "+proj=longlat"
crs_ll=CRS(proj4_ll)
crs_utm=CRS(proj4_utm)
proj4string(ab_ll)=proj4_ll
ab_utm=spTransform(ab_ll,crs_utm)
ab_utm@coords


mydem=get_aws_terrain(locations=ab_utm@coords,
                      z = 10, prj = proj4_utm,expand=1)
pdf("~/2023/BSE5304Labs09/R/Lab09/dem.pdf")
plot(mydem)
plot(ab_utm,add=T)
text(ab_utm, labels=ab_utm@data$site_no, cex=0.6, font=2,pos=1)
dev.off()

A=SpatialPoints(USGS03072655$site)# Up gradient site Lick Run
B=SpatialPoints(USGS03086000$site) # Down gradient site ROA River atNiagara
proj4string(A)=proj4_ll
proj4string(B)=proj4_ll
A_utm=spTransform(A,crs_utm)
B_utm=spTransform(B,crs_utm)
# Cut the DEM down to a more manageable size
#change?
cropmydem=crop(mydem,extend(extent(ab_utm),2000))
cropmydem=trim(cropmydem)
cropmydem=cropmydem*1000.0

pdf("~/2023/BSE5304Labs09/R/Lab09/cropped_dem.pdf")
plot(cropmydem)
plot(ab_utm,add=T)
dev.off()
# Set up the weighting functions
altDiff <- function(x){x[2] - x[1]}
hd <- transition(cropmydem, altDiff, 8, symm=FALSE)
slope <- geoCorrection(hd)
adj <- adjacent(cropmydem, cells=1:ncell(cropmydem), pairs=TRUE, directions=8)
speed <- slope
speed[adj] <- 6 * exp(-3.5 * abs(slope[adj] + 0.05))
Conductance <- geoCorrection(speed)
# Find and plot the flow path
AtoB <- shortestPath(Conductance, A_utm, B_utm, output="SpatialLines")
pdf("~/2023/BSE5304Labs09/R/Lab09/dem_streams_wrong.pdf")
plot(cropmydem)
plot(ab_utm,add=T)
plot(AtoB,add=T)
plot(streams_utm,col="blue",add=T)
plot(AtoB,add=T)
SpatialLinesLengths(AtoB)
dev.off()


#adjust spatial points
A=SpatialPoints(USGS03075070$site)# Up gradient site Lick Run
B=SpatialPoints(USGS03086000$site) # Down gradient site ROA River atNiagara
proj4string(A)=proj4_ll
proj4string(B)=proj4_ll
A_utm=spTransform(A,crs_utm)
B_utm=spTransform(B,crs_utm)
# Cut the DEM down to a more manageable size
#change?
cropmydem=crop(mydem,extend(extent(ab_utm),600))
cropmydem=trim(cropmydem)
cropmydem=cropmydem*1000.0

# pdf("~/2023/BSE5304Labs09/R/Lab09/cropped_dem.pdf")
plot(cropmydem)
plot(ab_utm,add=T)
# dev.off()
# Set up the weighting functions
altDiff <- function(x){x[2] - x[1]}
hd <- transition(cropmydem, altDiff, 8, symm=FALSE)
slope <- geoCorrection(hd)
adj <- adjacent(cropmydem, cells=1:ncell(cropmydem), pairs=TRUE, directions=8)
speed <- slope
speed[adj] <- 6 * exp(-3.5 * abs(slope[adj] + 0.05))
Conductance <- geoCorrection(speed)

pdf("~/2023/BSE5304Labs09/R/Lab09/dem_streams.pdf")
plot(cropmydem)
plot(ab_utm,add=T)
plot(AtoB,add=T)
plot(streams_utm,col="blue",add=T)
plot(AtoB,add=T)
SpatialLinesLengths(AtoB)
dev.off()

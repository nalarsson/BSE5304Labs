# 
# Since everything depends on the libraries you install
# it is worthwhile loading them at the beginning
#
objects()  # This will list the objects you have.
rm(list=objects()) # Removes ALL the objects… so be careful here.
#
# What is going to change from use case to use case 
LabNo="/Lab04a"
myflowgage_id="0205551460"  # Old Friendly Gage
#
# What needs to be loaded
#
if (!require("pacman")) install.packages("pacman")
pacman::p_load(ggplot2,dplyr,patchwork,rnoaa)
pacman::p_load(operators,topmodel,DEoptim,soilDB,sp,curl,httr,
               rnoaa,raster,shapefiles,rgdal,elevatr,terra,progress,lubridate)
#
# Getting our organization on for where we want to put
# Data, external programs, and our project files.
# Things are going to get messy if we don't start issolating
# our data files by Lab
#
myhomedir=Sys.getenv("HOME")
datadir=paste0(myhomedir,"/data",LabNo)
dir.create(datadir,recursive = T)
srcdir=paste0(myhomedir,"/src")
dir.create(srcdir,recursive = T)
# Setting the directory for where the GitHub project exists. 
# This depends on where you set up your git, and what you called it locally, 
# but when you start a new git project, it will be the first directory you 
# are placed in... or if later in the project:
# WOOO HOOO... took me a few hours to find this function!
# 
mygitdir=rstudioapi::getActiveProject()
mypdfdir=paste0(mygitdir,"/pdfs",LabNo)
dir.create(mypdfdir)
# 
setwd(mygitdir)
system("git config --global user.email 'nlarsson@vt.edu' ") 
system("git config --global user.name 'Natalie Larsson' ")
system("git config pull.rebase false")
#how did I blow this up???
#
# This week, we discovered some "features" that make removing and 
# re-installing the EcoHydrology Library necessary.
#
setwd(srcdir)
#detach("package:EcoHydRology", unload = TRUE)
# remove.packages("EcoHydRology", lib="~/R/x86_64-pc-linux-gnu-library/4.2")
system("svn checkout svn://scm.r-forge.r-project.org/svnroot/ecohydrology/"); 
install.packages(c("ecohydrology/pkg/EcoHydRology/"),repos = NULL)
pacman::p_load(EcoHydRology)

setwd(datadir)
#
# Should we do a gage that is easy, or deal with some reality?
#
myflowgage=get_usgs_gage(myflowgage_id,begin_date = "2015-01-01",
                         end_date = "2022-03-01")

#
# This is where some folks had issues... they forgot to check their 
# watershed areas per the homework... though there were ways to fix
# it later with lower resolution DEM pull
#
print(paste0("reported Area ",myflowgage$area))
# For most watershed modelling purposes we normalize Q in mm/day for basins
myflowgage$flowdata$Qmm = myflowgage$flowdata$flow/myflowgage$area/10^3

# In the Lab02, we introduced you to a way to quickly get your WX Data 
# for any location in the world way easier than traditional download and
# parsing methods most old people use.
#
source("https://raw.githubusercontent.com/Rojakaveh/FillMissWX/main/FillMissWX.R")
# Remove this if it syncs with the EcoHydrology Version.
WXData=FillMissWX(declat=myflowgage$declat, declon=myflowgage$declon,
                  StnRadius=30,minstns=10,date_min="2010-01-01",
                  date_max="2023-02-01",targElev=myflowgage$elev,
                  method = "IDW",alfa=2)

BasinData=merge(WXData,myflowgage$flowdata,by.x="date",by.y="mdate")
TMWB=BasinData



# Things we repeat can be saved and accessed without messing up our working area
# https://raw.githubusercontent.com/vtdrfuka/BSE5304Labs/main/R/Lab05SetupDRF.R
# becomes: 
url="https://raw.githubusercontent.com/vtdrfuka/BSE5304Labs/main/R/Lab05SetupDRF.R"
# This will grab the solution for last weeks Lab03 Homework
download.file(url,"Lab05SetupDRF.R")
file.edit("Lab05SetupDRF.R")

# Grab out models for Snow and TMWB
# https://raw.githubusercontent.com/vtdrfuka/BSE5304Labs/main/R/TMWBFuncs.R
# becomes: 
url="https://raw.githubusercontent.com/vtdrfuka/BSE5304Labs/main/R/TMWBFuncs.R"
# This will grab the solution for last weeks Lab03 Homework
download.file(url,"TMWBFuncs.R")
#file.edit("TMWBFuncs.R")
source(url)
# I actually am starting to trust my snow model
url="https://raw.githubusercontent.com/vtdrfuka/BSE5304Labs/main/R/TISnow.R"
# This will grab the solution for last weeks Lab03 Homework
source(url)

outTMWB = TMWBmodel(TMWBdf=TMWB)
NSE(Yobs=outTMWB$Qmm, Ysim=outTMWB$Qpred)

TMWBoptFunc <- function(x){
  x1 <- x[1]
  x2 <- x[2]
  x3 <- x[3]
  x4 <- x[4]
  outTMWB=TMWBmodel(TMWBdf = TMWB, fcres = x1, Z=x2, SFT=x3, bmlt6=x4)
  1-NSE(Yobs=outTMWB$Qmm, Ysim=outTMWB$Qpred)
}
x=c(0.3,1000)
TMWBoptFunc(x)
lower <- c(.01,300,-1,.1)
upper <-c(.95,3000,6,5)

#outDEMoptim=DEoptim(TMWBoptFunc,lower,upper)
#detach issues??
outDEoptim <- DEoptim(TMWBoptFunc, lower, upper, 
                      DEoptim.control(NP = 80, #80 iteration running 80*400 with different numbers
                                      itermax = 10, 
                                      F = 1.2, CR = 0.7))





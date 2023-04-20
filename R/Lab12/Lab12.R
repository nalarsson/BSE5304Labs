## INTRO AND EXAMPLES -----

if(!exists("LabNo")){
  LabNum <- as.numeric(readline(prompt="What Week is this? "))
  LabNo=paste0("/Lab",LabNum)
}

# Since everything depends on the libraries you install
# it is worthwhile loading them at the beginning
#
if (!require("pacman")) install.packages("pacman")
pacman::p_load(ggplot2,dplyr,patchwork,rnoaa)
pacman::p_load(operators,topmodel,DEoptim,soilDB,sp,curl,httr,
               rnoaa,raster,shapefiles,rgdal,elevatr,terra,progress,lubridate)
#
# Getting our organization on for where we want to put
# Data, external programs, and our project files.
# Things are going to get messy if we don't start isolating
# our data files by Lab
#
user=Sys.getenv("USER")
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
system(paste0("git config --global user.email 'nlarsson@vt.edu'")) 
system(paste0("git config --global user.name 'Natalie Larsson' "))
system("git config pull.rebase false")
#
# This week, we discovered some "features" that make removing and 
# re-installing the EcoHydrology Library necessary.
#
#setwd(srcdir)

#detach("package:EcoHydRology", unload = TRUE)
#remove.packages("EcoHydRology", lib="~/R/x86_64-pc-linux-gnu-library/4.2")
#system("svn checkout svn://scm.r-forge.r-project.org/svnroot/ecohydrology/"); 
#install.packages(c("ecohydrology/pkg/EcoHydRology/"),repos = NULL)
pacman::p_load(EcoHydRology)
setwd(datadir)



pacman::p_load(data.table,multisensi)
# KEEP OPEN AS YOU WILL BE WALKING THROUGH IT FOR LAB	
vignette("multisensi-vignette")

verhulst <- function(K, Y0, a, t) {
  output <- K/(1 + (K/Y0 - 1) * exp(-a * t))
  return(output)
}

T <- seq(from = 5, to = 100, by = 5)
verhulst2 <- function(X, t = T) {
  out <- matrix(nrow = nrow(X), ncol = length(t), NA)
  for (i in 1:nrow(X)) {
    out[i, ] <- verhulst(X$K[i], X$Y0[i], X$a[i], t)
  }
  out <- as.data.frame(out)
  names(out) <- paste("t", t, sep = "")
  return(out)
}

n <- 10
set.seed(1234)

X <- data.frame(K = runif(n, min = 100, max = 1000), Y0 = runif(n, min = 1,
                                                                max = 40), a = runif(n, min = 0.05, max = 0.2))
Y <- verhulst2(X)
par(cex.axis = 0.7, cex.lab = 0.8)
plot(T, Y[1, ], type = "l", xlab = "Time", ylab = "Population size",
     ylim = c(0, 1000))
for (i in 2:n) {
  lines(T, Y[i, ], type = "l", col = i)
}

library(multisensi)
verhulst.seq <- multisensi(model=verhulst2, reduction=NULL, center=FALSE,
                           design.args = list( K=c(100,400,1000), Y0=c(1,20,40), a=c(0.05,0.1,0.2)))

print(verhulst.seq, digits = 2)

plot(verhulst.seq, normalized = TRUE, color = terrain.colors, gsi.plot = FALSE)
title(xlab = "Time in half-decades.")
plot(verhulst.seq, normalized = FALSE, color = terrain.colors, gsi.plot = FALSE)
title(xlab = "Time in half-decades.")

X <- expand.grid(K=c(100,400,1000), Y0=c(1,20,40), a=c(0.05,0.1,0.2))
Y <- verhulst2(X) ## this part can be performed outside R if necessary
verhulst.seq <- multisensi(design=X, model=Y, reduction=NULL, center=FALSE)

verhulst.pca <- multisensi(design=X, model=Y, reduction=basis.ACP, scale=FALSE)

summary(verhulst.pca, digits = 2)

plot(verhulst.pca, graph = 1)

plot(verhulst.pca, graph = 2)

plot(verhulst.pca, graph = 3)

verhulst.poly <- multisensi(design = X, model = Y, reduction = basis.poly,
                            dimension = 0.99, center = FALSE, scale = FALSE, cumul = FALSE,
                            basis.args = list(degree=6, x.coord=T), analysis = analysis.anoasg,
                            analysis.args = list(formula=2, keep.outputs=FALSE))

summary(verhulst.poly, digits=2)

plot(verhulst.poly, nb.comp = 3, graph = 1)

## bsplines
verhulst.bspl <- multisensi(design=X, model=Y, reduction=basis.bsplines,
                            dimension=NULL, center=FALSE, scale=FALSE,
                            basis.args=list(knots=10, mdegree=3), cumul=FALSE,
                            analysis=analysis.anoasg,
                            analysis.args=list(formula=2, keep.outputs=FALSE))

plot(verhulst.bspl, nb.comp = 5, graph = 1)

library(sensitivity)

m <- 10000
Xb <- data.frame(K = runif(m, min = 100, max = 1000), Y0 = runif(m, min = 1,
                                                                 max = 40), a = runif(m, min = 0.05, max = 0.2))
verhulst.seq.sobol <- multisensi(design = sobol2007, model = verhulst2,
                                 reduction = NULL, analysis = analysis.sensitivity, center = TRUE,
                                 design.args = list(X1 = Xb[1:(m/2), ], X2 = Xb[(1 + m/2):m, ], nboot = 100),
                                 analysis.args = list(keep.outputs = FALSE))
## [*] Design
## [*] Response simulation
## [*] Analysis + Sensitivity Indices
print(verhulst.seq.sobol, digits = 2)
plot(verhulst.seq.sobol, normalized = TRUE, color = terrain.colors, gsi.plot = FALSE)
title(xlab = "Time in half-decades")





verhulst.seq.fast <- multisensi(design = fast99, model = verhulst2,
                                center = FALSE, reduction = NULL, analysis = analysis.sensitivity,
                                design.args=list( factors=c("K","Y0","a"), n=1000, q = "qunif",
                                                  q.arg = list(list(min=100, max=1000), list(min=1, max=40),
                                                               list(min = 0.05, max = 0.2))),
                                analysis.args=list(keep.outputs=FALSE))

print(verhulst.seq.fast,digits=2)

plot(verhulst.seq.fast, normalized = TRUE, color = terrain.colors, gsi.plot = FALSE)
title(xlab = "Time in half-decades")


## EcoHydrology Function --------


J <- seq(from = 1, to = 365, by = 5)
# Solar(lat, Jday, Tx, Tn, albedo=0.2, forest=0, slope=0, aspect = 0,
#      units="kJm2d")
# Note that the EcoHydRology::Solar() function is for specific days, 
# as such, we will want to create a function to loop through our period
# of interest:
Solar_Looped <- function(X, Jday = J) {
  out <- matrix(nrow = nrow(X), ncol = length(Jday), NA)
  for (i in 1:nrow(X)) {
    out[i, ] <- Solar(lat=X$lat[i],
                      Jday=Jday, Tx=X$Tx[i], 
                      Tn=(X$Tx[i]-X$Trange[i]), 
                      X$slope[i],X$aspect[i],units="Wm2")
  }
  out <- as.data.frame(out)
  names(out) <- paste("Jday", Jday, sep = "")
  return(out)
}

n <- 10
set.seed(1234)
X <- data.frame(Tx = runif(n, min = 5, max = 30), Trange = runif(n, min = 2,
                                                                   max = 16), slope = runif(n, min = 0.0, max = 0.2),
                  aspect = runif(n, min = 0.0, max = 0.2),
                  lat=runif(n, min = 0.0, max = 1.1))  # 1.1 radians lat is where?

# View(X)
#
Y <- Solar_Looped(X,Jday = J)
#
# You can ignore all the warnings, remember Errors=bad, warnings=not so much 
# So lets move on and build our summary graph
par(cex.axis = 0.7, cex.lab = 0.8)
plot(J, Y[1, ], type = "l", xlab = "Day of Year", 
       ylab = "Surface Short Wave Rad(W/m^2)")
for (i in 2:n) {
  lines(J, Y[i, ], type = "l", col = i)
}
# 
# Well, that is kewl, yet expected
#
# Multisensitivities
# 3 Sequential univariate sensitivity analyses
# 3.1 Calculation of sensitivity indices
Solar_Looped.seq <- multisensi(model=Solar_Looped, reduction=NULL, center=FALSE,
                                 design.args = list( Tx = c(5,15,25), 
                                                     Trange = c(2,9,16), 
                                                     slope = c(0.1,0.2,0.3),
                                                     aspect = c(0.1,.5,1.0),
                                                     lat=c(0.1,.77,1.1)))

print(Solar_Looped.seq, digits = 2)
#
# 3.2 Graphical representation of sensitivity indices
#
# dev.off() # Clean up previous par()
plot(Solar_Looped.seq, normalized = TRUE, color = terrain.colors, gsi.plot = FALSE)#normalized the upper subplot shows the extreme (tirets), #inter-quartile (grey) and median (bold line) output values
title(xlab = "Days of the Year.")
plot(Solar_Looped.seq, normalized = FALSE, color = terrain.colors, gsi.plot = FALSE)
title(xlab = "Days of the Year.")




# 4 Multivariate sensitivity analysis based on PCA
# read the vignette, though note we are using the multisensi() function
# to run our model (i.e. no “design” variable, and model=Solar_Looped)
Solar_Looped.pca <- multisensi(model=Solar_Looped, reduction=basis.ACP, scale=FALSE,
                                 design.args = list( Tx = c(5,15,25), 
                                                     Trange = c(2,9,16), 
                                                     slope = c(0.1,0.2,0.3),
                                                     aspect = c(0.1,.5,1.0),
                                                     lat=c(0.1,.77,1.1)))

summary(Solar_Looped.pca, digits = 2)
# 4.2 Graphical representation for PCA based analysis with 
# explanation in vignette. These graphs require the plot window to be larger
# and might give "Error in plot.new() : figure margins too large". 
# If so expand the plot window.
# dev.off()
plot(Solar_Looped.pca, graph = 1)
plot(Solar_Looped.pca, graph = 2)
plot(Solar_Looped.pca, graph = 3)
#
# 5.1 Polynomial reduction of the multivariate output
# Skip 5.1 Polynomial reduction for now and move on to
# 6 Alternative methods of sensitivity analysis
# 6.1 With Sobol2007 implemented in the package sensitivity
# 
library(sensitivity)
m <- 10000
Xb <- data.frame(Tx = runif(m, min = 5, max = 30), 
                   Trange = runif(m, min = 2,max = 16), 
                   slope = runif(m, min = 0.0, max = 0.2),
                   aspect = runif(m, min = 0.0, max = 0.2),
                   lat=runif(m, min = 0.0, max = 1.1))

Solar_Looped.seq.sobol <- multisensi(design = sobol2007, model = Solar_Looped,
                                       reduction = NULL, analysis = analysis.sensitivity, center = TRUE,
                                       design.args = list(X1 = Xb[1:(m/2), ], X2 = Xb[(1 + m/2):m, ], nboot = 100),
                                       analysis.args = list(keep.outputs = FALSE))
#
# Note, this is a good time time get a drink of water and/or pee as 
# it is running the function m=10,000 times (a few minutes).
#
print(Solar_Looped.seq.sobol, digits = 2)
# dev.off()
plot(Solar_Looped.seq.sobol, normalized = TRUE, color = terrain.colors)


## HOMEWORK QUESTION 1  - NETRAD -------

J <- seq(from = 1, to = 365, by = 5)

NetRad_Looped <- function(X, Jday = J) {
  out <- matrix(nrow = nrow(X), ncol = length(Jday), NA)
  for (i in 1:nrow(X)) {
    out[i, ] <- NetRad( lat=X$lat[i],
                              Jday=Jday, Tx=X$Tx[i],
                              Tn=(X$Tx[i]-X$Trange[i]),
                              X$slope[i],X$aspect[i],units="Wm2"
    )
  }
  out <- as.data.frame(out)
  names(out) <- paste("Jday", Jday, sep = "")
  return(out)
}



n <- 10
set.seed(1234)
X <- data.frame(Tx = runif(n, min = 5, max = 30), Trange = runif(n, min = 2,
                                                                 max = 16), slope = runif(n, min = 0.0, max = 0.2),
                aspect = runif(n, min = 0.0, max = 0.2),
                lat=runif(n, min = 0.0, max = 1.1))  # 1.1 radians lat is where?


Y <- NetRad_Looped(X,Jday = J)
#
# You can ignore all the warnings, remember Errors=bad, warnings=not so much 
# So lets move on and build our summary graph
par(cex.axis = 0.7, cex.lab = 0.8)
plot(J, Y[1, ], type = "l", xlab = "Day of Year", 
     ylab = "Net Surface Short Wave Rad(W/m^2)")
for (i in 2:n) {
  lines(J, Y[i, ], type = "l", col = i)
}


NetRad_Looped.seq <- multisensi(model=NetRad_Looped, reduction=NULL, center=FALSE,
                               design.args = list( Tx = c(5,15,25), 
                                                   Trange = c(2,9,16), 
                                                   slope = c(0.1,0.2,0.3),
                                                   aspect = c(0.1,.5,1.0),
                                                   lat=c(0.1,.77,1.1)))

print(NetRad_Looped.seq, digits = 2)
#
# 3.2 Graphical representation of sensitivity indices
#
# dev.off() # Clean up previous par()
plot(NetRad_Looped.seq, normalized = TRUE, color = terrain.colors, gsi.plot = FALSE)#normalized the upper subplot shows the extreme (tirets), #inter-quartile (grey) and median (bold line) output values
title("NetRad Sensitivity Analysis", xlab = "Days of the Year.")
plot(NetRad_Looped.seq, normalized = FALSE, color = terrain.colors, gsi.plot = FALSE)
title("NetRad Sensitivity Analysis", xlab = "Days of the Year.")




# 4 Multivariate sensitivity analysis based on PCA
# read the vignette, though note we are using the multisensi() function
# to run our model (i.e. no “design” variable, and model=Solar_Looped)
NetRad_Looped.pca <- multisensi(model=NetRad_Looped, reduction=basis.ACP, scale=FALSE,
                               design.args = list( Tx = c(5,15,25), 
                                                   Trange = c(2,9,16), 
                                                   slope = c(0.1,0.2,0.3),
                                                   aspect = c(0.1,.5,1.0),
                                                   lat=c(0.1,.77,1.1)))

summary(NetRad_Looped.pca, digits = 2)
# 4.2 Graphical representation for PCA based analysis with 
# explanation in vignette. These graphs require the plot window to be larger
# and might give "Error in plot.new() : figure margins too large". 
# If so expand the plot window.
# dev.off()
plot(NetRad_Looped.pca, graph = 1)
plot(NetRad_Looped.pca, graph = 2)
title("Most Sensitive Variables, NetRad")
plot(NetRad_Looped.pca, graph = 3)
#
# 5.1 Polynomial reduction of the multivariate output
# Skip 5.1 Polynomial reduction for now and move on to
# 6 Alternative methods of sensitivity analysis
# 6.1 With Sobol2007 implemented in the package sensitivity
# 
m <- 10000
Xb <- data.frame(Tx = runif(m, min = 5, max = 30), 
                 Trange = runif(m, min = 2,max = 16), 
                 slope = runif(m, min = 0.0, max = 0.2),
                 aspect = runif(m, min = 0.0, max = 0.2),
                 lat=runif(m, min = 0.0, max = 1.1))

NetRad_Looped.seq.sobol <- multisensi(design = sobol2007, model = NetRad_Looped,
                                     reduction = NULL, analysis = analysis.sensitivity, center = TRUE,
                                     design.args = list(X1 = Xb[1:(m/2), ], X2 = Xb[(1 + m/2):m, ], nboot = 100),
                                     analysis.args = list(keep.outputs = FALSE))
#
# Note, this is a good time time get a drink of water and/or pee as 
# it is running the function m=10,000 times (a few minutes).
#
print(NetRad_Looped.seq.sobol, digits = 2)
# dev.off()
plot(NetRad_Looped.seq.sobol, normalized = TRUE, color = terrain.colors)




## HOMEWORK QUESTION 1 - PET_FROMTEMP -----
# Solution for PET_fromTemp
# Trick is you have to notice that "lat_radians" has replaced "lat" and
# there is no "units" variable... and... notice that the function has to
# be fixed to allow Jday to be a vector of a different size than Tmax and Tmin
PET_fromTemp <- function (Jday, Tmax_C, Tmin_C, lat_radians, 
                          AvgT = (Tmax_C + Tmin_C)/2, albedo = 0.18, TerrestEmiss = 0.97, 
                          aspect = 0, slope = 0, forest = 0, PTconstant=1.26, AEparams=list(vp=NULL, opt="linear"))
{
  cloudiness <- EstCloudiness(Tmax_C, Tmin_C)
  DailyRad <- NetRad(lat_radians, Jday, Tmax_C, Tmin_C, albedo, forest, slope, aspect, AvgT, cloudiness, TerrestEmiss, AvgT, AEparams=AEparams)
  potentialET <- PTpet(DailyRad, AvgT, PTconstant)
  potentialET[which(potentialET < 0)] <- 0
  potentialET[which(Tmax_C == -999 | Tmin_C == -999)] <- (-999)
  return(potentialET)
}


PET_Looped <- function(X, Jday = J) {
  out <- matrix(nrow = nrow(X), ncol = length(Jday), NA)
  for (i in 1:nrow(X)) {
    out[i, ] <- PET_fromTemp(Jday=Jday, Tmax_C=X$Tx[i], 
                             Tmin_C=(X$Tx[i]-X$Trange[i]),
                             lat_radians=X$lat_radians[i]
                    
      
                      # lat=X$lat[i],
                      # Jday=Jday, Tx=X$Tx[i], 
                      # Tn=(X$Tx[i]-X$Trange[i]), 
                      # X$slope[i],X$aspect[i],units="Wm2"
                      )
  }
  out <- as.data.frame(out)
  names(out) <- paste("Jday", Jday, sep = "")
  return(out)
}

J <- seq(from = 1, to = 365, by = 5)

n <- 10
set.seed(1234)
X <- data.frame(Tx= runif(n, min = 5, max = 30), Trange = runif(n, min = 2,
                   max = 16), lat=runif(n, min = 0, max = 1.1))  # 1.1 radians lat is where?
X$lat_radians = X$lat*pi/180

Y <- PET_Looped(X,Jday = J)
#
# You can ignore all the warnings, remember Errors=bad, warnings=not so much 
# So lets move on and build our summary graph
par(cex.axis = 0.7, cex.lab = 0.8)
plot(J, Y[1, ], type = "l", xlab = "Day of Year", 
     ylab = "Potential Evapotranspiration (meters)")
for (i in 2:n) {
  lines(J, Y[i, ], type = "l", col = i)
}

lat=c(0.1,.77,1.1)
lat_radians.arg=lat*pi/180

PET_Looped.seq <- multisensi(model=PET_Looped, reduction=NULL, center=FALSE,
                                design.args = list( Tx = c(5,15,25), 
                                                    Trange = c(2,9,16),
                                                    lat_radians=lat_radians.arg))

print(PET_Looped.seq, digits = 2)
#
# 3.2 Graphical representation of sensitivity indices
#
# dev.off() # Clean up previous par()
plot(PET_Looped.seq, normalized = TRUE, color = terrain.colors, gsi.plot = FALSE)#normalized the upper subplot shows the extreme (tirets), #inter-quartile (grey) and median (bold line) output values
title(xlab = "Days of the Year.")
plot(PET_Looped.seq, normalized = FALSE, color = terrain.colors, gsi.plot = FALSE)
title(xlab = "Days of the Year.")



# 4 Multivariate sensitivity analysis based on PCA
# read the vignette, though note we are using the multisensi() function
# to run our model (i.e. no “design” variable, and model=Solar_Looped)
PET_Looped.pca <- multisensi(model=PET_Looped, reduction=basis.ACP, scale=FALSE,
                                design.args = list( Tx = c(5,15,25), 
                                                    Trange = c(2,9,16),
                                                    lat_radians=lat_radians.arg))

summary(PET_Looped.pca, digits = 2)
# 4.2 Graphical representation for PCA based analysis with 
# explanation in vignette. These graphs require the plot window to be larger
# and might give "Error in plot.new() : figure margins too large". 
# If so expand the plot window.
# dev.off()
plot(PET_Looped.pca, graph = 1)
plot(PET_Looped.pca, graph = 2)
title("Most Sensitive Variables, PET")
plot(PET_Looped.pca, graph = 3)
#
# 5.1 Polynomial reduction of the multivariate output
# Skip 5.1 Polynomial reduction for now and move on to
# 6 Alternative methods of sensitivity analysis
# 6.1 With Sobol2007 implemented in the package sensitivity
# 
m <- 10000
lat.rand=runif(m, min = 0.0, max = 1.1)
lat.rand.radians=lat.rand*pi/180
Xb <- data.frame(Tx = runif(m, min = 5, max = 30), 
                 Trange = runif(m, min = 2,max = 16), 
                 lat_radians=lat.rand.radians)

PET_Looped.seq.sobol <- multisensi(design = sobol2007, model = PET_Looped,
                                      reduction = NULL, analysis = analysis.sensitivity, center = TRUE,
                                      design.args = list(X1 = Xb[1:(m/2), ], X2 = Xb[(1 + m/2):m, ], nboot = 100),
                                      analysis.args = list(keep.outputs = FALSE))
#
# Note, this is a good time time get a drink of water and/or pee as 
# it is running the function m=10,000 times (a few minutes).
#
print(PET_Looped.seq.sobol, digits = 2)
# dev.off()
plot(PET_Looped.seq.sobol, normalized = TRUE, color = terrain.colors)



## GRAD STUDENT QUESTION - SOILSTORAGE ----

#example
SoilStorage(S_avg=120, field_capacity=0.2, soil_water_content=0.1, porosity=0.3)

# S=(1000/CN)-10



# Solar(lat, Jday, Tx, Tn, albedo=0.2, forest=0, slope=0, aspect = 0,
#      units="kJm2d")
# Note that the EcoHydRology::Solar() function is for specific days, 
# as such, we will want to create a function to loop through our period
# of interest:

CN <- seq(from = 1, to = 100, by = 5)
S = sapply(CN, makeS)

SoilStorage_Looped <- function(X, S_avg=S) {
  out <- matrix(nrow = nrow(X), ncol = length(S_avg), NA)
  for (i in 1:nrow(X)) {
    out[i, ] <- SoilStorage(
                      S_avg = S_avg,
                      field_capacity = X$fc[i],
                      soil_water_content = X$swc[i],
                      porosity = X$por[i]
      
                      # lat=X$lat[i],
                      # Jday=Jday, Tx=X$Tx[i], 
                      # Tn=(X$Tx[i]-X$Trange[i]), 
                      # X$slope[i],X$aspect[i],units="Wm2"
      )
  }
  out <- as.data.frame(out)
  names(out) <- paste("S_avg", S_avg, sep = "")
  return(out)
}

makeS <- function(CN) {
  S=1000/CN-10
  return(S)
}



n <- 10
set.seed(1234)

#need to adjust these values
X <- data.frame(fc = runif(n, min = 0.01, max = .35), swc = runif(n, min = 0.01,
                      max = 1), por = runif(n, min = 0.01, max = .25)
                )  

#View(X)
#
Y <- SoilStorage_Looped(X, S_avg = S)
#
# You can ignore all the warnings, remember Errors=bad, warnings=not so much 
# So lets move on and build our summary graph
par(cex.axis = 0.7, cex.lab = 0.8)
plot(CN, Y[1, ], type = "l", xlab = "Average S", 
     ylab = "Initial Abstraction (mm)")
title("Initial Abstraction Depth v. NRCS Curve Number")
for (i in 2:n) {
  lines(CN, Y[i, ], type = "l", col = i)
}
# 
# Well, that is kewl, yet expected
#
# Multisensitivities
# 3 Sequential univariate sensitivity analyses
# 3.1 Calculation of sensitivity indices
# X <- expand.grid(fc=c(0.01,0.15,0.25), swc=c(0.01,0.2,0.3), swc=c(0.01,0.2,0.3))
# Y <- SoilStorage_Looped(X, S_vec)
SoilStorage_Looped.seq <- multisensi(design=expand.grid, model=SoilStorage_Looped, reduction=NULL, center=FALSE,
                               design.args = list( S_avg=S,
                                                   fc=c(0.01,0.15,0.25),
                                 swc=c(0.01,0.2,0.3),
                                 por=c(0.01,0.2,0.25)))


print(SoilStorage_Looped.seq, digits = 2)
#
# 3.2 Graphical representation of sensitivity indices
#
# dev.off() # Clean up previous par()
plot(SoilStorage_Looped.seq, normalized = TRUE, color = terrain.colors, gsi.plot = FALSE)#normalized the upper subplot shows the extreme (tirets), #inter-quartile (grey) and median (bold line) output values
title(xlab = "S Average")
plot(SoilStorage_Looped.seq, normalized = FALSE, color = terrain.colors, gsi.plot = FALSE)
title(xlab = "S Average")




# 4 Multivariate sensitivity analysis based on PCA
# read the vignette, though note we are using the multisensi() function
# to run our model (i.e. no “design” variable, and model=Solar_Looped)
SoilStorage_Looped.pca <- multisensi(model=SoilStorage_Looped, reduction=basis.ACP, scale=FALSE,
                               design.args = list( 
                                 fc=c(0.1,0.2,0.25),
                                 swc=c(0.3,0.6,1),
                                 por=c(0.1,0.2,0.25)
                                 ))

summary(SoilStorage_Looped.pca, digits = 2)
# 4.2 Graphical representation for PCA based analysis with 
# explanation in vignette. These graphs require the plot window to be larger
# and might give "Error in plot.new() : figure margins too large". 
# If so expand the plot window.
# dev.off()
plot(SoilStorage_Looped.pca, graph = 1)
plot(SoilStorage_Looped.pca, graph = 2)
plot(SoilStorage_Looped.pca, graph = 3)
#
# 5.1 Polynomial reduction of the multivariate output
# Skip 5.1 Polynomial reduction for now and move on to
# 6 Alternative methods of sensitivity analysis
# 6.1 With Sobol2007 implemented in the package sensitivity
# 
library(sensitivity)
m <- 10000
Xb <- data.frame(
  fc = runif(n, min = 0.01, max = .35), swc = runif(n, min = 0.01,
  max = 1), por = runif(n, min = 0.01, max = .25)
)

SoilStorage_Looped.seq.sobol <- multisensi(design = sobol2007, model = SoilStorage_Looped,
                                     reduction = NULL, analysis = analysis.sensitivity, center = TRUE,
                                     design.args = list(X1 = Xb[1:(m/2), ], X2 = Xb[(1 + m/2):m, ], nboot = 100),
                                     analysis.args = list(keep.outputs = FALSE))
#
# Note, this is a good time time get a drink of water and/or pee as 
# it is running the function m=10,000 times (a few minutes).
#
print(SoilStorage_Looped.seq.sobol, digits = 2)
# dev.off()
plot(SoilStorage_Looped.seq.sobol, normalized = TRUE, color = terrain.colors)

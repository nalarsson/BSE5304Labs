PLossModel = function (RunoffModel=TIC05,tau=9.3, dt=1, kF=.015, TI=5, TIC_rast=TIC) {

  RunoffModel=TMWB
  RunoffModel = CNmodel(CNmodeldf = RunoffModel, CNavg=VSAParams$CN[TI], #MUSLE$TIclass[i]],
                  declat=myflowgage$declat,declon=myflowgage$declon)
  RunoffModel$qpeak=RunoffModel$Qpred/3600/24/1000*myflowgage$area/nTIclass*10^6    #MUSLE$areaSQKM*10^6 #m^3/sec
  RunoffModel$sed=(RunoffModel$Qpred*RunoffModel$qpeak*myflowgage$area/nTIclass)^.56*MUSLE$KCPLSCFRG118[TI]  #MUSLE$areaSQKM*100)^.56*MUSLE$KCPLSCFRG118[1] 
  

  tmpdf=data.frame(date=RunoffModel$date,
                    Rt=RunoffModel$Qpred/1000, 
                    Tavg=(RunoffModel$MaxTemp+RunoffModel$MinTemp)/2)
  tmpdf$MF=0
  tmpdf$DF=0

  attach(tmpdf)
  
#
# Loop to solve MF and DF
  for (i in 2:length(date)){
    if(MF[i]<=MF[i-1]){
    MF[i]=MF[i-1]*exp(-dt/tau)-DF[i-1]
    }
    DF[i]=MF[i]*(kF*MF[i]*Rt[i]/(1+kF*MF[i]*Rt[i]))
  }
  tmpdf$MF=MF
  tmpdf$DF=DF
  detach(tmpdf)
  rm(list=c("MF","DF")) # Clean up the environment 
  #dev.off() #reset graphics device
  plot(tmpdf$date,tmpdf$MF)
  plot(tmpdf$date,tmpdf$DF)

  TIclass = c(5,4,3,2,1)
  VSAlookup = c(1:5)
  
  index=which(TIclass==TI)
  muTS_tmp=(((520-0)/5)*VSAsol$TIClass[index]+0)
  MS_tmp=(((18.5-3.7)/5)*VSAsol$TIClass[index]+3.7)*2000 
  
  QS= 3.0 # A guess using the middle of the range 1-5
  TR=20   # reference Temperature from Table 2.
  tmpdf$muS= muTS_tmp*QS^((tmpdf$Tavg-TR)/10)  # Eq. 5
  tmpdf$DS=(tmpdf$muS*MS_tmp*tmpdf$Rt)/10^6          # Eq. 4
  plot(tmpdf$date,tmpdf$DS)
  
  
  DPLT=data.frame()
  DPLTtmp=data.frame(date=RunoffModel$date,
                  Rt=RunoffModel$Qpred/1000.0,
                  Tavg=(RunoffModel$MaxTemp+RunoffModel$MinTemp)/2)
  DPLTtmp$B=min(RunoffModel$Qmm)*myflowgage$area*1000*1000/1000 # m^3/day
  muTB=2.1*10^(-5) # Easton Table 2
  QB=2.2           # Easton Table 2
  TB=17            # Easton Table 2
  DPLTtmp$muB=muTB*QB^((DPLTtmp$Tavg-TB)/10)  # Easton eq. 10
  DPLTtmp$LB=DPLTtmp$muB*DPLTtmp$B     # Easton eq. 9
  plot(DPLTtmp$date,DPLTtmp$LB)
  
  rclconst=matrix(ncol=2)
  rclconst=rbind(rclconst,c(TI,mean(DPLTtmp$LB)))
  
  
  setwd("~/2023/BSE5304Labs11/pdfs/Lab11")
  pdfname=paste0("meanPlossTI",TI)
  png(pdfname)
  TIC_loss=reclassify(TIC_rast,rclconst)
  dev.off()
  
  DPTIvalue=paste0("DPTI",TI)
  assign(DPTIvalue, tmpdf)
  #assign(DPLT, DPLTtmp)
  
  return(TIC_loss)
  #return(DPTIvalue)
}

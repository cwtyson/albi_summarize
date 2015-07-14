#### THIS SCRIPT IS FOR DOING PERMUTATIONS TO FIND THE BEST BUFFER SIZE ####

library(foreach)
library(doParallel)
library(doSNOW)
library(sp)
library(rgdal)
library(raster)
library(geosphere)
library(adehabitatLT)
library(adehabitat)
library(adehabitatHR)
library(adehabitatHS)
library(adehabitatMA)
library(circular)
library(CircStats)
library(randomForest)
library(dplyr)


source("C:/Users/Grant/Dropbox/GrantHumphriesBackup/Projects/Albatross/Codes/albi_summarize/streamlined_functions.R")

###################################################################################
#### MODEL SETTINGS ####
###################################################################################

TmBuff<-1300
fptRad<-5000    
resTRad<-500
resTTIME<-50
NPoints <-9
SubSamp<-120
Theta<-0.7
StepSize<-2
TurnSens<-60
species<-"wAAL"          #BBAL
Timezone<-1              # 1 is for BBAL
outWS<-"C:/Users/Grant/Dropbox/GrantHumphriesBackup/Projects/Albatross/WAAL/analysis/"

Res<-10
#DistBuff<-5000
#EndDist<-5000


GETdir<-paste(outWS, "/resampled_", Res, "_sec/",sep="")

DistBuffs<-c(5000,8000,10000,12000,15000,20000,25000,30000,40000,50000)



filenames=list.files(path=GETdir, pattern="*.txt",full.names=TRUE)
datalist = lapply(filenames, function(x){read.table(file=x,sep=",",header=T)})
filenames=list.files(path=GETdir, pattern="*.txt")

cl<-makeCluster(7, outfile="")
registerDoSNOW(cl)


for(XX in 1:100){
  for(DistBuff in DistBuffs){
  
    
    EndDist<-DistBuff - (DistBuff/2)
    
    
    MATRIX<-matrix(ncol=50,nrow=0)
    
    count<-1
    for(k in datalist){
      
      Data<-filenames[count]
      
      print(Data)
      events<-which(k$X.event==1)
      
      Events<-foreach(x=events,.combine='rbind',.errorhandling='remove',.verbose=TRUE,
                      .packages=c('dplyr','raster','rgdal','sp','adehabitat','adehabitatLT','geosphere')) %dopar% {
  
        OBJ<-Summarize.at.point(x,Data,TmBuff,DistBuff,EndDist,fptRad,resTRad,resTTIME,
                                NPoints,SubSamp,Theta,StepSize,TurnSens,species,Timezone,outWS)
        return(OBJ)
      }
      
      Events<-cbind(Events,1)
      colnames(Events)[length(colnames(Events))]<-"Event"
    
      nevn.i<-which(k$X.event==0)
      nevets<-sample(nevn.i,size=length(events),replace=FALSE)
      
      NEvents<-foreach(x=nevets,.combine='rbind',.errorhandling='remove',.verbose=TRUE,
                       .packages=c('dplyr','raster','rgdal','sp','adehabitat','adehabitatLT','geosphere')) %dopar% { 
        OBJ2<-Summarize.at.point(x,Data,TmBuff,DistBuff,EndDist,fptRad,resTRad,resTTIME,
                                 NPoints,SubSamp,Theta,StepSize,TurnSens,species,Timezone,outWS)
        return(OBJ2)
      }
      
      NEvents<-cbind(NEvents,0)
      colnames(NEvents)[length(colnames(NEvents))]<-"Event"
      
      try({
        OUTPUT<-rbind(Events,NEvents)
        MATRIX<-rbind(MATRIX,OUTPUT)
      })
      
      
      count<-count+1
    }
    
    
    
    DF<-data.frame(MATRIX)
    
    DF$T_before<-as.double(as.character(DF$T_before))
    DF$T_after<-as.double(as.character(DF$T_after))
    DF$Sp_before<-as.double(as.character(DF$Sp_before))
    DF$Sp_after<-as.double(as.character(DF$Sp_after))
    DF$Bath<-as.double(as.character(DF$Bath))
    DF$sst<-as.double(as.character(DF$sst))
    DF$slp_before<-as.double(as.character(DF$slp_before))
    DF$slp_after<-as.double(as.character(DF$slp_after))
    DF$ul_before<-as.double(as.character(DF$ul_before))
    DF$ul_after<-as.double(as.character(DF$ul_after))
    DF$mwp_before<-as.double(as.character(DF$mwp_before))
    DF$mwp_after<-as.double(as.character(DF$mwp_after))
    DF$chn_before<-as.double(as.character(DF$chn_before))
    DF$chn_after<-as.double(as.character(DF$chn_after))
    DF$FPT_before<-as.double(as.character(DF$FPT_before))
    DF$FPT_after<-as.double(as.character(DF$FPT_after))
    DF$ResT_before<-as.double(as.character(DF$ResT_before))
    DF$ResT_after<-as.double(as.character(DF$ResT_after))
    DF$ResT_during<-as.double(as.character(DF$ResT_during))
    DF$NumTurnsBefore<-as.double(as.character(DF$NumTurnsBefore))
    DF$NumTurnsAfter<-as.double(as.character(DF$NumTurnsAfter))
    DF$SPslope<-as.double(as.character(DF$SPslope))
    DF$Fnspd<-as.double(as.character(DF$Fnspd))
    DF$Sp_around_event<-as.double(as.character(DF$Sp_around_event))
    DF$Sp_at_event<-as.double(as.character(DF$Sp_at_event))
    DF$torSeg1<-as.double(as.character(DF$torSeg1))
    DF$AccelSeg1<-as.double(as.character(DF$AccelSeg1))
    DF$torSeg2<-as.double(as.character(DF$torSeg2))
    DF$AccelSeg2<-as.double(as.character(DF$AccelSeg2))
    
    DF<-DF[,4:50]
    DF<-DF[complete.cases(DF),]
    
    
    rf1 <- randomForest(Event~.,ntree=10000,mtry=15,data = DF, importance=TRUE)
    
    Conf.Mat<-rf1$confusion
    PCCzero<-round(1-Conf.Mat[5],2)
    PCCone<-round(1-Conf.Mat[6],2)
    
    
    print(paste(DistBuff,PCCzero,PCCone,sep=","))
    
    ROW<-c(DistBuff,PCCzero,PCCone,XX)
    write(ROW,file = "C:/Temp/LOG4.txt",ncolumns=4,append=TRUE,sep = ",")
  }

  
}

stopCluster(cl)

Data.Tab<-tbl_df(read.table("C:/Temp/LOG4.txt",sep=",",header=F))

names(Data.Tab)<-c("resolution","pcczero","pccone","run")

GR<-Data.Tab %>% group_by(resolution) %>% summarise(pcczeroMEAN=mean(pcczero),pcconeMEAN=mean(pccone),pcczeroSD=sd(pcczero),pcconeSD=sd(pccone))
GR$total<-(GR$pcczeroMEAN + GR$pcconeMEAN)/2


rowMeans()













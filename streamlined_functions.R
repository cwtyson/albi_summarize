library(dplyr)
library(sp)
library(rgdal)
library(geosphere)
library(CircStats)
library(circular)
library(adehabitat)
library(adehabitatLT)
library(adehabitatHR)
library(adehabitatHS)
library(adehabitatMA)
library(foreach)
library(doSNOW)
library(doParallel)
library(raster)
library(zoo)


resample <- function(Data, setTZ, Res, WD, outWD,species,fptRad,resTRad,resTTIME){
  
  VARLIST<-c("X.Latitude","X.Longitude","X.Speed","Uplift","SWH","Charn","Bathy","Dist_to_coast","X10MWind",
             "MWD","SST","SLP","MWP","MetWind","WindDifferential","WindDirToFlight","Period")
  
  X<-tbl_df(read.table(paste(WD, "/", Data, sep=""), sep=",",header=T))
  #X$WindDirToFlight<-X$WindDirs
  
  #   GMTtime<-as.POSIXct(strptime(paste(X$X.Date,X$X.Time),format="%d-%m-%Y %H:%M:%S",tz="GMT"),format ="%Y-%m-%d %H:%M:%S",tz="GMT")
  #   attributes(GMTtime)$tzone<-"Indian/Maldives"
  #   
  #   X$LocalDateTime<-GMTtime
  #   
  Tm<-as.POSIXct(strptime(X$LocalDateTime, format="%Y-%m-%d %H:%M:%S",tz=setTZ),format ="%Y-%m-%d %H:%M:%S",tz=setTZ)
  
  Tm<-as.POSIXct(Tm,format="%Y-%m-%d %H:%M:%S",tz=setTZ)
  
  Z<-project(cbind(X$X.Longitude,X$X.Latitude), "+proj=utm +zone=42 ellps=WGS84")    
  M<-as.ltraj(Z,Tm,X$X.AnimalID)
  
  refda<-strptime("00:00", "%H:%M")
  M2<-setNA(M,refda,10,tol=10,units="sec")
  M2<-sett0(M2,refda,10,tol=10,units="sec")
  
  ## Resample to desired resolution
  M2<-adehabitat::subsample(M2,Res)
  
  X2<-tbl_df(data.frame(M2))
  X2time<-as.POSIXct(X2$date,"%Y-%m-%d %H:%M:%S",tz=setTZ)
  X2$date<-format(X2time,"%Y-%m-%d %H:%M:%S",tz=setTZ)
  
  X.T<-strptime(X$LocalDateTime,"%Y-%m-%d %H:%M:%S",tz=setTZ)
  
  X3<-matrix(nrow=0,ncol=19,byrow=TRUE)
  colnames(X3)<-c("LocalDateTime","X.event",VARLIST)
  
  ##### Creates a new matrix, merging the old data to the new data (within 10 seconds of a resampled point)
  print("-----------------------------------------------------------------------")
  print("-----------------------------------------------------------------------")
  print(Data)
  print("-----------------------------------------------------------------------")
  print("Merging data into a matrix... computing turn angles")
  print("-----------------------------------------------------------------------")
  
  
  for(tt in 1:nrow(X2)){
    
    ### Will run every 100 rows as a progress meter
    if(tt %% 100 == 0){print(paste(tt," of ",nrow(X2)," rows",sep=""))}
    
    X2.T<-strptime(X2$date[tt],"%Y-%m-%d %H:%M:%S",tz=setTZ)
        
    ## Get points within the average resolution of points in the track

    alldiff<-abs(difftime(X2.T,X.T))
    Index<-which(alldiff==min(alldiff))
    Index<-Index[1]

    ROW<-data.frame(LocalDateTime=X2.T,X.event=0,X[Index,VARLIST])
    
    X3<-rbind(X3,ROW)
  }
  
  print("-----------------------------------------------------------------------")
  print("Calculating the bearing between successive points")
  print("-----------------------------------------------------------------------")
  
  
  Trck<-tbl_df(data.frame(X3)) %>% filter(X.Longitude!="NA")
  
  ####### Gets the bearing between successive points and calculates the turn angles    
  for(i in 1:nrow(Trck)){
    if(i == nrow(Trck)){
      print("end of file")
      break
    }
    B<-tryCatch({bearing(cbind(Trck$X.Longitude[i],Trck$X.Latitude[i]),cbind(Trck$X.Longitude[i+1],Trck$X.Latitude[i+1]))},
                error = function(e){return(NA)})
    
    D<-tryCatch({distCosine(cbind(Trck$X.Longitude[i],Trck$X.Latitude[i]),cbind(Trck$X.Longitude[i+1],Trck$X.Latitude[i+1]))},
                error = function(e){return(NA)})        
    Trck$X.bearing[i]<-B
    Trck$X.dist[i]<-D
  }
  
  Tm3<-as.POSIXct(strptime(Trck$LocalDateTime,format="%Y-%m-%d %H:%M:%S",tz=setTZ),format ="%Y-%m-%d %H:%M:%S",tz=setTZ)
  
  Z3<-project(cbind(Trck$X.Longitude,Trck$X.Latitude), "+proj=utm +zone=42 ellps=WGS84")    
  
  Trck$ID<-species
  M3<-as.ltraj(Z3,Tm3,id=Trck$ID)
  
  print("-----------------------------------------------------------------------")
  print("Now calculating first passage time and residence time")
  print("-----------------------------------------------------------------------")
  
  # Calculates first passage time for the track
  N<-fpt(M3,fptRad)
  
  # Calcluates residence time for the track
  A<-residenceTime(M3,radius=resTRad,maxt=resTTIME,units="seconds")
  
  ## Appends FPT and RESIDENCE TIME to the original data frame
  # This command is in here because FPT and RES can be sensitive to changes
  
  P<-c(Trck,tbl_df(N[[1]]))
  Q<-c(P,tbl_df(A[[1]][2]))
  
  Trck<-tbl_df(as.data.frame(Q))
  
  names(Trck)[length(names(Trck))]<-"resT"
  names(Trck)[length(names(Trck))-1]<-"fpt"
  
  ### We add a spline here to help remove NAs where they occu
    
  Trck[,"resT"]<-abs(na.spline(Trck[,"resT"]))
  Trck[,"fpt"]<-abs(na.spline(Trck[,"fpt"]))
  
  ### Change the last 10 points to remove the error (they are over-inflated)
  Trck[(nrow(Trck)-10):nrow(Trck),"resT"]<-mean(Trck$resT,na.rm=TRUE)
  Trck[(nrow(Trck)-10):nrow(Trck),"fpt"]<-mean(Trck$fpt,na.rm=TRUE)
  
  M3<-data.frame(M3)
  Trck$TurnAngles=round(M3$rel.angle*180/pi, 2)
  
  
  
  #### This will flag the track for points where the bird is either sitting or flying
  print("-----------------------------------------------------------------------")
  print("Now calculating when birds are sitting or flying")
  print("-----------------------------------------------------------------------")
  
  
  Tspeed<-dplyr::select(Trck,X.Speed)
  Tspeed<-as.vector(unlist(Tspeed[,1]))
  
  SampSps<-sapply(seq(1,nrow(Trck),by=1),function(x){mean(Tspeed[x:x-10],na.rm=T)})
  SampSps[is.nan(SampSps)]<-0
  Trck$Sit.Fly<-sapply(SampSps,function(x){if(x<6){x="Sit"}else{x="Fly"}})
  
  ################################################################
  ### Now we have to get the events associated with the points ###
  ################################################################
  print("Now associating events with tracks")
  print("----------------------------------------------------------------------")
  
  EvIndices<-which(X$X.event==1)
  X2.Tm<-strptime(Trck$LocalDateTime,"%Y-%m-%d %H:%M:%S",tz=setTZ)
  
  for(j in EvIndices){
    #print(j)
    EvTime<-X.T[j]
    alldiff<-abs(difftime(EvTime,X2.Tm))
    Event<-which(alldiff==min(alldiff))
    Event<-Event[1]
    Trck$X.event[Event]<-1
  }
  
    
  
  ################################################################################
  ###### Write the output table ######
  ################################################################################
  
  ## Create folder (if it doesn't exist) for the current resampled resolution
  if (!file.exists(paste(outWD, "/resampled_", Res, "_sec/", sep = "")))
  {
    dir.create(paste(outWD, "/resampled_", Res, "_sec/", sep = "")) # Create directory if it does not exist   
  }
  
  ## Name of output file
  Outname<-paste(outWD, "/resampled_", Res, "_sec/", substr(Data,1,nchar(Data)-4),"_Output.txt", sep="")      
  
  print(paste("Table for ",substr(Data,1,(nchar(Data)-10))," now being written...",sep=""))
  print("-------")
  write.table(Trck,Outname,sep=",",row.names=FALSE)
  print("complete...")
  
}

WindDifferential<-function(a,b){
  
  test<-is.na(a)
  if(test == TRUE){
    Flight<-"No Flight"
  }
  if(test == FALSE){
    x<-(a-b + 180) %% 360 - 180
    if(abs(x) < 45){Flight<-"Headwind"}
    if(abs(x) > 45 && abs(x) < 135){Flight<-"Crosswind"}
    if(abs(x) > 135){Flight<-"Tailwind"}
  }  
  return(Flight)
}



WindDiffDirectionality<-function(a,b){
  
  test<-is.na(a)
  if(test == TRUE){
    WindDirToFlight<-"No Flight"
  }
  if(test == FALSE){ 
    x<-(a-b + 180) %% 360 - 180
    if(x < 0){WindDirToFlight<-"Right"}
    if(x > 0){WindDirToFlight<-"Left"}
  }  
  return(WindDirToFlight)
}



circMean<-function(input){
  if(input<0){
    circMean<- 360+input
  } else {
    circMean<- input
  }
  return(circMean)
}


Summarize.at.point<-function(Index,Data,TmBuff,DistBuff,EndDist,NPoints,SubSamp,Theta,StepSize,TurnSens,species,timezone,outWD){
    
  ## Set timezone
  setTZ <- ifelse(timezone == 1,"Indian/Maldives",ifelse(timezone == 2,"Indian/Mahe", stop("No time zone set") ))
  
  if(Index<NPoints) stop("Index too close to beginning of track...")
  
  Eindex<-Index
  
  SumMat<-matrix(ncol=41+(StepSize*5))
  
  BirdName<-paste("Bird_",substr(Data,gregexpr(pattern="_",Data)[[1]][1]+1,gregexpr(pattern="_",Data)[[1]][2]-1),sep="")
  
  
  FREAD<-paste(outWD, "/resampled_", Res, "_sec/",Data, sep="")
  
  
  X<-tbl_df(read.table(FREAD,sep=",",header=T))
  
  
  ########################################################################
  ###### This section sets up the 10 sec tracks for data extraction ######
  ########################################################################
  
  a1<-cbind(X$X.Longitude[Eindex],X$X.Latitude[Eindex])
  b1<-cbind(X$X.Longitude,X$X.Latitude)      
  DistanceToEvent<-pointDistance(a1,b1,longlat=TRUE)
  
  Tim<-as.POSIXct(strptime(X$LocalDateTime,format="%Y-%m-%d %H:%M:%S",tz=setTZ),format ="%Y-%m-%d %H:%M:%S",tz=setTZ)  
  TmDiff<-difftime(Tim,Tim[Eindex])
  
  EF<-cbind(TmDiff,DistanceToEvent)
  
  ### Selects points that are DistBuff Meters away from event and between -TmBuff and 0
  BeforeIndices<-which(EF[,2] < DistBuff)
  BeforeIndices<-BeforeIndices[which(EF[BeforeIndices,1] > -TmBuff & EF[BeforeIndices,1] < 0)]
  BeforeEvent<-X[BeforeIndices,]
  
  
  ### Selects points that are DistBuff Meters away from event and between TmBuff and 0
  AfterIndices<-which(EF[,2] < EndDist)
  AfterIndices<-AfterIndices[which(EF[AfterIndices,1] < TmBuff & EF[AfterIndices,1] > 0)]
  AfterEvent<-X[AfterIndices,]
  
  
  
  if(length(BeforeIndices)>2 & length(AfterIndices)> 2){
    
    ####################################################################################################
    ######### This code summarizes the track ############
    ####################################################################################################
    
    
    #Tortuosity before and after (Non-subbed data)
    T_before<-(distCosine(c(BeforeEvent$X.Longitude[1],BeforeEvent$X.Latitude[1]),c(X$X.Longitude[Eindex],X$X.Latitude[Eindex])))/(sum(BeforeEvent$X.dist,na.rm=TRUE))
    T_after<-(distCosine(c(X$X.Longitude[Eindex],X$X.Latitude[Eindex]),c(AfterEvent$X.Longitude[nrow(AfterEvent)],AfterEvent$X.Latitude[nrow(AfterEvent)])))/(sum(AfterEvent$X.dist,na.rm=TRUE))
    
    
    #Straightness by T_before / T_after
    
    if(T_before > Theta){
      Approach<-"Straight"
    }else{Approach<-"Tortuous"}
    
    if(T_after > Theta){
      Depart<-"Straight"
    }else{Depart<-"Tortuous"}
    
    #average speed around event  
    Sp_around_event<-mean(X$X.Speed[Eindex-3:Eindex+3],na.rm=TRUE)
    Sp_at_event<-X$X.Speed[Eindex]
    
    #Mean speed before and after
    Sp_before<-mean(BeforeEvent$X.Speed,na.rm=TRUE)
    Sp_after<-mean(AfterEvent$X.Speed,na.rm=TRUE)
    
    #Bathymetry at feeding location 
    Bath<-X$Bathy[Eindex]
    
    #Mean wind speed
    #wsp_before<-mean(BeforeEvent$X10MWind,na.rm=TRUE)
    #wsp_after<-mean(AfterEvent$X10MWind,na.rm=TRUE)
    
    #Mean SST at feeding location
    sst<-X$SST[Eindex]
    
    #Mean Sea level pressure  SLP
    slp_before<-mean(BeforeEvent$SLP,na.rm=TRUE)
    slp_after<-mean(AfterEvent$SLP,na.rm=TRUE)
    
    
    #Mean Uplift  Uplift
    ul_before<-mean(BeforeEvent$Uplift,na.rm=TRUE)
    ul_after<-mean(AfterEvent$Uplift,na.rm=TRUE)
    
    
    #Mean Wave Period  MWP
    mwp_before<-mean(BeforeEvent$MWP,na.rm=TRUE)
    mwp_after<-mean(AfterEvent$MWP,na.rm=TRUE)
    
    #Mean charnock  Charn
    chn_before<-mean(BeforeEvent$Charn,na.rm=TRUE)
    chn_after<-mean(AfterEvent$Charn,na.rm=TRUE)
    
    
    #### Most common WindDifferential category
    WD_before<-names(which.max(table(BeforeEvent$WindDifferential)))
    WD_after<-names(which.max(table(AfterEvent$WindDifferential)))
    
    #### Most common WindDireToFlight category
    WDTF_before<-names(which.max(table(BeforeEvent$WindDirToFlight)))
    WDTF_after<-names(which.max(table(AfterEvent$WindDirToFlight)))
    
    
    #### Defines the wind differential from the index to NPoints distance
    
    PHeadW<-WindDifferential(bearing(cbind(X$X.Longitude[Eindex-NPoints],X$X.Latitude[Eindex-NPoints]),
                                     cbind(X$X.Longitude[Eindex],X$X.Latitude[Eindex])),
                             circMean(deg(circ.mean(rad(X$MetWind[Eindex-NPoints:Eindex])))))
    
    if(PHeadW=="Headwind"){
      HeadWind<-"Yes"
    }else{HeadWind<-"No"}
    
    #####################################################
    ########## THIS IS THE OLD WAY OF DEFINING HEAD WIND
    ########## UPDATED MAY 6, 2015 ###########
    #####################################
    
    
    #TABLE<-table(X$WindDifferential[Eindex-NPoints:Eindex])
    #headn<-which(names(TABLE)=="Headwind")
    #PecHeadw<-as.double(TABLE[[headn]])/sum(TABLE)
    
    #     # If bird is in headwind > 40% of the time, call it a headwind
    #     if(PecHeadw > 0.2){
    #       HeadWind<-"Yes"
    #     }else{HeadWind<-"No"}
    ##################################################################    
    
    #### Get mean FPT prior to and after events
    FPT_before<-mean(BeforeEvent$fpt,na.rm=TRUE)
    FPT_after<-mean(AfterEvent$fpt,na.rm=TRUE)
    
    #### Get mean Residence time prior to and after events
    ResT_before<-mean(BeforeEvent$resT,na.rm=TRUE)
    ResT_after<-mean(AfterEvent$resT,na.rm=TRUE)
    
    #### Get Residence time OF event
    ResT_during<-X$resT[Eindex]
    
    #### Extract Times from events ####
    Date<-as.character(X$LocalDateTime[Eindex])
    
    #### Extract Day / Night for the particular event ####
    Period<-as.character(X$Period[Eindex]) 
    
    
    #### Extract Prey mass for the particular event ####
    #PreyMass<-as.character(X$Mass[Eindex])
    
    
    #### Number of turns before event
    NumTurnsBefore<-length(which(BeforeEvent$TurnAngles>TurnSens))    
    
    #### Is there a major turn prior to the event? ####
    MTB<-length(which(BeforeEvent$TurnAngles>90))
    
    if(MTB>0){MajTurnsBefore = "Yes"}else{MajTurnsBefore = "No"}
    
    
    
    ### Number of turns > 60 deg after event
    NumTurnsAfter<-length(which(AfterEvent$TurnAngles>TurnSens))
    
    
    ### Is there a major turn after the event?
    MTA<-length(which(AfterEvent$TurnAngles>90))
    if(MTA>0){
      MajTurnsAfter = "Yes"
    }else{MajTurnsAfter = "No"}
    
    
    
    #### Calculate if a line is "straight" ####
    
    #### We subsample track, calculate tortuosity of that subsample ####
    #### In order to fill in the blanks where "NAs" occur due to subsampling
    #### we have to convert to a track again in order to get the new distances
    #### Also, we have to add the distance from the last point to the event to get total distance
    
    
    ### Is the event a multiple feeding event? ###
    if(sum(BeforeEvent$X.event)+sum(AfterEvent$X.event) > 0){
      
      Multi<-"Yes"
      
    }else{Multi <- "No"}
    
    
    
    #### Pull acceleration from last few points prior to event ####
    
    SPslope<-summary(lm(X$X.Speed[Eindex-5:Eindex]~c(1:length(X$X.Speed[Eindex-5:Eindex]))))$coefficients[2]
    
    Fnspd<-mean(X$X.Speed[Eindex-2:Eindex],na.rm=T)
    
    #### Extract the behaviour type ####
    #Beh_Type<-as.character(X$Type[Eindex])      
    
    #print(paste("Data extracted for index",Eindex,sep=" "))
    #######################################################################################################
    EventData<-cbind(Eindex,
                     BirdName,
                     Date,
                     #Beh_Type,
                     T_before,
                     T_after,
                     Sp_before,
                     Sp_after,
                     Bath,
                     sst,
                     #wsp_before,
                     #wsp_after,
                     slp_before,
                     slp_after,
                     ul_before, 
                     ul_after,
                     mwp_before,
                     mwp_after,
                     chn_before,
                     chn_after,
                     WD_before,
                     WD_after,
                     WDTF_before,
                     WDTF_after,
                     HeadWind,
                     FPT_before,
                     FPT_after,
                     ResT_before,
                     ResT_after,
                     ResT_during,
                     Approach,
                     Depart,
                     NumTurnsBefore,
                     NumTurnsAfter,
                     MajTurnsBefore,
                     MajTurnsAfter,
                     #PreyMass,
                     Period,
                     SPslope,
                     Fnspd,
                     Multi,
                     #App_before,
                     Sp_around_event,
                     Sp_at_event)
    
    
    ### Create a sliding window with a step analysis (e.g. 5 steps = total distance / 5)
    
    start<-Eindex-length(BeforeIndices)
    StepLen<-floor(length(BeforeIndices)/StepSize)
    
    MidLen<-round(StepLen/StepSize)
    
    
    for(st in 1:StepSize){
      
      end<-start + StepLen 
      if(end > Eindex){
        end<-Eindex
      }
      
      midend<-start+MidLen
      
      midstart<-midend+1
      
      
      if(start==midend){
        nstart<-start-1
        seg1bear<-bearing(cbind(X$X.Longitude[nstart],X$X.Latitude[nstart]),
                          cbind(X$X.Longitude[midend],X$X.Latitude[midend]))
      }else{
        seg1bear<-bearing(cbind(X$X.Longitude[start],X$X.Latitude[start]),
                          cbind(X$X.Longitude[midend],X$X.Latitude[midend]))
      }
      
      if(midstart==end){
        nend<-end+1
        seg2bear<-bearing(cbind(X$X.Longitude[midstart],X$X.Latitude[midstart]),
                          cbind(X$X.Longitude[nend],X$X.Latitude[nend]))
      }else{
        seg2bear<-bearing(cbind(X$X.Longitude[midstart],X$X.Latitude[midstart]),
                          cbind(X$X.Longitude[end],X$X.Latitude[end]))
      }
      
      
      
      SegDiff<-(seg1bear-seg2bear + 180) %% 360 - 180
      
      ### Is there a drastic path change along this segment? 
      
      if(SegDiff > TurnSens){
        Pathchange<-"Yes"
        
        
        #         PCbrg<-bearing(cbind(X$X.Longitude[end-10],X$X.Latitude[end-10]),
        #                        cbind(X$X.Longitude[end],X$X.Latitude[end]))
        
        PCbrg<-seg2bear
        
        PCcrmn<-circMean(deg(circ.mean(rad(X$MetWind[end-10:end]))))
        
        ### If there is a path change, does it go into a headwind?
        
        PCWinDif<-WindDifferential(PCbrg,PCcrmn)
        if(PCWinDif=="Headwind"){
          PCheadwind<-"Yes"
        }else{PCheadwind<-"No"}
        
        
      }else{Pathchange<-"No"
            PCheadwind<-"No change"}
      
      
      ### What is the wind differential across the whole segment? 
      brg<-bearing(cbind(X$X.Longitude[start],X$X.Latitude[start]),
                   cbind(X$X.Longitude[end],X$X.Latitude[end]))
      
      crmn<-circMean(deg(circ.mean(rad(X$MetWind[start:end]))))
      
      WinDifSeg<-WindDifferential(brg,crmn)
      
      
      #### Tortuosity and Acceleration across each segment
      
      torSeg<-(distCosine(c(X$X.Longitude[start],X$X.Latitude[start]),
                          c(X$X.Longitude[end],X$X.Latitude[end])))/(sum(X$X.dist[start:end],na.rm=TRUE))
      
      AccelSeg<-summary(lm(X$X.Speed[start:end]~c(1:length(X$X.Speed[start:end]))))$coefficients[2]
      
      ###################################################################################
      
      EventData<-cbind(EventData,torSeg,AccelSeg,Pathchange,PCheadwind,WinDifSeg)
      
      tsegname<-paste("torSeg",as.character(st),sep="")
      asegname<-paste("AccelSeg",as.character(st),sep="")
      pthchngname<-paste("Pathchange",as.character(st),sep="")
      PChwname<-paste("PathChange_headwind",as.character(st),sep="")
      Wdiffname<-paste("WinDifSeg",as.character(st),sep="")
      
      colnames(EventData)[ncol(EventData)-4]<-tsegname
      colnames(EventData)[ncol(EventData)-3]<-asegname
      colnames(EventData)[ncol(EventData)-2]<-pthchngname
      colnames(EventData)[ncol(EventData)-1]<-PChwname
      colnames(EventData)[ncol(EventData)]<-Wdiffname
      
      
      start<-end+1
      #remove(nstart)
      #remove(nend)
    }
    
  }
  
  #print(EventData)
  return(EventData)
  
}


streamlined <- function (WD, outWD, Res, timezone, species, resmp = 1, summarize = 1, filePattern) {
  
  ## Argument descriptions:
  ## WD - workspace (character string, not ending with with a "/")
  ## outWD - output directory (character string, not ending with with a "/")
  ## Resamp - desired resampling interval in seconds (integer)
  ## timezone - 1 = "Indian/Maldives", 2 = "Indian/Mahe" (integer)  ### BBAL = Maldives
  ## species - species code to be added to the data frame (e.g. BBAL, WAAL) (character string)
  ## resmp - resample the data? (1 = yes)
  ## summarize - summarize the data? (1 = yes)
  ## filePattern - character sequence used to pick out files in workspace (e.g. "GPS") (character string)
  ## Hard coded arguments
  TmBuff<-1300
  DistBuff<-25000   #### AT 120 sec resolution, only capturing a few points because birds are flying fast! 
  EndDist<-20000
  fptRad<-800    
  resTRad<-500
  resTTIME<-50
  NPoints <-9
  SubSamp<-120
  Theta<-0.7
  StepSize<-2
  TurnSens<-60
 
  
  ## Set working directory and get files to process
  setwd(WD)
  BirdList<-dir(pattern = filePattern)
  
  ## Set timezone
  setTZ <- ifelse(timezone == 1,"Indian/Maldives",ifelse(timezone == 2,"Indian/Mahe", stop("No time zone set") ))
  
  if(resmp == 1)
  {  
    ## For each data.frame resample at current resolution
    for(Data in BirdList){
      ## For each data.frame, resample
      resample(Data, setTZ, Res, WD, outWD,species,fptRad,resTRad,resTTIME)
    }
  }
  
  ## If summarize option has been selected, then summarize each resampled data frame
  if(summarize == 1)
  {
    
    ## Get output directory for current resolution
    resDir <- paste(outWD, "/resampled_", Res, "_sec/", sep ="")
    
    BirdList2<-dir(resDir, pattern = filePattern)
    
    ## Establish and register clusters
    cl<-makeCluster(5, outfile="")
    registerDoSNOW(cl)
    
    ## For each data.frame in current resampled folder: summarize at point
    for(Data in BirdList2){
      
      #     Data<-BirdList2[1]
      
      print(Data)
      bird<-read.table(paste(resDir, Data, sep = ""), sep=",",header=T)
      #writeLines(c(""),"C:/Temp/log.txt")
      
      SEQ<-seq(10,nrow(bird),by=1)
      

      #
      # Call summarize.at.point function
      d<-foreach(Index=SEQ,.combine='rbind',.verbose=TRUE,.errorhandling = 'remove',.packages=c('dplyr','raster','rgdal','sp','adehabitat','adehabitatLT','geosphere'))%dopar%{    #
        #for(Index in SEQ){
        #print(Index)
        #writeLines(as.character(Index),"C:/Temp/log.txt")
        cat(as.character(Index),"\n")
        Summarize.at.point(Index,Data,TmBuff,DistBuff,EndDist,NPoints,SubSamp,Theta,StepSize,TurnSens,species, timezone,outWD)
      }
      
      NewDat<-bird[d[,"Eindex"],]
      
      MM<-cbind(NewDat,d)
      
      ## Create "summarized" folder in output directory if it doesn't exist
      if (!file.exists(paste(resDir, "/summarized/", sep="")))
      {
        dir.create(paste(resDir, "/summarized/", sep=""))
      }
      
      ## Name of output directory.
      name<-paste(resDir, "/summarized/", substr(Data,1,nchar(Data)-4),"_Output.txt", sep="")
      
      ## Write to folder
      write.table(MM,name,sep=",",row.names=F)
      
    }
    
    ## Terminate SNOW cluster
    stopCluster(cl)
    
    
    
  }
  
}
library(zoo)
library(FSA)
library(psych)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(reshape2)
library(dplyr)
library(openair)
??strptime
####HV####
HV<-read.delim("Pool300NoOutliers_timeadded.txt")
headTail(HV)
HV$DateTime<-strptime(HV$DateTime, format="%m/%d/%y %R")
#HV<-HV[order(HV$DateTime),]
HV$Day<-format(HV$DateTime,"%D")
headTail(HV)
HV$DateTime<-as.POSIXct(HV$DateTime)
names(HV)[1] <- "date"
HV$Day<-format(HV$date,"%D")
#subsetting HV df into 1999-2001,2003-2014 to create climatology avoiding bleaching years
HV$year<-format(HV$date,'%Y')
nbHV<-subset(HV, year !='2002' & year !='2003'  & year !='2015' & year !='2016' & year !='2017')
nbHV$year<-factor(nbHV$year)
names(nbHV)[1] <- "date"
nbHV$Day<-format(nbHV$date,"%D")
headTail(nbHV)
#create Night df to calculate daily, then monthly max mean from dd/mm/yyy
Night1<-selectByDate(nbHV, start="14/11/2000",end="31/12/2014",hour = 20:23)
Night2<-selectByDate(nbHV, start="14/11/2000",end="31/12/2014",hour = 00:05)
Night<-rbind(Night1[1:2],Night2[1:2])
head(Night2)
head(Night)
Night<-Night[order(Night$date),]
Night$Day<-format(Night$date,"%D")
Night$year<-format(Night$date,'%Y')
Night<-subset(Night, year !='2002')
Night$year<-factor(Night$year)
levels(Night$year)
#HVdaily calcs from Night df to find MMM
HVdaily<-data.frame("DayRange"=tapply(Night$Pool300, Night$Day, function(x) range(x)[2]-range(x)[1]),"DayMin"=tapply(Night$Pool300, Night$Day, min),"DayMax"=tapply(Night$Pool300, Night$Day, max), "DayMean"=tapply(Night$Pool300, Night$Day, mean))
#HVdaily calcs from nb years to find MMM
HVdaily<-data.frame("DayRange"=tapply(nbHV$Pool300, nbHV$Day, function(x) range(x)[2]-range(x)[1]),"DayMin"=tapply(nbHV$Pool300, nbHV$Day, min),"DayMax"=tapply(nbHV$Pool300, nbHV$Day, max), "DayMean"=tapply(nbHV$Pool300, nbHV$Day, mean))

HVdaily$Date<-strptime(row.names(HVdaily),format="%m/%d/%y")
levels(HVmonthlymean$year)
HVdaily<-HVdaily[order(HVdaily$Date),]
HVdaily$month<-format(HVdaily$Date, "%B-%Y")
HVmonthlymean<-data.frame("MonthlyMeans"=tapply(HVdaily$DayMean,HVdaily$month, mean))
HVmonthlymean$month<-row.names(HVmonthlymean)
HVmonthlymean$year<-format(strptime(paste0("01-",HVmonthlymean$month),format="%d-%B-%Y"),"%Y")
HVmonthlymean$mo<-format(strptime(paste0("01-",HVmonthlymean$month),format="%d-%B-%Y"),"%B")
Summarize(MonthlyMeans~mo, data=HVmonthlymean, digits=3)
#MMM January = 29.153 night nb, 29.432 nb
#Climatology = maximum of the monthly mean SST climatology, mean temperature of the climatologically warmest month at the location
#Hotspot = difference between a nighttime SST value and the corresponding climatology value
#HV df to calculate Hotspots and DHW DHD
HVdaily<-data.frame("DayRange"=tapply(HV$Pool300, HV$Day, function(x) range(x)[2]-range(x)[1]),"DayMin"=tapply(HV$Pool300, HV$Day, min),"DayMax"=tapply(HV$Pool300, HV$Day, max), "DayMean"=tapply(HV$Pool300, HV$Day, mean))
HVdaily$Date<-strptime(row.names(HVdaily),format="%m/%d/%y")
headTail(HVdaily)
HVdaily<-HVdaily[order(HVdaily$Date),]
HVdaily$month<-format(HVdaily$Date, "%B-%Y")
HVmonthlymean<-data.frame("MonthlyMeans"=tapply(HVdaily$DayMean,HVdaily$month, mean))
HVmonthlymean$month<-row.names(HVmonthlymean)
HVmonthlymean$year<-format(strptime(paste0("01-",HVmonthlymean$month),format="%d-%B-%Y"),"%Y")
head(HVmonthlymean)
HVmonthlymean$mo<-format(strptime(paste0("01-",HVmonthlymean$month),format="%d-%B-%Y"),"%B")
Summarize(MonthlyMeans~mo, data=HVmonthlymean, digits=3)

Hotspots<-NULL
for(i in 1:nrow(HVdaily)){
  if((HVdaily[i,"DayMean"]-28.9)>=1){
    z<-HVdaily[i,"DayMean"]-28.9
  }else{
    z<-0
  }
  Hotspots<-c(Hotspots,z)
}
HVdaily$hotspot<-Hotspots
#DHW = the accumulation of HotSpots at that location over a rolling 12-week time period,84 is no of days
HVdaily$DHW<-(1/7)*c(rep(0,83),rollapply(HVdaily$hotspot,84,sum))
HVdaily$DHD<-c(rep(0,83),rollapply(HVdaily$hotspot,84,sum))
HVdaily$year<-format(strptime(paste0("01-",HVdaily$month),format="%d-%B-%Y"),"%Y")
headTail(HVdaily)
write.table(HVdaily,file="2000-2017_CRWMMM_HVdaily_DHW.txt",sep="\t")

#Trying out degree heating minutes? based off of rolling window by row (30min increments)
Hotspots<-NULL
for(i in 1:nrow(HV)){
  if((HV[i,"Pool300"]-29.454)>=1){
    z<-HV[i,"Pool300"]-29.454
  }else{
    z<-0
  }
  Hotspots<-c(Hotspots,z)
}
HV$hotspot<-Hotspots
HV$DHM<-(1/336)*c(rep(0,4031),rollapply(HV$hotspot,4032,sum))
#sum over day and add to HVdaily
DHM<-data.frame(tapply(HV$DHM, HV$Day, sum))
HVdaily<-cbind(HVdaily,DHM)
names(HVdaily)
names(HVdaily)[names(HVdaily) == "tapply.HV.DHM..HV.Day..sum."] <- "DHM"
#monthly stats
HVmonthly<-data.frame("MonthlyMin"=tapply(HVdaily$DayMin,HVdaily$month, min),"MonthlyMax"=tapply(HVdaily$DayMax,HVdaily$month, max),"MonthlyRange"=tapply(HVdaily$DayRange, HVdaily$month, max), "MonthlyMean"=tapply(HVdaily$DayMean,HVdaily$month, mean), "MonthlyDHW"=tapply(HVdaily$DHW ,HVdaily$month, max), "DHWmean"=tapply(HVdaily$DHW ,HVdaily$month, mean))
headTail(HVmonthly)
write.table(HVmonthly,file="2000-2017_nbMMM_HVmonthly_inDHW.txt",sep="\t")

####MV####
MV<-read.delim("Pool400NoOutliers_timeadded.txt")
headTail(MV)
MV$DateTime<-strptime(MV$DateTime, format="%m/%d/%y %R")
MV$Day<-format(MV$DateTime,"%D")
MV$DateTime<-as.POSIXct(MV$DateTime)
#subsetting MV df into 2000-2001,2003-2014 to create climatology avoiding bleaching years
MV$year<-format(MV$DateTime,'%Y')
nbMV<-subset(MV, year !='2002' & year !='2003' & year !='2015' & year !='2016' & year !='2017')
nbMV$year<-factor(nbMV$year)
names(nbMV)[1] <- "date"
nbMV$Day<-format(nbMV$date,"%D")
headTail(nbMV)
#create Night df to calculate daily, then monthly max mean from
MNight1<-selectByDate(nbMV, start="14/11/2000",end="31/12/2014",hour = 20:23)
MNight2<-selectByDate(nbMV, start="14/11/2000",end="31/12/2014",hour = 00:05)
MNight<-rbind(MNight1[1:2],MNight2[1:2])
head(MNight2)
headTail(MNight)
MNight<-MNight[order(MNight$date),]
MNight$Day<-format(MNight$date,"%D")
MNight$year<-format(MNight$date,'%Y')
MNight<-subset(MNight, year !='2002')
MNight$year<-factor(MNight$year)

#MVdaily calcs from Night df to find MMM
MVdaily<-data.frame("DayRange"=tapply(MNight$Pool400, MNight$Day, function(x) range(x)[2]-range(x)[1]),"DayMin"=tapply(MNight$Pool400, MNight$Day, min),"DayMax"=tapply(MNight$Pool400, MNight$Day, max), "DayMean"=tapply(MNight$Pool400, MNight$Day, mean))
#MVdaily calcs from nb years to find MMM
MVdaily<-data.frame("DayRange"=tapply(nbMV$Pool400, nbMV$Day, function(x) range(x)[2]-range(x)[1]),"DayMin"=tapply(nbMV$Pool400, nbMV$Day, min),"DayMax"=tapply(nbMV$Pool400, nbMV$Day, max), "DayMean"=tapply(nbMV$Pool400, nbMV$Day, mean))
MVdaily$Date<-strptime(row.names(MVdaily),format="%m/%d/%y")
headTail(MVdaily)
MVdaily<-MVdaily[order(MVdaily$Date),]
MVdaily$month<-format(MVdaily$Date, "%B-%Y")
MVmonthlymean<-data.frame("MonthlyMeans"=tapply(MVdaily$DayMean,MVdaily$month, mean))
MVmonthlymean$month<-row.names(MVmonthlymean)
MVmonthlymean$year<-format(strptime(paste0("01-",MVmonthlymean$month),format="%d-%B-%Y"),"%Y")
tail(MVmonthlymean)
MVmonthlymean$mo<-format(strptime(paste0("01-",MVmonthlymean$month),format="%d-%B-%Y"),"%B")
Summarize(MonthlyMeans~mo, data=MVmonthlymean, digits=3)
#MMM January= 29.324 nb night, 29.494 nb years 
#MV df to calculate Hotspots and DHW DHD
MVdaily<-data.frame("DayRange"=tapply(MV$Pool400, MV$Day, function(x) range(x)[2]-range(x)[1]),"DayMin"=tapply(MV$Pool400, MV$Day, min),"DayMax"=tapply(MV$Pool400, MV$Day, max), "DayMean"=tapply(MV$Pool400, MV$Day, mean))
MVdaily$Date<-strptime(row.names(MVdaily),format="%m/%d/%y")
headTail(MVdaily)
MVdaily<-MVdaily[order(MVdaily$Date),]
MVdaily$month<-format(MVdaily$Date, "%B-%Y")
MVmonthlymean<-data.frame("MonthlyMeans"=tapply(MVdaily$DayMean,MVdaily$month, mean))
MVmonthlymean$month<-row.names(MVmonthlymean)
MVmonthlymean$year<-format(strptime(paste0("01-",MVmonthlymean$month),format="%d-%B-%Y"),"%Y")
head(MVmonthlymean)
MVmonthlymean$mo<-format(strptime(paste0("01-",MVmonthlymean$month),format="%d-%B-%Y"),"%B")
Summarize(MonthlyMeans~mo, data=MVmonthlymean, digits=3)

Hotspots<-NULL
for(i in 1:nrow(MVdaily)){
  if((MVdaily[i,"DayMean"]-28.9 )>=1){
    z<-MVdaily[i,"DayMean"]-28.9 
  }else{
    z<-0
  }
  Hotspots<-c(Hotspots,z)
}
MVdaily$hotspot<-Hotspots
#DHW = the accumulation of HotSpots at that location over a rolling 12-week time period,84 is no of days
MVdaily$DHW<-(1/7)*c(rep(0,83),rollapply(MVdaily$hotspot,84,sum))
MVdaily$DHD<-c(rep(0,83),rollapply(MVdaily$hotspot,84,sum))
MVdaily$year<-format(strptime(paste0("01-",MVdaily$month),format="%d-%B-%Y"),"%Y")
headTail(MVdaily)
write.table(MVdaily,file="2000-2017_CRWMMM_MVdaily_DHW.txt",sep="\t")

#monthly stats
MVmonthly<-data.frame("MonthlyMin"=tapply(MVdaily$DayMin,MVdaily$month, min),"MonthlyMax"=tapply(MVdaily$DayMax,MVdaily$month, max),"MonthlyRange"=tapply(MVdaily$DayRange, MVdaily$month, max), "MonthlyMean"=tapply(MVdaily$DayMean,MVdaily$month, mean), "MonthlyDHW"=tapply(MVdaily$DHW ,MVdaily$month, max), "DHWmean"=tapply(MVdaily$DHW ,MVdaily$month, mean))
headTail(MVmonthly)
write.table(MVmonthly,file="2000-2017_nbMMM_MVmonthly_inDHW.txt",sep="\t")

####LV####
LV<-read.delim("Pool500NoOutliers_timeadded.txt")
headTail(LV)
LV$DateTime<-strptime(LV$DateTime, format="%m/%d/%y %R")
LV$Day<-format(LV$DateTime,"%D")
LV$DateTime<-as.POSIXct(LV$DateTime)
#subsetting LV df into 2000-2001,2003-2014 to create climatology avoiding bleaching years
LV$year<-format(LV$DateTime,'%Y')
nbLV<-subset(LV, year !='2002' & year !='2003' & year !='2015' & year !='2016' & year !='2017')
nbLV$year<-factor(nbLV$year)
names(nbLV)[1] <- "date"
nbLV$Day<-format(nbLV$date,"%D")
headTail(nbLV)
#create Night df to calculate daily, then monthly max mean from
LNight1<-selectByDate(nbLV, start="14/11/2000",end="31/12/2014",hour = 20:23)
LNight2<-selectByDate(nbLV, start="14/11/2000",end="31/12/2014",hour = 00:05)
LNight<-rbind(LNight1[1:2],LNight2[1:2])
headTail(LNight)
LNight<-LNight[order(LNight$date),]
LNight$Day<-format(LNight$date,"%D")
LNight$year<-format(LNight$date,'%Y')
LNight<-subset(LNight, year !='2002')
LNight$year<-factor(LNight$year)

#LVdaily calcs from Night df to find MMM
LVdaily<-data.frame("DayRange"=tapply(LNight$Pool500, LNight$Day, function(x) range(x)[2]-range(x)[1]),"DayMin"=tapply(LNight$Pool500, LNight$Day, min),"DayMax"=tapply(LNight$Pool500, LNight$Day, max), "DayMean"=tapply(LNight$Pool500, LNight$Day, mean))
#LVdaily calcs from nb years df to find MMM
LVdaily<-data.frame("DayRange"=tapply(nbLV$Pool500, nbLV$Day, function(x) range(x)[2]-range(x)[1]),"DayMin"=tapply(nbLV$Pool500, nbLV$Day, min),"DayMax"=tapply(nbLV$Pool500, nbLV$Day, max), "DayMean"=tapply(nbLV$Pool500, nbLV$Day, mean))
LVdaily$Date<-strptime(row.names(LVdaily),format="%m/%d/%y")
headTail(LVdaily)
LVdaily<-LVdaily[order(LVdaily$Date),]
LVdaily$month<-format(LVdaily$Date, "%B-%Y")
LVmonthlymean<-data.frame("MonthlyMeans"=tapply(LVdaily$DayMean,LVdaily$month, mean))
LVmonthlymean$month<-row.names(LVmonthlymean)
LVmonthlymean$year<-format(strptime(paste0("01-",LVmonthlymean$month),format="%d-%B-%Y"),"%Y")
tail(LVmonthlymean)
LVmonthlymean$mo<-format(strptime(paste0("01-",LVmonthlymean$month),format="%d-%B-%Y"),"%B")
Summarize(MonthlyMeans~mo, data=LVmonthlymean, digits=3)
#MMM January= 29.338 night nb, 29.499 for nb year MMM 
#LV df to calculate Hotspots and DHW DHD
LVdaily<-data.frame("DayRange"=tapply(LV$Pool500, LV$Day, function(x) range(x)[2]-range(x)[1]),"DayMin"=tapply(LV$Pool500, LV$Day, min),"DayMax"=tapply(LV$Pool500, LV$Day, max), "DayMean"=tapply(LV$Pool500, LV$Day, mean))
LVdaily$Date<-strptime(row.names(LVdaily),format="%m/%d/%y")
headTail(LVdaily)
LVdaily<-LVdaily[order(LVdaily$Date),]
LVdaily$month<-format(LVdaily$Date, "%B-%Y")
LVmonthlymean<-data.frame("MonthlyMeans"=tapply(LVdaily$DayMean,LVdaily$month, mean))
LVmonthlymean$month<-row.names(LVmonthlymean)
LVmonthlymean$year<-format(strptime(paste0("01-",LVmonthlymean$month),format="%d-%B-%Y"),"%Y")
head(LVmonthlymean)
LVmonthlymean$mo<-format(strptime(paste0("01-",LVmonthlymean$month),format="%d-%B-%Y"),"%B")
Summarize(MonthlyMeans~mo, data=LVmonthlymean, digits=3)

Hotspots<-NULL
for(i in 1:nrow(LVdaily)){
  if((LVdaily[i,"DayMean"]-28.9)>=1){
    z<-LVdaily[i,"DayMean"]-28.9
  }else{
    z<-0
  }
  Hotspots<-c(Hotspots,z)
}
LVdaily$hotspot<-Hotspots
#DHW = the accumulation of HotSpots at that location over a rolling 12-week time period,84 is no of days
LVdaily$DHW<-(1/7)*c(rep(0,83),rollapply(LVdaily$hotspot,84,sum))
LVdaily$DHD<-c(rep(0,83),rollapply(LVdaily$hotspot,84,sum))
LVdaily$year<-format(strptime(paste0("01-",LVdaily$month),format="%d-%B-%Y"),"%Y")
headTail(LVdaily)
write.table(LVdaily,file="2000-2017_CRWMMM_LVdaily_DHW.txt",sep="\t")

#monthly stats
LVmonthly<-data.frame("MonthlyMin"=tapply(LVdaily$DayMin,LVdaily$month, min),"MonthlyMax"=tapply(LVdaily$DayMax,LVdaily$month, max),"MonthlyRange"=tapply(LVdaily$DayRange, LVdaily$month, max), "MonthlyMean"=tapply(LVdaily$DayMean,LVdaily$month, mean), "MonthlyDHW"=tapply(LVdaily$DHW ,LVdaily$month, max), "DHWmean"=tapply(LVdaily$DHW ,LVdaily$month, mean))
headTail(LVmonthly)
write.table(LVmonthly,file="2000-2017_nbMMM_LVmonthly_inDHW.txt",sep="\t")

####CRW####
time_series_CRW5kmv3.1_1985_2017_AS_OfuIsland_lat14177949s_lon169654364w_col0206_row2083.txt
CRW<-read.delim("time_series_CRW5kmv3.1_2000_2017_AS_OfuIsland_lat14177949s_lon169654364w_col0206_row2083.txt")
headTail(CRW)
#subsetting CRW df into 2000-2001,2003-2014 to create climatology avoiding bleaching years
nbCRW<-subset(CRW, YEAR !='2002' & YEAR !='2003' & YEAR !='2015' & YEAR !='2016' & YEAR !='2017')
headTail(nbCRW)
nbCRW$month<-interaction(nbCRW$M,nbCRW$YEAR)
nbCRW$date<-interaction(nbCRW$M,nbCRW$D,nbCRW$YEAR)
CRWMlymean<-data.frame("MonthlyMeans"=tapply(nbCRW$SST,nbCRW$month, mean))
CRWMlymean$month<-row.names(CRWMlymean)
CRWMlymean$mo<-rep(1:12,rep=13)
Summarize(MonthlyMeans~mo, data=CRWMlymean, digits=3) #max is 29.475 for Apr

#CRW df to calculate Hotspots and DHW DHD
CRW$month<-interaction(CRW$M,CRW$YEAR)
CRW$date<-interaction(CRW$M,CRW$D,CRW$YEAR)
CRWmonthmean<-data.frame("MonthlyMeans"=tapply(CRW$SST,CRW$month, mean))
CRWmonthmean$month<-row.names(CRWmonthmean)
CRWmonthmean$mo<-rep(1:12,rep=18)
CRWmonthmean$year<-rep(2000:2017,each=12)
Summarize(MonthlyMeans~mo, data=CRWmonthmean, digits=3) #

Hotspots<-NULL
for(i in 1:nrow(CRW)){
  if((CRW[i,"SST"]-29.475)>=1){
    z<-CRW[i,"SST"]-29.475
  }else{
    z<-0
  }
  Hotspots<-c(Hotspots,z)
}
CRW$hotspot<-Hotspots
#DHW = the accumulation of HotSpots at that location over a rolling 12-week time period,84 is no of days
CRW$newdhw<-(1/7)*c(rep(0,83),rollapply(CRW$hotspot,84,sum))
CRW$newdhd<-c(rep(0,83),rollapply(CRW$hotspot,84,sum))
headTail(CRW)
CRWmonthly<-data.frame("MonthlyMin"=tapply(CRW$SST,CRW$month, min),"MonthlyMax"=tapply(CRW$SST,CRW$month, max), "MonthlyMean"=tapply(CRW$SST,CRW$month, mean), "Monthlynewdhw"=tapply(CRW$newdhw ,CRW$month, max),"Monthlynewdhwmean"=tapply(CRW$newdhd,CRW$month, mean), "MonthlyDHW"=tapply(CRW$DHW ,CRW$month, max))
CRWdaily<-data.frame("Date"=CRW$date,"DayMean"=CRW$SST, "Dailynewdhw"=CRW$newdhw,"Dailynewdhd"=CRW$newdhd, "DHW"=CRW$DHW)

write.table(CRWdaily,file="2000-2017_dailyCRW5km_DHW.txt",sep="\t")

####All data####
# AllData<-list(HV,MV,LV)
# AllDataSummary<-data.frame(sapply(AllData,function(x) summary(x$Temp)))
# names(AllDataSummary)<-Depths
# AllDataSummary[6,]-AllDataSummary[1,]

Sites<-c("HV","MV","LV")
headTail(LVdaily)
AllDailys<-list(HVdaily,MVdaily,LVdaily)
summary(AllDailys)
DayRangeSummary<-data.frame(sapply(AllDailys,function(x) summary(x$DayRange)))
names(DayRangeSummary)<-Sites
DayMeanSummary<-data.frame(sapply(AllDailys,function(x) summary(x$DayMean)))
names(DayMeanSummary)<-Sites
DayMinSummary<-data.frame(sapply(AllDailys,function(x) summary(x$DayMin)))
names(DayMinSummary)<-Sites
DayMaxSummary<-data.frame(sapply(AllDailys,function(x) summary(x$DayMax)))
names(DayMaxSummary)<-Sites
#need to not summarize but cumulative count
DayDHWSummary<-data.frame(sapply(AllDailys,function(x) summary(x$DHW)))
names(DayDHWSummary)<-Sites
#need to not summarize but cumulative count
DayDHDSummary<-data.frame(sapply(AllDailys,function(x) summary(x$DHD)))
names(DayDHDSummary)<-Sites

DayRangeSummary
DayMeanSummary
DayMinSummary
DayMaxSummary
DayDHWSummary
DayDHDSummary

AllDailysMerge<-cbind(AllDailys[[1]],"Site"=Sites[1])
for(i in 2:5){
  AllDailysMerge<-rbind(AllDailysMerge,cbind(AllDailys[[i]],"Site"=Sites[i]))	
}

names(AllDailysMerge)[6] <- "date"
AllDailysMerge$date<-as.POSIXct(AllDailysMerge$date)
AllDailyM<-selectByDate(AllDailysMerge, start="01/01/2010",end="31/12/2017")
str(AllDailyM)
write.table(AllDailysMerge,file="2000-2017_nbMMMinsituDHW_daily_allsites.txt",sep="\t")

AllMonthlys<-list(HVmonthly,MVmonthly,LVmonthly)
AllMonthlysMerge<-cbind(AllMonthlys[[1]],"Sites"=Sites[1])
names(AllMonthlysMerge)


####Wilcox-rank sum####
dat<-read.delim("2000-2017_CRWinsituDHW_daily.txt")
summary(dat)
#separate into sites
HV<-dat[dat$Sites=='HV',]
dim(HV) #5713 13
MV<-dat[dat$Sites=='MV',]
dim(MV) #5536 13
LV<-dat[dat$Sites=='LV',]
dim(LV) #5476 13
#keep rows with matching dates since each df is different in length, all at 5053 final rownum
HV.5476<-match_df(HV,LV,on="date")
MV.5476<-match_df(MV,LV,on="date")
MV.5262<-match_df(MV.5476,HV.5476,on="date")
HV.5262<-match_df(HV.5476,MV.5262,on="date")
LV.5262<-match_df(LV,MV.5262,on="date")
HVvMV<-rbind.data.frame(HV.5262,MV.5262)
MVvLV<-rbind.data.frame(MV.5262,LV.5262)
HVvLV<-rbind.data.frame(HV.5262,LV.5262)
#create the parameters for the sliding window
slidingWindow=seq(1,1433,by=7)
#The following function takes the data and melts it so that the quantitative factor of choice can be selected 
#for a two week time period both years and then tested with a Wilcox sign-rank test and finally return a dataframe 
#with all the p-values for the tests for each window.
wilcoxvalues <- lapply(slidingWindow, function(i) {
  HVvMV.sub=HVvMV[i:(i+13), ]
  HVvMV.sub=melt(HVvMV.sub, id.vars=c("date","month","year"))
  HVvMV.sub=na.omit(HVvMV.sub)})
  HVvMV.sub$variable=as.factor(HVvMV$variable)
  data.frame(week=paste0("Week: ", i%/%7+1, "-", i%/%7+2),
             p.values=wilcox.test(inDHW~Sites, HVvMV.sub)$p.value)
})

wilcoxdf.inDHW<-do.call(rbind, wilcoxvalues.inDHW)


HVwind<-rollapply(HV$inDHW, width = 14, by = 7, FUN = sum, na.rm = TRUE)
dim(LV)
5713/7

MVwind<-rollapply(MV$inDHW, width = 14, by = 7, FUN = sum, na.rm = TRUE)

LVwind<-rollapply(LV$inDHW, width = 14, by = 7, FUN = sum, na.rm = TRUE)
?rollapply


WindowAnalysis<-data.frame("Window"=numeric(), "HVsum"=numeric, "MVsum"=numeric(), "Pvalue"=numeric())
for(i in c(1:781)){
  StartIdx=i*7
  EndIdx=i*7+14
  HVMVWilTest<-wilcox.test(HVDHWs ~ MV$DHWs, data=DHWs[StartIdx:EndIdx,])
  HVLVWilTest<-wilcox.test(HVDHWs ~ MV$DHWs, data=DHWs[StartIdx:EndIdx,])
  MVLVWilTest<-wilcox.test(HVDHWs ~ MV$DHWs, data=DHWs[StartIdx:EndIdx,])
  WindowAnalysis[i,]<-data.frame("Window"=i, "HVsum"=sum(DHWs$HVDHWs[StartIdx:EndIdx]), "MVsum"=sum(DHWs$MVDHWs[StartIdx:EndIdx]), "LVsum"=, "Pvalue"=WilTest$pvalue)
#WindowAnalysis<-rbind(WindowAnalysis, c())
}

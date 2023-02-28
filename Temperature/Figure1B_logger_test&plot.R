#### This is an R script that manipulates insitu HOBO temp data for descriptive statistics, ANoVA, correlations, and plots
#Written by Courtney Klepac (with edits/assistance from Dan Barshis)

library(FSA)
library(dplyr)
library(lubridate)
library(tidyr)
library(zoo)
library(expss)
library(emmeans)
library(psych)
library(PerformanceAnalytics)
library(cowplot)
library(ggplot2)

#read in hobotemp file
temp<-read.delim("AStemp2016-17_avg2.txt",header=T)
#temp<-read.delim(pipe("pbpaste"))
head(temp)
# to convert the DateTime column to and R DateTime format see ?strptime() for more explanation
temp$DateTime<-strptime(temp$DateTime, format="%m/%d/%y %H:%M")
#make month and year-month column to do summarize stats
temp$year=year(temp$DateTime)
temp$month=month(temp$DateTime)
temp$day<-format(temp$DateTime,"%D")
temp<-unite(temp, 'date', c('month','year'),remove=F)
temp<-unite(temp, 'd', c('day','month','year'),remove=F)
head(temp)
####daily descriptive ststs####
# use the day to calculate daily stats (range, min, max, mean)
tempdaily<-data.frame("HVdayRange"=tapply(temp$HV, temp$day, function(x) range(x)[2]-range(x)[1]),"HVdayMin"=tapply(temp$HV, temp$day, min),"HVdayMax"=tapply(temp$HV, temp$day, max), "HVdayMean"=tapply(temp$HV, temp$day, mean),
                      "MVdayRange"=tapply(temp$MV, temp$day, function(x) range(x)[2]-range(x)[1]),"MVdayMin"=tapply(temp$MV, temp$day, min),"MVdayMax"=tapply(temp$MV, temp$day, max), "MVdayMean"=tapply(temp$MV, temp$day, mean),       
                      "LVdayRange"=tapply(temp$LV, temp$day, function(x) range(x)[2]-range(x)[1]),"LVdayMin"=tapply(temp$LV, temp$day, min),"LVdayMax"=tapply(temp$LV, temp$day, max), "LVdayMean"=tapply(temp$LV, temp$day, mean))
head(tempdaily)

#summarize max, min, mean, etc for each day
# HVsumm<-Summarize(HV~d,data=temp,digits=3)
# MVsumm<-Summarize(MV~d,data=temp,digits=3)
# LVsumm<-Summarize(LV~d,data=temp,digits=3)
# write.table(HVsumm,"ASdailytempsumm2-HV_2016-17.txt",sep = "\t")
# write.table(MVsumm,"ASdailytempsumm2-MV_2016-17.txt",sep = "\t")
# write.table(LVsumm,"ASdailytempsumm2-LV_2016-17.txt",sep = "\t")

# print out summary stats for daily measurements
tempdaily$date<-row.names(tempdaily)
headTail(tempdaily)
tempdaily$date<-strptime(tempdaily$date, tz="GMT", format='%m/%d/%y')
tempdaily<-tempdaily[order(tempdaily$date),]
headTail(tempdaily)

write.table(tempdaily,"ASdailytempsumm2_2016-17.txt",sep = "\t")

#boxplot of daily stats, do for each site
par(mfrow=c(1,3))
boxplot(tempdaily$LVdayRange, main = "Daily temperature Range °C",ylab = "temperature °C",notch=T, col="blue",cex.lab=1.5,cex.axis=1.5)
boxplot(tempdaily$LVdayMean, main = "Daily temperature Mean °C", ylab = "temperature °C",notch=T, col="blue",cex.lab=1.5,cex.axis=1.75)
boxplot(tempdaily$LVdayMax, main = "Daily temperature Max °C", ylab = "temperature °C",notch=T, col="blue",cex.lab=1.5,cex.axis=1.8)

#plot of entire temperature trace
#use the tickpos to set the axis labels you want to plot (note these need to match your first date-time file)
tickpos<-seq(as.POSIXct("2016-01-01"),as.POSIXct("2017-01-29"),by="1 month")
quartz()
pdf("AS2016-17_meantemptrace.pdf",14, 7)
plot(tempdaily$date,tempdaily$HVdayMean, type="l", lwd= 1.5, ylab="Mean Water temperature °C", xlab="Date", xaxt='n',ylim=c(27,32), col="red")
# use points to add other temp traces to the same plot, will automatically use the date-time to plot them together
points(tempdaily$date,tempdaily$MVdayMean, type="l", lwd=1.5, col="gold")
points(tempdaily$date, tempdaily$LVdayMean, type="l", lwd=1.5, col="blue")
abline(h=30.2,lty=2,col="black")
legend("bottomright",c("HV","MV","LV"),lty=1,lwd=2.5,col=c("red","gold","blue"),bty="n")
#to add your specific axis format to the x-axis
axis.POSIXct(side=1, at=tickpos, format="%Y-%b")
dev.off()

#plot Ofu climatology (SST and DHW) and site temps
tempdaily<-read.delim("AS_climatology_2016-17.txt")
tempdaily$date<-strptime(tempdaily$date, tz="GMT", format='%m/%d/%y')
tempdaily<-tempdaily[order(tempdaily$date),]
head(tempdaily)
pdf("AS2016-17_climatology.pdf",14, 7)
par(mar=c(5,5,3,5))
plot(tempdaily$date,tempdaily$HVDayMean, type="n", cex.lab=1.5,cex.axis=1.5, ylab="Water temperature (°C)", xlab="Month", xaxt='n',ylim=c(25,32), main="Ofu Pool temperature & DHW")
# use points to add other temp traces to the same plot, will automatically use the date-time to plot them together
points(tempdaily$date,tempdaily$HVDayMean, type="l", lwd=2, col="red")
points(tempdaily$date,tempdaily$MVDayMean, type="l", lwd=2, col="gold")
points(tempdaily$date, tempdaily$LVDayMean, type="l", lwd=2, col="blue")
#points(tempdaily$date, tempdaily$SST, type="l", lwd=1.5, col="gray55")
abline(h=30.2,lty=3,col="black")
legend("right",c("HV","MV","LV","Ofu SST","Ofu DHW"),lty=c(1,1,1,1,2),lwd=2.5,col=c("red","gold","blue","gray55","black"),bty="n")
#add in DHW with separate yaxis
par(new=T)
plot(tempdaily$date,tempdaily[,"DHW"],type="l",lty=2, col="black", lwd=2, xaxt='n',yaxt='n',xlab='',ylab='', ylim=c(0,15))
axis(side=4, cex.axis=1.5)
mtext("DHW (°C week)",cex=1.5,side=4,line=3)
#to add your specific axis format to the x-axis
tickpos<-seq(as.POSIXct("2016-01-01"),as.POSIXct("2017-05-30"),by="1 month")
axis.POSIXct(side=1, at=tickpos, format="%Y-%b")
dev.off()


####monthly stats####
# make a new column for just the month-year in R format
temp$MonthYear<-format(temp$DateTime,"%B-%Y")
dim(tempMonthly)
head(temp)
# use the month to calculate monthly stats (range, min, max, mean)
tempMonthly<-data.frame("HVMonthYearRange"=tapply(temp$HV, temp$MonthYear, function(x) range(x)[2]-range(x)[1]),"HVMonthYearMin"=tapply(temp$HV, temp$MonthYear, min),"HVMonthYearMax"=tapply(temp$HV, temp$MonthYear, max), "HVMonthYearMean"=tapply(temp$HV, temp$MonthYear, mean),
                        "MVMonthYearRange"=tapply(temp$MV, temp$MonthYear, function(x) range(x)[2]-range(x)[1]),"MVMonthYearMin"=tapply(temp$MV, temp$MonthYear, min),"MVMonthYearMax"=tapply(temp$MV, temp$MonthYear, max), "MVMonthYearMean"=tapply(temp$MV, temp$MonthYear, mean),
                        "LVMonthYearRange"=tapply(temp$LV, temp$MonthYear, function(x) range(x)[2]-range(x)[1]),"LVMonthYearMin"=tapply(temp$LV, temp$MonthYear, min),"LVMonthYearMax"=tapply(temp$LV, temp$MonthYear, max), "LVMonthYearMean"=tapply(temp$LV, temp$MonthYear, mean))
# #summarize max, min, mean, etc for each month
# HVmo<-Summarize(HV~date,data=temp,digits=3)
# MVmo<-Summarize(MV~date,data=temp,digits=3)
# LVmo<-Summarize(LV~date,data=temp,digits=3)

# print out summary stats for daily measurements
headTail(tempMonthly)
write.table(tempMonthly,"ASmonthlysumm_2016-17.txt",sep = "\t")

####AOV of in situ temp####
#anova of max temp
maxmodel<-lm(max~site*season,data=temp)
anova(maxmodel)

#least-sq means comparisons
max <- emmeans(maxmodel, pairwise ~ site|season, weights = "proportional", adjust="none")
rbind(max$contrasts)
max <- emmeans(maxmodel, pairwise ~ season|site, weights = "proportional", adjust="none")
rbind(max$contrasts)
max <- emmeans(maxmodel, pairwise ~ site, weights = "proportional", adjust="none")
rbind(max$contrasts)

##old way of pairwise comparisons, likely underestimates since computes all possible comparisons, even nonsense ones
#lsmeans(model,pairwise~site*season,adjust='tukey')
#lsmeans(model,pairwise~season,adjust='tukey')

#anova of mean temp
meanmodel=lm(mean~site*season,data=temp)
anova(meanmodel)

mean <- emmeans(meanmodel, pairwise ~ site|season, weights = "proportional", adjust="none")
rbind(mean$contrasts)
mean <- emmeans(meanmodel, pairwise ~ season|site, weights = "proportional", adjust="none")
summary(mean$emmeans)
rbind(mean$contrasts)
mean <- emmeans(meanmodel, pairwise ~ site, weights = "proportional", adjust="none")
summary(mean$emmeans)
rbind(mean$contrasts)

#anova of min temp
minmodel=lm(min~site*season,data=temp)
anova(minmodel)

min <- emmeans(minmodel, pairwise ~ site|season, weights = "proportional", adjust="none")
rbind(min$contrasts)
min <- emmeans(minmodel, pairwise ~ season|site, weights = "proportional", adjust="none")
rbind(min$contrasts)
min <- emmeans(minmodel, pairwise ~ site, weights = "proportional", adjust="none")
rbind(min$contrasts)

#anova of dtr
dtrmodel=lm(dtr~site*season,data=temp)
anova(dtrmodel)
dtr <- emmeans(dtrmodel, pairwise ~ site|season, weights = "proportional", adjust="none")
rbind(dtr$contrasts)
dtr <- emmeans(dtrmodel, pairwise ~ season|site, weights = "proportional", adjust="none")
rbind(dtr$contrasts)
dtr <- emmeans(dtrmodel, pairwise ~ site, weights = "proportional", adjust="none")
rbind(dtr$contrasts)

####calculate DR90 and days over threshold####
#Excel - df manipulation from dailysumm temp file, cbind sites and add in season column#
temp2<-read.delim("AS_tempsumday2_2016-17.txt",header=T)
headTail(temp2)
#subset new df based on site
HV<-temp2[temp2$site=='HV',]
MV<-temp2[temp2$site=='MV',]
LV<-temp2[temp2$site=='LV',]

#calculate DR90%, the 90th quartile for daily temp ranges over rolling 10day
#add to tempsumday in Excel, have to copy values for rolling 10day
head(HVDR90)
HVDR90<-rollapply(HV$dtr,width=10,FUN="quantile",p=0.9,by=10,na.rm=TRUE)
MVDR90<-rollapply(MV$dtr,width=10,FUN="quantile",p=0.9,by=10)
LVDR90<-rollapply(LV$dtr,width=10,FUN="quantile",p=0.9,by=10)

#loop over number of rows by date to count number of days over 30.1
#also have to copy values for each row
sum(count_row_if(gt(30.2),HV$max)) #122
sum(count_row_if(gt(30.2),MV$max)) #128
sum(count_row_if(gt(30.2),LV$max)) #106

sum(count_row_if(gt(31),HV$max)) #69
sum(count_row_if(gt(31),MV$max)) #72
sum(count_row_if(gt(31),LV$max)) #38

sum(count_row_if(gt(32),HV$max)) #33
sum(count_row_if(gt(32),MV$max)) #27
sum(count_row_if(gt(32),LV$max)) #8

sum(count_row_if(gt(33),HV$max)) #5
sum(count_row_if(gt(33),MV$max)) #0
sum(count_row_if(gt(33),LV$max)) #0

#what about in each season? Subrract summer sum from total to get winter values
#subset new df based on site
HVsumm<-HV[HV$season=='Summer',]
MVsumm<-MV[MV$season=='Summer',]
LVsumm<-LV[LV$season=='Summer',]

sum(count_row_if(gt(30.2),HVsumm$max)) #93
sum(count_row_if(gt(30.2),MVsumm$max)) #93
sum(count_row_if(gt(30.2),LVsumm$max)) #85

sum(count_row_if(gt(31),HVsumm$max)) #60
sum(count_row_if(gt(31),MVsumm$max)) #58
sum(count_row_if(gt(31),LVsumm$max)) #35

sum(count_row_if(gt(32),HVsumm$max)) #29
sum(count_row_if(gt(32),MVsumm$max)) #24
sum(count_row_if(gt(32),LVsumm$max)) #8

sum(count_row_if(gt(33),HVsumm$max)) #5
sum(count_row_if(gt(33),MVsumm$max)) #0
sum(count_row_if(gt(33),LVsumm$max)) #0


####correlation#####
#use modified tempsumday file from above which has DR90 and days over threshold added in
temp3<-read.delim(pipe("pbpaste"))

names(temp3)
temp.num=dplyr::select(temp,mean,max,min,dtr,dr90,days.31,days.32)
headTail(temp.num)
corr.test(temp.num,
          use = "pairwise",
          method="pearson",
          adjust="none",     # Can adjust p-values; see ?p.adjust for options
          alpha=.05)

chart.Correlation(temp.num,
                  method="pearson",
                  histogram=TRUE,
                  pch=16)

# scatter<-ggscatter(temp, x = "heatchl", y = "lightchl", 
#                    add = "reg.line", conf.int = TRUE, 
#                    cor.coef = TRUE, cor.method = "pearson",
#                    xlab = "", ylab = "")
# ggsave("./cortest-natVexp_20191212.pdf", width = 10, height = 6)

####plotting summary of mean, max, dtr ####
summary(temp2)
temp2$site=factor(temp2$site,levels(temp2$site)[c(1,3,2)])
#seasonal 
#boxplot of mean
pdf("20200701_2016-17_OfuSeasonboxplot-outliersmean.pdf",10,8)

p3<-ggplot(temp2,aes(x=season, y=mean, fill=site)) + geom_boxplot() + stat_summary(fun.y=mean, geom="point", shape=18, size=3, position = position_dodge(0.75)) +       
  ylab("Mean Daily Temperature (°C)") + ggtitle("C") + scale_fill_manual(guide=FALSE,values=c("red","gold","blue")) + theme_bw(base_size=14) +theme(axis.title.x=element_blank())

p1<-ggplot(temp2,aes(x=season, y=max, fill=site)) + geom_boxplot() +  stat_summary(fun.y=mean, geom="point", shape=18, size=3, position = position_dodge(0.75)) +              
  ylab("Max Daily Temperature (°C)") + ggtitle("A") + scale_fill_manual(guide=FALSE,values=c("red","gold","blue")) + theme_bw(base_size=14) +theme(axis.title.x=element_blank())  

p2<-ggplot(temp2,aes(x=season, y=dtr, fill=site)) + geom_boxplot() +  stat_summary(fun.y=mean, geom="point", shape=18, size=3, position = position_dodge(0.75)) +              
  ylab("Daily Temperature Range (°C)") + ggtitle("B")+ scale_fill_manual(guide=FALSE,values=c("red","gold","blue")) + theme_bw(base_size=14) +theme(axis.title.x=element_blank()) 
plot_grid(p1, p2, p3, nrow=1)    
#save_plot("20190221_2015-16_OfuSeasonboxplot-outliers-6.pdf", gg_all, base_height = 6)
dev.off()   

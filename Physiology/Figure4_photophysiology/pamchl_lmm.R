# Ofu Island por2 CBASS - written by Courtney Klepac
# Physiological measurements - PAM & chl

#libraries
library(car)
library(lmerTest)
library(emmeans)
library(ggplot2)
library(ggpubr) 
library(Hmisc)
library(dplyr)
library(FSA)
library(lme4)
library(ggplot2)
library(car)
library(multcomp)
library(gtools)
library(plyr)
library(PerformanceAnalytics)
library(plotrix)

#### data input and organize####
dat=read.csv("AS_por2.csv", header=T)
colnames(dat)
dat$tank <- as.factor(dat$tank)
dat$origin_dest=as.factor(dat$origin_dest)
dat$origin=as.factor(dat$origin)
dat$dest=as.factor(dat$dest)
dat$time=as.factor(dat$time)
dat$trt=as.factor(dat$trt)
dat$colony=as.factor(dat$colony)


#set factor levels
dat$origin_dest=factor(dat$origin_dest,levels(dat$origin_dest)[c(1,5,3,2,6,4)])
dat$dest=factor(dat$dest,levels(dat$dest)[c(1,3,2)])
dat$origin=factor(dat$origin,levels(dat$origin)[c(1,3,2)])
dat$time=factor(dat$time,levels(dat$time)[c(1,3,2)])

levels(dat$time)[levels(dat$time)=="1mo"] <- "Aug-2016"
levels(dat$time)[levels(dat$time)=="6mo"] <- "Feb-2017"
levels(dat$time)[levels(dat$time)=="24mo"] <- "June-2018"
summary(dat)
tail(dat)

#removing 24mo from df
dat<-subset(dat, time !='24mo')
dat$time<-factor(dat$time)
#removing non-assayed samples
dat<-subset(dat, trt !="NA")
dat$trt<-factor(dat$trt)

#creating interactions for post hoc comparisons
dat$grtime<-interaction(dat$time,dat$origin_dest) 
dat$oritime<-interaction(dat$origin,dat$time)
dat$desttime<-interaction(dat$dest,dat$time)
dat$grtrt<-interaction(dat$origin_dest,dat$trt)
names(dat)

#creating heat and cont datasets
heat<-dat[dat$trt=="heat",]
summary(heat)
cont<-dat[dat$trt=="cont",]
#calculate max quant yield, Fv/Fm retained after heat stress
heat$mqy<-(1-heat$ylossnorm)
dat$mqy<-(1-dat$ylossnorm)
#calculated total chl retained after heat stress
heat$chlrat<-(1-(cont$totchl-heat$totchl)/cont$totchl)

####PAM####
# Statistical tests for maximum quantum yield retained at end of assay 
# Analysis of variance (2-way) ANOVA dat (mixed linear model)
# Mixed Effects model origin_dest- random effect of colony nested in tank
# Since the model is completely balanced the SS type I, II or III will provide the same result
pam.model <- lmer(mqy ~ origin_dest * time  + (1|tank/colony) , data = heat)
summary(pam.model)
anova(pam.model)
rand(pam.model)
confint(pam.model,oldNames=F)

# Model fitting and assumptions diagnostic 
x = residuals(pam.model)
plot(fitted(pam.model), x) 
#leveneTest(x ~ grtrt * time * tank* colony, data=dat, center=mean) # formal statistical test for homogeinity of variance (not recommended due the small sample size)
plot(pam.model2) # Residual vs Fitted values
qqnorm(x); qqline(x) # qq plot to check for normal distribution of residuals
hist(x) # histogram of residuals to check for normal distribution of residuals
shapiro.test(x) # formal statistical test (not recommended due the small sample size)

# comparing between reef sites within each time treatment
pam.emms.reef <- emmeans(pam.model, pairwise ~ origin, weights = "proportional", adjust="none")
summary(pam.emms.reef$emmeans)
# P.value adjustment of the Bonferroni
rbind(pam.emms.reef$contrasts, adjust="bonferroni")

# Statistical tests for fv/fm yloss normalized ((0-21/0)) including treatment
# Analysis of variance (2-way) ANOVA dat (mixed linear model)
# Mixed Effects model origin_dest- random effect of colony nested in tank
pam.model2 <- lmer(ylossnorm ~ origin_dest *time *trt + (1|tank/colony), data = dat)
summary(pam.model2)
anova(pam.model2)
rand(pam.model2)
confint(pam.model2,oldNames=F)

# Model fitting and assumptions diagnostic 
x = residuals(pam.model2)
plot(fitted(pam.model2), x) 
#leveneTest(x ~ grtrt * time * tank* colony, data=dat, center=mean) # formal statistical test for homogeinity of variance (not recommended due the small sample size)
plot(pam.model2) # Residual vs Fitted values
qqnorm(x); qqline(x) # qq plot to check for normal distribution of residuals
hist(x) # histogram of residuals to check for normal distribution of residuals
shapiro.test(x) # formal statistical test (not recommended due the small sample size)

# comparing between reef sites within each time treatment
pam.emms.reef <- emmeans(pam.model2, pairwise ~ trt|time, weights = "proportional", adjust="none")
# P.value adjustment of the Bonferroni
rbind(pam.emms.reef$contrasts, adjust="bonferroni")
# comparing between reef sites within each time treatment
pam.emms.reef <- emmeans(pam.model2, pairwise ~ time|origin_dest|trt, weights = "proportional", adjust="none")
# P.value adjustment of the Bonferroni
rbind(pam.emms.reef$contrasts, adjust="bonferroni")

# hr21 Fv/Fm only
# Mixed Effects model origin and dest- random effect of colony nested in tank
# Since the model is completely balanced the SS type I, II or III will provide the same result
pam.model3 <- lmer(hr21 ~ origin_dest * time *trt + (1|tank/colony), data = dat)
summary(pam.model3)
anova(pam.model3)
rand(pam.model3)
confint(pam.model3,oldNames=F)

# Model fitting and assumptions diagnostic 
x = residuals(pam.model3)
plot(fitted(pam.model3), x) 
#leveneTest(x ~ grtrt * time * tank* colony, data=dat, center=mean) # formal statistical test for homogeinity of variance (not recommended due the small sample size)
plot(pam.model3) # Residual vs Fitted values
qqnorm(x); qqline(x) # qq plot to check for normal distribution of residuals
hist(x) # histogram of residuals to check for normal distribution of residuals
shapiro.test(x) # formal statistical test (not recommended due the small sample size)

# comparing between reef sites within each time treatment
pam.emms.reef <- emmeans(pam.model, pairwise ~ trt, weights = "proportional", adjust="none")
# P.value adjustment of the Bonferroni
rbind(pam.emms.reef$contrasts, adjust="none")

####PAMplots#####
#Dan's summary stat calcs
#pammean=tapply(dat$hr21,list(dat$grtrt,dat$time), mean)
#pamSD=tapply(dat$hr21,list(dat$grtrt,dat$time), sd)
#pamN=tapply(dat$hr21,list(dat$grtrt,dat$time), length)
#pamSE=pamSD/sqrt(pamN)

###hr21
sum<-Summarize(hr21~grtrt+time, data=dat, digits=3)
sum$se = sum$sd/sqrt(sum$n)
sum$CI_lower2 <- sum$mean - qnorm(0.975)*sum$se
sum$CI_upper2 <- sum$mean + qnorm(0.975)*sum$se
sum$Origin<-factor(c("HV","MV","LV","HV","MV","LV"))
sum$dest<-factor(c("HV","HV","HV","MV","MV","LV"))
sum$Time=factor(c("Aug 2016","Aug 2016","Aug 2016","Aug 2016","Aug 2016","Aug 2016","Aug 2016","Aug 2016","Aug 2016","Aug 2016","Aug 2016","Aug 2016","Feb 2017","Feb 2017","Feb 2017","Feb 2017","Feb 2017", "Feb 2017","Feb 2017","Feb 2017","Feb 2017","Feb 2017","Feb 2017", "Feb 2017"))
sum$Origin=factor(sum$Origin,levels=unique(sum$Origin))
sum$dest=factor(sum$dest,levels=unique(sum$dest))
sum$Treatment=factor(c("cont","cont","cont","cont","cont","cont","heat","heat","heat","heat","heat","heat"))

#points with sd or ci 95% bars, comment out which you want
pd=position_dodge(.75)
plot.pam<-ggplot(sum,aes(x=dest,y=mean,fill=Origin,color=Treatment)) + 
  geom_point(aes(shape=Origin),size=3,position=pd) +
  #geom_errorbar(aes(ymin=mean - sd,ymax=mean + sd), width=0.2,position=pd) +
  geom_errorbar(aes(ymin=CI_lower2, ymax=CI_upper2, color=Treatment), width=.2, position=pd) + 
  facet_wrap(~time) +
  theme_bw() + 
  theme(plot.margin=margin(0.5,0.5,0.5,0.5,'cm'), panel.grid.minor=element_blank()) +
  scale_colour_manual(values=c("gray","black")) +
  coord_cartesian(xlim=c(0.5,3.4), ylim=c(0,0.6), expand=F) + 
  ylab("Max Quantum Yield (Fv/Fm)\n after Acute Heat Stress") + 
  xlab("Transplant Site")
save_plot("./PAM/plots/20201006_porII_fvfm-trt-ci.pdf", plot.pam,
          base_aspect_ratio = 1.3)

#plotting Fv/Fm normalized to t=0, both heat and control, with sdev bars
sum=Summarize(ylossnorm~grtrt+time, data=dat, digits=3)
sum$se = sum$sd / sqrt(sum$n)
sum$se = signif(sum$se, digits=3)
sum$time=factor(sum$time,levels=unique(sum$time))
sum$Origin<-factor(c("HV","MV","LV","HV","MV","LV","HV","MV","LV","HV","MV","LV","HV","MV","LV","HV","MV","LV","HV","MV","LV","HV","MV","LV"))
sum$dest<-factor(c("HV","HV","HV","MV","MV","LV","HV","HV","HV","MV","MV","LV","HV","HV","HV","MV","MV","LV","HV","HV","HV","MV","MV","LV"))
sum$Origin=factor(sum$Origin,levels=unique(sum$Origin))
sum$dest=factor(sum$dest,levels=unique(sum$dest))
sum$trt=factor(c("cont","cont","cont","cont","cont","cont","heat","heat","heat","heat","heat","heat","cont","cont","cont","cont","cont","cont","heat","heat","heat","heat","heat","heat"))
sum

pd=position_dodge(.75)
plot.pam<-ggplot(sum,aes(x=dest,y=mean,fill=Origin,color=trt)) + 
  geom_point(aes(shape=Origin),size=3,position=pd) +
  geom_errorbar(aes(ymin=mean - sd,ymax=mean + sd), width=0.2,position=pd) + 
  facet_wrap(~time) +
  theme_bw() + 
  theme(plot.margin=margin(0.5,0.5,0.5,0.5,'cm'), panel.grid.minor=element_blank()) +
  scale_colour_manual(values=c("gray","black")) +
  coord_cartesian(xlim=c(0.5,3.4), ylim=c(0,0.6), expand=F) + 
  ylab("Loss in Max Quantum Yield (Fv/Fm)\n after Acute Heat Stress") + 
  xlab("Transplant Site")
save_plot("./plots/20200305_porII_ylossnorm-trt-sd.pdf", plot.pam,
          base_aspect_ratio = 1.3)

#boxplot
pam <- ggplot(data=dat, 
              aes(x=dest, y=mqy, label= time, fill=origin, color=trt)) +
  scale_fill_manual(values = c ("red", "gold", "blue"), name = "Reef site") +
  scale_color_manual(values = c('gray','black'), name = "Treatment") +
  stat_boxplot(geom ='errorbar', width = 0.5, lwd=0.5)+
  geom_boxplot(width=0.5, lwd=0.5, fatten=1) +
  expand_limits(y = 0)+
  facet_grid(~time, space = "free", scales = "free")+ #this can also be used but rows and columns have to be specified -> facet_grid(. ~ experiment)
  theme_bw() +
  #annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, colour = "black", size=2) +
  theme(line= element_line(size = 0.75),
        axis.line = element_line(colour = "grey20"),
        axis.ticks.length = unit(0.2 , "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(), 
        strip.text.x = element_text(color = "grey20", size = 16, angle = 0, hjust = 0, vjust = 0.5, face = "plain"),
        panel.spacing = unit(2, "lines"))
pam + xlab(label = "Transplant Site") + ylab(label = "Proportion of Photosynthetic efficiency (Fv/Fm) Retained")+
  theme(axis.text.x = element_text(color = "grey20", size = 13, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 13, angle = 0, hjust = .5, vjust = .5, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 14, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 14, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        legend.title = element_text(colour="grey20", size=14, face="bold"),
        legend.text = element_text(colour="grey20", size=13, face="plain"),
        legend.position="bottom")
ggsave("./PAM/plots/PAMmqy-trt_boxplot_16122019.pdf", width = 10, height = 6)


####CHL####
# Statistical tests for total chl ratio
# Analysis of variance (2-way) ANOVA (mixed linear model)
# Mixed Effects model origin and dest separate- random effect of colony nested in tank
# Since the model is completely balanced the SS type I, II or III will provide the same result
chl.model <- lmer(chlrat ~ origin_dest * time + (1|tank/colony), data = heat)
summary(chl.model)
anova(chl.model)
rand(chl.model)
confint(chl.model,oldNames=F)

# Model fitting and assumptions diagnostic 
x = residuals(chl.model)
plot(fitted(chl.model), x) 
for(i in c("Aug-2016","Feb-2017")){
  print(paste0("Time",i))
  print(bartlett.test(chlrat ~ origin_dest, data=heat[heat$time==i,]))
}
plot(chl.model) # Residual vs Fitted values
qqnorm(x); qqline(x) # qq plot to check for normal distribution of residuals
hist(x) # histogram of residuals to check for normal distribution of residuals
shapiro.test(x) # formal statistical test (not recommended due the small sample size)

# comparing between reef sites within each time treatment
chl.emms.reef <- emmeans(chl.model, pairwise ~ origin_dest, weights = "proportional", adjust="none")
# P.value adjustment of the Bonferroni
rbind(chl.emms.reef$contrasts, adjust="bonferroni")

# comparing between reef sites within each time treatment
chl.emms.reef <- emmeans(chl.model, pairwise ~ time, weights = "proportional", adjust="none")
# P.value adjustment of the Bonferroni
rbind(chl.emms.reef$contrasts, adjust="none")

# Mixed Effects model origin_dest - random effect of colony nested in tank
# Since the model is completely balanced the SS type I, II or III will provide the same result
chl.model2 <- lmer(totchl ~ origin_dest * time *trt + (1|tank/colony), data = dat)
summary(chl.model2)
anova(chl.model2)
rand(chl.model2)
confint(chl.model2,oldNames=F)

# Model fitting and assumptions diagnostic 
x = residuals(chl.model2)
plot(fitted(chl.model2), x) 
#leveneTest(x ~ grtrt * time * tank* colony, data=dat, center=mean) # formal statistical test for homogeinity of variance (not recommended due the small sample size)
plot(chl.model2) # Residual vs Fitted values
qqnorm(x); qqline(x) # qq plot to check for normal distribution of residuals
hist(x) # histogram of residuals to check for normal distribution of residuals
shapiro.test(x) # formal statistical test (not recommended due the small sample size)
for(i in c("Aug-2016","Feb-2017")){
  print(paste0("Time",i))
  print(bartlett.test(totchl ~ origin_dest, data=dat[dat$time==i,]))
}

# comparing between reef sites within each time treatment
chl.emms.reef <- emmeans(chl.model2, pairwise ~ origin_dest|time|trt, weights = "proportional", adjust="none")
# P.value adjustment of the Bonferroni
rbind(chl.emms.reef$contrasts, adjust="bonferroni")

# comparing between reef sites within each time treatment
names(dat)
chl.emms.reef <- emmeans(chl.model2, pairwise ~time, weights = "proportional", adjust="none")
# P.value adjustment of the Bonferroni
rbind(chl.emms.reef$contrasts, adjust="none")

####CHLplots####
#chlrat#
sum<-Summarize(chlrat~origin_dest+time, data=heat, digits=3)
sum$se = sum$sd/sqrt(sum$n)
sum$CI_lower2 <- sum$mean - qnorm(0.975)*sum$se
sum$CI_upper2 <- sum$mean + qnorm(0.975)*sum$se
sum$Origin<-factor(c("HV","MV","LV","HV","MV","LV"))
sum$dest<-factor(c("HV","HV","HV","MV","MV","LV"))
sum$Time=factor(c("Aug 2016","Aug 2016","Aug 2016","Aug 2016","Aug 2016","Aug 2016","Feb 2017","Feb 2017","Feb 2017","Feb 2017","Feb 2017", "Feb 2017"))
sum$Origin=factor(sum$Origin,levels=unique(sum$Origin))
sum$dest=factor(sum$dest,levels=unique(sum$dest))

#points with sd bars
pd=position_dodge(.75)
plot.chlrat<-ggplot(sum,aes(x=dest,y=mean,color=Origin)) + 
  geom_point(aes(shape=Origin),size=3.5,position=pd) + 
  geom_errorbar(aes(ymin=mean - sd,ymax=mean + sd), width=.2,position=pd) +
  facet_wrap(~Time) +
  theme_bw() + 																		
  theme(plot.margin=margin(0.5,0.5,0.5,0.5,'cm'), panel.grid.minor=element_blank()) + 	
  scale_colour_manual(values=c("black","black","black")) + 						
  coord_cartesian(xlim=c(0.5,3.25), ylim=c(0,1.5), expand=F) + 
  ylab("Total Chlorophyll (ug/cm-2) Retained") + 
  xlab("Transplant Site")
save_plot("./chl/plots/20200527_por2_chlratio-sd-black.pdf", plot.chlrat,
          base_aspect_ratio = 1.3)

#end of ramp total chl points with sd or ci 95% bars, comment out which you want
sum=Summarize(totchl~grtrt+time, data=dat, digits=3)
sum$se = sum$sd/sqrt(sum$n)
sum$CI_lower2 <- sum$mean - qnorm(0.975)*sum$se
sum$CI_upper2 <- sum$mean + qnorm(0.975)*sum$se
sum$time=factor(sum$time,levels=unique(sum$time))
sum$Origin<-factor(c("HV","MV","LV","HV","MV","LV","HV","MV","LV","HV","MV","LV","HV","MV","LV","HV","MV","LV","HV","MV","LV","HV","MV","LV"))
sum$dest<-factor(c("HV","HV","HV","MV","MV","LV","HV","HV","HV","MV","MV","LV","HV","HV","HV","MV","MV","LV","HV","HV","HV","MV","MV","LV"))
sum$Origin=factor(sum$Origin,levels=unique(sum$Origin))
sum$dest=factor(sum$dest,levels=unique(sum$dest))
sum$Treatment=factor(c("cont","cont","cont","cont","cont","cont","heat","heat","heat","heat","heat","heat","cont","cont","cont","cont","cont","cont","heat","heat","heat","heat","heat","heat"))
sum

pd=position_dodge(.75)
plot.totalchl<-ggplot(sum,aes(x=dest,y=mean,fill=Origin,color=Treatment)) + 
  geom_point(aes(shape=Origin),size=3,position=pd) +
  #geom_errorbar(aes(ymin=mean - sd,ymax=mean + sd), width=0.2,position=pd) + 
  geom_errorbar(aes(ymin=CI_lower2, ymax=CI_upper2, color=Treatment), width=.2, position=pd) + 
  facet_wrap(~time) +
  theme_bw() + 
  theme(plot.margin=margin(0.5,0.5,0.5,0.5,'cm'), panel.grid.minor=element_blank()) +
  scale_colour_manual(values=c("gray","black")) +
  coord_cartesian(xlim=c(0.5,3.4), ylim=c(0,26), expand=F) + 
  ylab("Total Chlorophyll (ug/cm-2)\n after Acute Heat Stress") + 
  xlab("Transplant Site")
save_plot("./chl/plots/20210329_porII_correctedtotalchl-trt-ci.pdf", plot.totalchl,
          base_aspect_ratio = 1.3)

#boxplot#
chl <- ggplot(data=heat, 
              aes(x=dest, y=chlrat, label= time, fill=origin)) +
  scale_fill_manual(values = c ("red", "gold", "blue"), name = "Reef site") +
  stat_boxplot(geom ='errorbar', width = 0.5, lwd=0.5)+
  geom_boxplot(width=0.5, lwd=0.5, fatten=1) +
  expand_limits(y = 0)+
  facet_grid(~time, space = "free", scales = "free")+ #this can also be used but rows and columns have to be specified -> facet_grid(. ~ experiment)
  theme_bw() +
  #annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, colour = "black", size=2) +
  theme(line= element_line(size = 0.75),
        axis.line = element_line(colour = "grey20"),
        axis.ticks.length = unit(0.2 , "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(), 
        strip.text.x = element_text(color = "grey20", size = 16, angle = 0, hjust = 0, vjust = 0.5, face = "plain"),
        panel.spacing = unit(2, "lines"))
chl + xlab(label = "Transplant Site") + ylab(label = "Proportion of Total Chlorophyll Retained (ug cm-2)")+
  theme(axis.text.x = element_text(color = "grey20", size = 13, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 13, angle = 0, hjust = .5, vjust = .5, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 14, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 14, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        legend.title = element_text(colour="grey20", size=14, face="bold"),
        legend.text = element_text(colour="grey20", size=13, face="plain"),
        legend.position="bottom")
ggsave("./chl/plots/20210330_chlratio_boxplot.pdf", width = 10, height = 6)


####PAM&CHLcombo####
#Fv/Fm against total chl scatter
names(heat)
ggplot() +
  geom_point(heat,mapping=aes(x=chlrat,y=mqy,colour=origin_dest,shape=time),size=5) +
  theme_bw() + 
  theme(plot.margin=margin(0.5,0.5,0.5,0.5,'cm'),
        panel.grid.minor = element_blank()) +
  geom_abline(slope=1,intercept=0) +
  scale_colour_manual(values=c("red","gold","blue","red","gold","blue")) + 	
  ylab("Fv/Fm retention") + 
  xlab("Total chl retention") +
  theme(legend.position="bottom", legend.text = element_text(size=8), legend.title=element_text(size=9)) +
  labs(shape = "Time",color = "Origin_Dest") +
  guides(shape = guide_legend(override.aes = list(size = 3)),color = guide_legend(override.aes = list(size = 3)))

cor.test(heat$mqy,heat$chlrat,
         use = "complete.obs",
         method="pearson",
         adjust="none",     # Can adjust p-values; see ?p.adjust for options
         alpha=.05)
scatter<-ggscatter(heat, x = "mqy", y = "chlrat", 
                   add = "reg.line", conf.int = TRUE, 
                   cor.coef = TRUE, cor.method = "pearson",
                   xlab = "Fv/Fm Retained", ylab = "Total Chlorophyll Retained")
ggsave("./20210330_cortest-mqyVSchlrat.pdf", width = 10, height = 6)


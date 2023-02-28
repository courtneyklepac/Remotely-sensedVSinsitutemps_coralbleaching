library(Hmisc)
library(MASS)
library(ggplot2)
library(lme4)
library(FSA)
library(qqplotr)
library(rcompanion)
require(multcomp)
require(nlme)
library(dplyr)
library(lmerTest)
library(psych)
library(PerformanceAnalytics)
library(cowplot)

#read in file
growth<-read.csv("AS_por2.csv", header=T)
#set factor levels
growth$origin_dest=as.factor(growth$origin_dest)
growth$origin=as.factor(growth$origin)
growth$dest=as.factor(growth$dest)
growth$time=as.factor(growth$time)
growth$tank=as.factor(growth$tank)

growth$origin_dest=factor(growth$origin_dest,levels(growth$origin_dest)[c(1,5,3,2,6,4)])
growth$dest=factor(growth$dest,levels(growth$dest)[c(1,3,2)])
growth$origin=factor(growth$origin,levels(growth$origin)[c(1,3,2)])
growth$time=factor(growth$time,levels(growth$time)[c(1,3,2)])
levels(growth$time)[levels(growth$time)=="1"] <- "Aug 2016"
levels(growth$time)[levels(growth$time)=="6"] <- "Feb 2017"
levels(growth$time)[levels(growth$time)=="24"] <- "Jun 2018"

#removing 24mo from df
growth<-subset(growth, time !='Jun 2018')
growth$time<-factor(growth$time)
#removing non-assayed samples
growth<-subset(growth, trt !="NA")
growth$trt<-factor(growth$trt)

#creating interactions for post hoc comparisons
growth$grtime<-interaction(growth$time,growth$origin_dest) 
growth$oritime<-interaction(growth$origin,growth$time)
growth$desttime<-interaction(growth$dest,growth$time)

# Statistical tests for growth rate (norm to weeks)
# Analysis of variance (2-way) ANOVA pam.growtha (mixed linear model)
# Mixed Effects model origin_dest - random effect of colony
# Since the model is completely balanced the SS type I, II or III will provide the same result
gr.model2 <- lmer(wkgrate ~ origin_dest * time  + (1|colony) , growtha = growth)
summary(gr.model2)
anova(gr.model2)
rand(gr.model2)
confint(gr.model2,oldNames=F)

# Model fitting and assumptions diagnostic 
x = residuals(gr.model2)
#leveneTest(x ~ grtrt * time * tank* colony, growtha=growth, center=mean) # formal statistical test for homogeinity of variance (not recommended due the small sample size)
plot(gr.model2) # Residual vs Fitted values
qqnorm(x); qqline(x) # qq plot to check for normal distribution of residuals
hist(x) # histogram of residuals to check for normal distribution of residuals
shapiro.test(x) # formal statistical test (not recommended due the small sample size)

# comparing between reef sites within each time treatment
gr.emms.reef <- emmeans(gr.model2, pairwise ~ origin_dest, weights = "proportional", adjust="none")
# P.value adjustment of the Bonferroni
rbind(gr.emms.reef$contrasts, adjust="bonferroni")

# summary descriptive stats for table in manuscript
sum<-Summarize(wkgrate~origin_dest+time, data=growth, digits=3)
sum$se = sum$sd/sqrt(sum$n)
sum$CI_lower2 <- sum$mean - qnorm(0.975)*sum$se
sum$CI_upper2 <- sum$mean + qnorm(0.975)*sum$se
#sum$time=factor(sum$time,levels= c("1mo", "6mo", "24mo"))
sum$Origin<-factor(c("HV","MV","LV","HV","MV","LV","HV","MV","LV","HV","MV","LV"))
sum$dest<-factor(c("HV","HV","HV","MV","MV","LV","HV","HV","HV","MV","MV","LV"))
sum$Origin=factor(sum$Origin,levels=unique(sum$Origin))
sum$Time=factor(c("Aug 2016","Aug 2016","Aug 2016","Aug 2016","Aug 2016","Aug 2016","Feb 2017","Feb 2017","Feb 2017","Feb 2017","Feb 2017","Feb 2017"))
sum$dest=factor(sum$dest,levels=unique(sum$dest))

#plot with CI error bars
pd=position_dodge(.75)
plot<-ggplot(sum,aes(x=dest,y=mean,colour=Origin)) + 
  geom_point(aes(shape=Origin),size=3,position=pd) +
  geom_errorbar(aes(ymin=CI_lower2, ymax=CI_upper2, color=Origin), width=.2, position=pd) + 
  facet_wrap(~time) +
  theme_bw() + 
  theme(plot.margin=margin(0.5,0.5,0.5,0.5,'cm'), panel.grid.minor=element_blank()) +
  scale_colour_manual(values=c("red","gold","blue")) +
  coord_cartesian(xlim=c(0.5,3.4), ylim=c(0,0.4), expand=F) + 
  ylab("Weekly Growth (g/wk)") + 
  xlab("Transplant Site")
save_plot("./growthsurvival/plots/202010_porII_wkgrate_ci.pdf", plot,
          base_aspect_ratio = 1.3)
          

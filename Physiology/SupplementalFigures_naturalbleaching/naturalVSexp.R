#######################
# NAT V EXP BLEACHING #

library(lme4)
library(lmerTest)
library(stats)
library(Hmisc)
library(ggpubr)
library(FSA)

##Read in data
data <- read.csv("natVexp-donoravg.csv",header=T)
summary(data)
data$originsample=as.factor(data$originsample)
data$sample=as.factor(data$sample)
data$site=as.factor(data$site)
data$colony=as.factor(data$colony)
data$originsample=as.factor(data$originsample)
levels(data$sample)[levels(data$sample)=="light"] <- "donor"
data$site=factor(data$site,levels(data$site)[c(1,3,2)])
#data$sample=factor(data$sample,levels(data$sample)[c(1,3,2)])

sum<-Summarize(perbleached~site, data=data, digits=3)

####CHL aov####
chl <- lmer(totalchl ~ site* sample  + (1|colony), data = data)
anova(chl)
# comparing among sample types
chl.emms.reef <- emmeans(chl, pairwise ~ sample, weights = "proportional", adjust="none")
# P.value adjustment of the Bonferroni
rbind(chl.emms.reef$contrasts, adjust="bonferroni")

####a priori contrasts####
levels(data$sample) # cont heat light
#standard contrast of sample
contrasts(data$sample)
model <- lmer(totalchl~sample +(1|colony), data = data)
summary(model)
#planned orthogonal contrasts
#contrast exp 0, 0, 1, -1
#contrast nat 1, -1, 0, 0
#contrast donor baseline diff than exp 0, 1, 0, -1
#contrast heat stress -1, 0, 1, 0

contrastmatrix<-cbind( c(0, 0, 1, -1),
                       c(1, -1, 0, 0),
                       c(-1, 0, 1, 0),
                       c(0, -1, 0, 1)
                        )
contrastmatrix
#contrasts(data$sample)<-contrastmatrix
data$sample
#aov w contrasts
contrast_mod<-aov(totalchl ~ sample + Error(colony), data)
summary(contrast_mod, split = list(sample = list("Exp" = 1, "Nat" = 2, "Heat" = 3,"Baseline" = 4))) 
#lmer w contrasts
contrast_lme<-lmer(totalchl ~ sample + (1|colony), data)
summary(contrast_lme) 

####CHL plot####
#Totalchl across sample type boxplot
bleach <- ggplot(data=data, 
                aes(x=sample, y=totalchl, fill=site)) +
  scale_fill_manual(values = c ("red", "gold", "blue"), name = "") +
  stat_boxplot(geom ='errorbar', width = 0.7, lwd=0.5)+
  geom_boxplot(width=0.7, lwd=0.5, fatten=1) +
  expand_limits(y = 0)+
  theme_bw() +
xlab(label = "") + ylab(label = "Total Chlorophyll (µg cm-2)")+
  theme(axis.text.x = element_text(color = "grey20", size = 13, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 13, angle = 0, hjust = .5, vjust = .5, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 14, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 14, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        legend.title = element_text(colour="grey20", size=14, face="bold"),
        legend.text = element_text(colour="grey20", size=13, face="plain"),
        legend.position="bottom")

ggsave("./natVexp-boxplot-20210520.pdf", width = 10, height = 6)

#% bleaching boxplot
sum<-Summarize(totalchl~sample+site, data=data, digits=3)
sum$se = sum$sd/sqrt(sum$n)
sum$CI_lower2 <- sum$mean - qnorm(0.975)*sum$se
sum$CI_upper2 <- sum$mean + qnorm(0.975)*sum$se
sum$site=factor(sum$site,levels=unique(sum$site))

pd=position_dodge(.75)
plot.totalchl<-ggplot(sum,aes(x=sample,y=mean,colour=site)) + 
  geom_point(aes(shape=site),size=3,position=pd) +
  geom_errorbar(aes(ymin=CI_lower2, ymax=CI_upper2), width=.2, position=pd) + 
  theme_bw() + 
  theme(plot.margin=margin(0.5,0.5,0.5,0.5,'cm'), panel.grid.minor=element_blank()) +
  scale_colour_manual(values=c("red","gold","blue")) +
  coord_cartesian(xlim=c(0.5,3.4), ylim=c(0,20), expand=F) + 
  ylab("Total Chlorophyll (µg cm-2)") + 
  xlab("Sample Type")
ggsave("./20210520_sampletype-totalchl-ci.pdf", width = 10, height = 6)


bleach<-ggplot(data=data, 
       aes(x=site, y=perbleached, fill=site)) +
  scale_fill_manual(values = c ("red", "gold", "blue")) +
  stat_boxplot(geom ='errorbar', width = 0.7, lwd=0.5)+
  geom_boxplot(width=0.7, lwd=0.5, fatten=1) +
  expand_limits(y = 0)+
  theme_bw() +
  xlab(label = "Reef Site") + ylab(label = "% Bleaching")+
  theme(axis.text.x = element_text(color = "grey20", size = 13, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 13, angle = 0, hjust = .5, vjust = .5, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 14, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 14, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        legend.title = element_text(colour="grey20", size=14, face="bold"),
        legend.text = element_text(colour="grey20", size=13, face="plain"),
        legend.position="bottom")

ggsave("./perbleach-boxplot-202105132.pdf", width = 10, height = 6)


#####correlation#####
#subset new df based on site
HV<-data[data$origin=='HV',]
MV<-data[data$origin=='MV',]
LV<-data[data$origin=='LV',]
cont<-data[data$sample=='cont',]
donor<-data[data$sample=='donor',]
heat<-data[data$sample=='heat',]
HV.cont<-HV[HV$trt=='cont',]
HV.heat<-HV[HV$trt=='heat',]

total <- merge(donor,cont,by=c("site","colony"))
total$site=factor(total$site,levels(total$site)[c(1,3,2)])
summary(total2)
names(total)
names(total)[8] <- "donorchl"
names(total)[17] <- "controlchl"
total2 <- merge(total,heat,by=c("site","colony"))
total2<-subset(total2, select = c(site,colony,donorchl,controlchl,totalchl,days31,days32,inDHW) )

shapiro.test(total2$donorchl)
shapiro.test(total2$controlchl)
shapiro.test(total2$totalchl)#heat

#individual Pearson correlations, change to variables of interest
cor.test(heat$totalchl,heat$inDHW,
         use = "complete.obs",
         method="pearson",
         adjust="none",     # Can adjust p-values; see ?p.adjust for options
         alpha=.05)
scatter<-ggscatter(heat, x = "inDHW", y = "totalchl", 
                   add = "reg.line", conf.int = TRUE, 
                   cor.coef = TRUE, cor.method = "pearson",
                   xlab = "in situ DHW", ylab = "Heat Ramet Total Chlorophyll (µg cm-2)")
ggsave("./20210513_heatchl-inDHW_corrscatter.pdf", width = 10, height = 6)


###comparing donor and heat chl##
cor.test(total$donorchl,total$heatchl,
         use = "complete.obs",
         method="pearson",
         adjust="none",     # Can adjust p-values; see ?p.adjust for options
         alpha=.05)

scatter<-ggscatter(total, x = "donorchl", y = "heatchl", 
                   add = "reg.line", conf.int = TRUE, 
                   cor.coef = TRUE, cor.method = "pearson",
                   xlab = "Donor Colony Total Chlorophyll (µg cm-2)", ylab = "Heat Ramet Total Chlorophyll (µg cm-2)")
ggsave("./cortest-natVexp_20210513.pdf", width = 10, height = 6)

###comparing donor and control chl##
cor.test(total$donorchl,total$controlchl,
         use = "complete.obs",
         method="pearson",
         adjust="none",     # Can adjust p-values; see ?p.adjust for options
         alpha=.05)

scatter<-ggscatter(total, x = "donorchl", y = "controlchl", 
                   add = "reg.line", conf.int = TRUE, 
                   cor.coef = TRUE, cor.method = "pearson",
                   xlab = "Donor Colony Total Chlorophyll (µg cm-2)", ylab = "Control Ramet Total Chlorophyll (µg cm-2)")
ggsave("./cortest-contnatVexp_20210513.pdf", width = 10, height = 6)

#correlation matrix of chl and environment
C<-cor(total2[,c(3:8)],use = "pairwise.complete.obs")
round(C,2)
C
# mat : is a matrix of data
# ... : further arguments to pass to the native R cor.test function
cor.mtest <- function(C, ...) {
  mat <- as.matrix(C)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}
# matrix of the p-value of the correlation
p.mat <- cor.mtest(total2[,c(3:8)])
p.mat
pdf(file = "20210520_natbleaching_corrplot.pdf")
corrplot(C, method="color",type="upper", col=brewer.pal(n=8, name="RdYlBu"), addCoef.col = "black",diag=FALSE, p.mat=p.mat,sig.level = 0.05,tl.col="black")
dev.off()

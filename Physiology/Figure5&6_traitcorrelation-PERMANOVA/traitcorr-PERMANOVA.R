#correlations between traits and temp metrics

#libraries
library(car)
library(vegan)
library(lmerTest)
library(emmeans)
library(ggplot2)
library(ggpubr) 
library(Hmisc)
library(dplyr)
library(FSA)
library(FactoMineR)
require(factoextra)
require(tidyr)
require(dplyr)
require(MASS)
require(reshape2)
require(cowplot)
library(car)
library(corrplot)
library(factoextra)
library(RColorBrewer)

####data input and organize####
#ratios
corr.data <- read.csv("AS_por2-2.csv", header = T)
#'raw'*use this file
corr.data <- read.csv("AS_por2-num.csv", header = T, na.strings=c("","NA"))

colnames(corr.data)
summary(corr.data)

#set factor levels
corr.data$origin <- as.factor(corr.data$origin) 
corr.data$origin_dest <- as.factor(corr.data$origin_dest) 
corr.data$dest <- as.factor(corr.data$dest)
corr.data$colony <- as.factor(corr.data$colony)

corr.data$origin_dest=factor(corr.data$origin_dest,levels(corr.data$origin_dest)[c(1,5,3,2,6,4)])
corr.data$dest=factor(corr.data$dest,levels(corr.data$dest)[c(1,3,2)])
corr.data$origin=factor(corr.data$origin,levels(corr.data$origin)[c(1,3,2)])
corr.data$trt <- as.factor(corr.data$trt) #'0' is control and '1' is heat
corr.data$time <- as.factor(corr.data$time) #number of months

head(corr.data)

#removing 24mo and trt NA's from df
corr.data<-subset(corr.data, time !='24')
corr.data$time<-factor(corr.data$time)
corr.data<-subset(corr.data, trt !="NA")
corr.data$trt<-factor(corr.data$trt)

#calculate max quant yield, Fv/Fm retained after heat stress
corr.data$mqy<-(1-corr.data$ylossnorm)
#subset into trts
corr.heat<-corr.data[corr.data$trt=='1',]
corr.cont<-corr.data[corr.data$trt=='0',]
#calculated total chl retained after heat stress
corr.heat$chlrat<-(1-(corr.cont$totchl-corr.heat$totchl)/corr.cont$totchl)
summary(corr.heat)
corr.heat$mqyrat<-(1-(corr.cont$mqy-corr.heat$mqy)/corr.cont$mqy)

####stepwise AIC####
#which temp metrics are best to use in overall model
corr.data.num=dplyr::select(corr.data,origin_dest,time,trt,wkgrate,mqy,totalchl,trt,mean,min,mmm,dtr,maxdtr,dr90,days31,days32)

sum(is.na(corr.data.num))
corr.data.num=na.omit(corr.data.num)
multreg.gr<- lm(wkgrate ~origin_dest+time+trt+mean+min+mmm+dtr+maxdtr+dr90+days31+days32,data=corr.data.num)
stepAIC(multreg.gr, direction='both')

multreg.mqy<- lm(mqy~origin_dest+time+trt+mean+min+mmm+dtr+maxdtr+dr90+days31+days32,data=corr.data.num)
stepAIC(multreg.mqy, direction='both')

multreg.chl<- lm(totalchl ~origin_dest+time+trt+mean+min+mmm+dtr+maxdtr+dr90+days31+days32,data=corr.data.num)
stepAIC(multreg.chl, direction='both')

####correlations#####
##entire trait data (all, heat, cont) and temp metrics for all sites##
names(corr.data)
#PCA include trt
corr.data.num=dplyr::select(corr.data,time,trt,wkgrate,ylossnorm,totchl,trt,mean,min,mmm,dtr,dr90,inDHW,nightDHW,CRWDHW,inDHM)
summary(corr.data.num)
#set trt (cont=0, heat=1) and time to numeric values
str(corr.data.num)

corr.data.num[,'trt'] <- as.numeric(as.character(corr.data.num[,'trt']))
corr.data.num[,'time'] <- as.numeric(as.character(corr.data.num[,'time']))

C<-cor(corr.data.num[,c(3:14)],use = "pairwise.complete.obs")
round(C,2)

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
p.mat <- cor.mtest(corr.data.num[,c(3:14)])
p.mat
pdf(file = "20200622_all_traittemps_corrplot.pdf")
corrplot(C, method="color",type="upper", col=brewer.pal(n=8, name="RdYlBu"), addCoef.col = "black",diag=FALSE, p.mat=p.mat,sig.level = 0.05,tl.col="black")
dev.off()

#scatterplot matrix
chart.Correlation(corr.data.num[,c(3:14)], histogram=TRUE, pch=19)
colnames(corr.data.num)


##looking w/in individual treats- control most similar to field transplants##
#corr.heat.num=dplyr::select(corr.heat,time,wkgrate,mqy,totalchl,mean,min,mmm,dtr,maxdtr,dr90,days31,days32)
corr.cont.num=dplyr::select(corr.cont,wkgrate,mqy,totchl,mean,min,mmm,dtr,dr90,CRWDHW,nightDHW,inDHW,inDHM)
str(corr.cont.num)
corr.cont.num[,'time'] <- as.numeric(as.character(corr.cont.num[,'time']))

C<-cor(corr.cont.num,use = "pairwise.complete.obs")
round(C,2)

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
p.mat <- cor.mtest(corr.cont.num)
head(p.mat[, 1:5])
pdf(file = "20210411_cont_traittempsDHW_corrplot.pdf")
corrplot(C, method="color",type="upper", col=brewer.pal(n=8, name="RdYlBu"), addCoef.col = "black",diag=FALSE, p.mat=p.mat,sig.level = 0.05,tl.col="black")
dev.off()

#scatterplot matrix
chart.Correlation(corr.data.num, histogram=TRUE, pch=19)
colnames(corr.data.num)

####PCA####
#create log transformed dataset
#first create dataframe of elements I want
log.data<-(corr.cont.num+1)
names(log.data)
dim(log.data)#120 13
#to add origin_dest from original to df, but need to remove NA rows to match dim with log.data
dim(corr.cont.num)
#which(is.na(corr.cont.num$totalchl)) #23, 63, 108
#which(is.na(corr.cont.num$ylossnorm)) #27
#corr.data2 <- corr.data[-c(23, 27, 63, 108), ]
#take natural log of numerical values and append factor columns to dataframe
nl<-log(log.data)
dim(nl)
nl[is.nan(nl)] <- 0

nl<-cbind(nl,corr.cont$origin_dest,corr.cont$colony,corr.cont$time)
names(nl)[names(nl) == "corr.cont$origin_dest"] <- "origin_dest"
names(nl)[names(nl) == "corr.cont$time"] <- "time"
names(nl)[names(nl) == "corr.cont$colony"] <- "colony"

res.pca <- prcomp(nl[,c(1:3)], center = TRUE,scale. = TRUE)
fviz_eig(res.pca)
eig.val <- get_eigenvalue(res.pca)
eig.val
var <- get_pca_var(res.pca)
var
var$contrib
var$cor
fviz_pca_var(res.pca, col.var = "black")
# Total cos2 of variables on Dim.1 and Dim.2
fviz_cos2(res.pca, choice = "var", axes = 1:2)
# Color by contribution: quality on the factor map
fviz_pca_var(res.pca, col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07")
)
# Contributions of variables to PC1
fviz_contrib(res.pca, choice = "var", axes = 1, top = 10)
# Contributions of variables to PC2
fviz_contrib(res.pca, choice = "var", axes = 2, top = 10)

#this plot combines ind and variables into one plot
cols <- c("HV_HV" = "red", "MV_HV" = "gold","LV_HV" = "blue", "HV_MV" = "darkgoldenrod","MV_MV"="tomato","LV_LV"="dodgerblue")

pcabiplot<-fviz_pca_biplot(res.pca,
             geom.ind = "point", mean.point = FALSE,# show points only (but not "text")
             habillage = nl$origin_dest, # color by groups
             pointsize = 3,
             palette = cols,
             addEllipses = FALSE, # ellipses
             label = "var",
             col.var = "black", repel = TRUE,
             legend.title = "Groups") + theme_minimal() +
             scale_shape_manual(values=c(16,17,15,16,17,15))

#export
ggexport(pcabiplot, filename = "20210414_PCA-shapes_noenviro.pdf")

#predict with qualitative variables
time <- as.factor(corr.cont$time)
pqual<-fviz_pca_ind(res.pca,
             col.ind = time, # color by groups
             palette = c("#00AFBB",  "#FC4E07"),
             addEllipses = TRUE, # Concentration ellipses
             ellipse.type = "confidence",
             legend.title = "Time",
             repel = TRUE)
ggexport(pqual, filename = "20210414_PCA-time_noenviro.pdf")

#predict with quantitative variables
quanti.sup <- corr.cont.num[, 4:12, drop = FALSE]
head(quanti.sup)
quanti.sup$time<-as.numeric(quanti.sup$time)
# Predict coordinates and compute cos2
quanti.coord <- cor(quanti.sup, res.pca$x)
quanti.cos2 <- quanti.coord^2
# Graph of variables including supplementary variables
p <- fviz_pca_var(res.pca)
pfiz<-fviz_add(p, quanti.coord, color ="blue", geom="arrow")
ggexport(pfiz, filename = "20210414_qualpredictors-wenviroPCA.pdf")

####PERMANOVA####
library(devtools)
install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis",force=TRUE)
library(pairwiseAdonis)
corr.cont.num=dplyr::select(corr.cont,origin_dest,time,colony,wkgrate,mqy,totchl,mean,min,mmm,dtr,dr90,CRWDHW,nightDHW,inDHW,inDHM)
corr.cont.num$grtime<-interaction(corr.cont.num$origin_dest,corr.cont.num$time)
corr.cont.num.complete<-na.omit(corr.cont.num)
corr.cont.scale<-corr.cont.num.complete %>% mutate_if(is.numeric, scale)
# Creating separate response (com) and predictor (meta) variable files and compute distance matrices for each
corr.cont.com<-dplyr::select(corr.cont.scale,wkgrate,mqy,totchl)
corr.cont.meta<-dplyr::select(corr.cont.scale,origin_dest,time,grtime,colony,mean,min,mmm,dtr,dr90,CRWDHW,nightDHW,inDHW,inDHM)
# compute distance matrix
corr.cont.dist<- vegdist(corr.cont.com, method="euclidean",na.rm = TRUE)
# permanova
adonis(corr.cont.dist ~ origin_dest*time, data = corr.cont.meta, strata=corr.cont.meta$colony,na.rm = TRUE,method="eu")
##start copy here for function pairwise.adonis() just have to run once
pairwise.adonis <- function(x,factors, sim.function = 'vegdist', sim.method = 'bray', p.adjust.m ='holm')
{
  library(vegan)
  
  co = combn(unique(as.character(factors)),2)
  pairs = c()
  F.Model =c()
  R2 = c()
  p.value = c()
  
  
  for(elem in 1:ncol(co)){
    if(sim.function == 'daisy'){
      library(cluster); x1 = daisy(x[factors %in% c(co[1,elem],co[2,elem]),],metric=sim.method)
    } else{x1 = vegdist(x[factors %in% c(co[1,elem],co[2,elem]),],method=sim.method)}
    
    ad = adonis(x1 ~ factors[factors %in% c(co[1,elem],co[2,elem])] );
    pairs = c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    F.Model =c(F.Model,ad$aov.tab[1,4]);
    R2 = c(R2,ad$aov.tab[1,5]);
    p.value = c(p.value,ad$aov.tab[1,6])
  }
  p.adjusted = p.adjust(p.value,method=p.adjust.m)
  sig = c(rep('',length(p.adjusted)))
  sig[p.adjusted <= 0.05] <-'.'
  sig[p.adjusted <= 0.01] <-'*'
  sig[p.adjusted <= 0.001] <-'**'
  sig[p.adjusted <= 0.0001] <-'***'
  
  pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted,sig)
  print("Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1")
  return(pairw.res)
  
} 
## end copy here
pairwise.adonis(corr.cont.com,factors=corr.cont.meta$time,sim.function='vegdist',sim.method='euclidean',p.adjust.m='bonferroni')
pairwise.adonis(corr.cont.com,factors=corr.cont.meta$grtime,sim.function='vegdist',sim.method='euclidean',p.adjust.m='bonferroni')

####HVMV-RTE####
corr.data <- read.delim("HVMVcorrtest.txt", header = T)
colnames(corr.data)

M<-cor(corr.data,use = "complete.obs")
head(round(M,2))

corrplot(M, method="color",type="upper", col=brewer.pal(n=8, name="YlOrRd"),sig.level = 0.05,tl.col="black")
# mat : is a matrix of data
# ... : further arguments to pass to the native R cor.test function
cor.mtest <- function(M, ...) {
  mat <- as.matrix(M)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], method="pearson")
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}
# matrix of the p-value of the correlation
p.mat <- cor.mtest(corr.data)
head(p.mat[, 1:5])
pdf(file = "20191220_HVMV_corrplot.pdf")
corrplot(M, method="color",type="upper", col=brewer.pal(n=8, name="YlOrRd"), addCoef.col = "black",diag=FALSE, p.mat=p.mat,sig.level = 0.05,tl.col="black")
dev.off()
#scatterplot matrix
chart.Correlation(corr.data, histogram=TRUE, pch=19)

####extra####
#additional graph option, makes PCA but variable vectors are separate figure and ellipses are split over time
# extract pc scores for first two component and add to dat dataframe
nl$pc1 <- res.pca$ind$coord[, 1] # indexing the first column

nl$pc2 <- res.pca$ind$coord[, 2]  # indexing the second column
#extract the data for the variable contributions to each of the pc axes
pca.vars <- res.pca$var$coord %>% data.frame
pca.vars$vars <- rownames(pca.vars)
pca.vars.m <- melt(pca.vars, id.vars = "vars")
#variable contribution plot has a circle around the variables that has a radius of 1
circleFun <- function(center = c(0,0),diameter = 1, npoints = 100){
  r = diameter / 2
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}

circ <- circleFun(c(0,0),2,npoints = 500)

nl<- transform(nl, time = as.factor(time))
str(nl)
p <- ggplot(data = nl, aes(x = pc1, y = pc2, colour = origin_dest, shape = time)) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_point(alpha = 0.8,size=2) + 
  scale_colour_manual(values = c("red", "gold", "blue", "darkgoldenrod","tomato", "dodgerblue")) 
p <- p + stat_ellipse(geom="polygon", aes(fill = origin_dest), 
                      alpha = 0.2, 
                      show.legend = FALSE, 
                      level = 0.95) +
  xlab("PC 1 (54.75%)") + 
  ylab("PC 2 (21.14%)") +
  theme_minimal() +
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill= "transparent"))

vars.p <-  ggplot() +
  geom_path(data = circ,aes(x,y), lty = 2, color = "white", alpha = 0.7) +
  geom_hline(yintercept = 0, lty = 2, color = "white", alpha = 0.9) +
  geom_vline(xintercept = 0, lty = 2, color = "white", alpha = 0.9) +  
  geom_segment(data = pca.vars, aes(x = 0, xend = Dim.1, y = 0, yend = Dim.2),
               arrow = arrow(length = unit(0.025, "npc"), type = "open"), 
               lwd = 1) + 
  geom_text(data = pca.vars, 
            aes(x = Dim.1*1.15, y =  Dim.2*1.15, 
                label = c("wkgrate", "mqy", "totalchl", "dtr", "min", "mean", "mmm")), 
            check_overlap = T, size = 3) +
  xlab("") + 
  ylab("") +
  coord_equal() +
  theme_minimal() +
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill= "transparent"))

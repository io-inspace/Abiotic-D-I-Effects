#####METADATA#####
#RESULTS 
#edited May 2021

#####call libraries and functions#####
########################################################################

library(lme4)
library(MuMIn)
library(ggplot2)
library(lsmeans)
library(car)
library(emmeans)
library(sjPlot)
library(DHARMa)
library(MASS)
library(ordinal)
library(sp)

diagnostics<- function(x,print=TRUE){
  plot(x) #using Cook's distance to check for outliers
  hist(x$residuals,breaks=3) #to check for normalcy
  a<-coef(x)#to get model coefficients
  b<-summary(x)#to get coefficients and stats
  c<-AIC(x)#self-explanatory
  d<-r.squaredGLMM(x)#r2 values for lme models only!  Use r2 stats from summary for lm models. 
  output<-list(a,b,c,d)
  return(output)
}
ctrl<-lmerControl(optimizer='bobyqa')#ctrl not implemented, just in case.
options(na.action=na.fail)

#####define habitat quality when species partition resources along a gradient#####
########################################################################

setwd("~/Documents/submissions/NichePart2020/for Dryad/landscapes_C")
bp<-read.csv('breakpoints.csv',header=FALSE,sep=',') #where species curves 
#intersect. i.e. defines the range of abiotic conditions for which 
#species x performs better than any other species.
#np<-read.csv('np.csv',header=FALSE,sep=',') #niche optima for each species

#####get abiotic predictors#####
########################################################################

smallpred<-read.csv('smallpred.csv',header = TRUE,sep=',')
largepred<-read.csv('largepred.csv',header = TRUE,sep = ',')

#####get S and comphet by simulation#####
########################################################################

models<-c('neutRdiffGneutC','neutRneutGdiffC')#'neutRneutGneutC','neutRdiffGneutC','neutRneutGdiffC','neutRdiffGdiffC')
ddstrength<-c('_0_5','_1_0','_1_5')
dispersal<-c('adj','int','uni')#'adj','int','uni')


#VIOLIN PLOTS:

#D:
d = 1
s = 3
m = 1

setwd(paste("~/Documents/submissions/NichePart2020/for Dryad/2020/code_editedsumm2019/neutral/",dispersal[d],'/',ddstrength[s],'/',models[m],sep=''))
#setwd(paste('~/Documents/MATLAB/simulation2019Jan/neutral/',dispersal[d],'/',ddstrength[s],'/',models[m],sep=''))
results_small<-read.csv(paste('results_small.csv',sep=''), header = TRUE, sep = ',')
results_large<-read.csv(paste('results_large.csv',sep=''), header = TRUE, sep = ',')

results_small<-cbind(smallpred,results_small)
results_large<-cbind(largepred,results_large)

#remove NaNs, round habitat
results_small<-results_small[!is.nan(results_small$comphet),]
results_small$meanC<-round(results_small$meanC,1)

#Violins:
#the A/S:
#you need to work from 2 data.frames, one with binned x and one with continous x:
binnedresults<-results_small
binnedresults$varC<-as.factor(round(binnedresults$varC,0))

#mod<-lmer(S~varC+(1|landscape/meanC),data=results_small,control=ctrl)
binnedmod<-lmer(S~varC+(1|landscape:meanC),data=binnedresults,control=ctrl)

DplotAS<-ggplot(binnedresults, aes(x=varC, y=S))+
  geom_violin(alpha=0,color=rgb(215/255,48/255,31/255),size=1.25)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black"),
            axis.text=element_text(size=12))+
            stat_summary(aes(group=1),fun=mean, geom="point", 
                         color="black", size=3)+ylim(0,4)+
  geom_abline(intercept=summary(binnedmod)$coefficients[1,1],slope=summary(binnedmod)$coefficients[2,1],size=2)

#points:

# DplotAS_points<-ggplot(results_small, aes(x=varC, y=S))+
#   geom_point(alpha=.25,color='blue')+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#                                            panel.background = element_blank(),axis.line = element_line(colour = "black"),axis.text=element_text(size=12))+
# ylim(0,4)+geom_abline(intercept=summary(mod)$coefficients[1,1],slope=summary(mod)$coefficients[2,1],size=3)

#the A/C:

mod<-lmer(comphet~varC+(1|landscape:meanC),data=binnedresults,control=ctrl)

DplotAC<-ggplot(binnedresults, aes(x=varC, y=comphet))+
  geom_violin(alpha=0,color=rgb(215/255,48/255,31/255),size=1.25)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              axis.line = element_line(colour = "black"),
              axis.text=element_text(size=12))+
  stat_summary(aes(group=1),fun=mean, geom="point", 
               color="black", size=3)+
  geom_abline(intercept=summary(mod)$coefficients[1,1],slope=summary(mod)$coefficients[2,1],size=2)

# DplotAC_points<-ggplot(results_small, aes(x=varC, y=comphet))+
#   geom_point(alpha=.25,color='blue')+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#                                            panel.background = element_blank(),axis.line = element_line(colour = "black"),axis.text=element_text(size=12))+
#  geom_abline(intercept=summary(mod)$coefficients[1,1],slope=summary(mod)$coefficients[2,1],size=3)


#the C/S:

binnedresults$comphet<-as.factor(round(binnedresults$comphet,1))

mod<-lmer(S~comphet+(1|landscape:meanC),data=binnedresults,control=ctrl)

DplotCS<-ggplot(binnedresults, aes(x=comphet, y=S))+
  geom_violin(alpha=0,color=rgb(215/255,48/255,31/255),size=1.25)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               axis.line = element_line(colour = "black"),
                axis.text=element_text(size=12))+
  stat_summary(aes(group=1),fun=mean, geom="point", 
               color="black", size=3)+
  ylim(0,4)+geom_abline(intercept=summary(mod)$coefficients[1,1],slope=summary(mod)$coefficients[2,1],size=2)

# DplotCS_points<-ggplot(results_small, aes(x=comphet, y=S))+
#   geom_point(alpha=.25,color='blue')+ylim(0,4)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#                                                      panel.background = element_blank(),axis.line = element_line(colour = "black"),axis.text=element_text(size=12))+
# ylim(0,4)+geom_abline(intercept=summary(mod)$coefficients[1,1],slope=summary(mod)$coefficients[2,1],size=3)
# 

#I
d = 1
s = 3
m = 2

setwd(paste("~/Documents/submissions/NichePart2020/for Dryad/2020/code_editedsumm2019/neutral/",dispersal[d],'/',ddstrength[s],'/',models[m],sep=''))
#setwd(paste('~/Documents/MATLAB/simulation2019Jan/neutral/',dispersal[d],'/',ddstrength[s],'/',models[m],sep=''))
results_small<-read.csv(paste('results_small.csv',sep=''), header = TRUE, sep = ',')
results_large<-read.csv(paste('results_large.csv',sep=''), header = TRUE, sep = ',')

results_small<-cbind(smallpred,results_small)
results_large<-cbind(largepred,results_large)

#remove NaNs, round habitat
results_small<-results_small[!is.nan(results_small$comphet),]
results_small$meanC<-round(results_small$meanC,1)

#the A/S:
#you need to work from 2 data.frames, one with binned x and one with continous x:
binnedresults<-results_small
binnedresults$varC<-as.factor(round(binnedresults$varC,0))

mod<-lmer(S~varC+(1|landscape:meanC),data=binnedresults,control=ctrl)

IplotAS<-ggplot(binnedresults, aes(x=varC, y=S))+
  geom_violin(alpha=0,color=rgb(253/255,204/255,138/255),size=1.25)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                    panel.background = element_blank(),
                                                    axis.line = element_line(colour = "black"),
                                                    axis.text=element_text(size=12))+
  stat_summary(aes(group=1),fun=mean, geom="point", 
               color="black", size=3)+ylim(0,4)+
  geom_abline(intercept=summary(mod)$coefficients[1,1],slope=summary(mod)$coefficients[2,1],size=2)

#points:
# 
# IplotAS_points<-ggplot(results_small, aes(x=varC, y=S))+
#   geom_point(alpha=.25,color='red')+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#                                            panel.background = element_blank(),axis.line = element_line(colour = "black"),axis.text=element_text(size=12))+
#   ylim(0,4)+geom_abline(intercept=summary(mod)$coefficients[1,1],slope=summary(mod)$coefficients[2,1],size=3)

#the A/C:
mod<-lmer(comphet~varC+(1|landscape:meanC),data=binnedresults,control=ctrl)

IplotAC<-ggplot(binnedresults, aes(x=varC, y=comphet))+
  geom_violin(alpha=0,color=rgb(253/255,204/255,138/255),size=1.25)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                    panel.background = element_blank(),
                                                    axis.line = element_line(colour = "black"),
                                                    axis.text=element_text(size=12))+
  stat_summary(aes(group=1),fun=mean, geom="point", 
               color="black", size=3)+ylim(0,.5)+
  geom_abline(intercept=summary(mod)$coefficients[1,1],slope=summary(mod)$coefficients[2,1],size=2)

# IplotAC_points<-ggplot(results_small, aes(x=varC, y=comphet))+
#   geom_point(alpha=.25,color='red')+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#                                            panel.background = element_blank(),axis.line = element_line(colour = "black"),axis.text=element_text(size=12))+
#   ylim(0,.5)+
#   geom_abline(intercept=summary(mod)$coefficients[1,1],slope=summary(mod)$coefficients[2,1],size=3)

#the C/S:
binnedresults$comphet<-as.factor(round(binnedresults$comphet,2))

mod<-lmer(S~comphet+(1|landscape:meanC),data=binnedresults,control=ctrl)

IplotCS<-ggplot(binnedresults, aes(x=comphet, y=S))+
  geom_violin(alpha=0,color=rgb(253/255,204/255,138/255),size=1.25)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size=12))+
  stat_summary(aes(group=1),fun=mean, geom="point", 
               color="black", size=3)+
  ylim(0,4)+geom_abline(intercept=summary(mod)$coefficients[1,1],slope=summary(mod)$coefficients[2,1],size=2)

# IplotCS_points<-ggplot(results_small, aes(x=comphet, y=S))+
#   geom_point(alpha=.25,color='red')+ylim(0,4)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#                                                      panel.background = element_blank(),axis.line = element_line(colour = "black"),axis.text=element_text(size=12))+
#   ylim(0,4)+geom_abline(intercept=summary(mod)$coefficients[1,1],slope=summary(mod)$coefficients[2,1],size=3)





#LMERs:

#now do stats on that shit and store the coefficients/pvals 
#to display as bargraphs or whatever:

#how do you want your data organized:
AS<-matrix(data=NA,ncol=6,nrow=length(models)*9,byrow=TRUE)
colnames(AS)<-c('slope','p','ChiSq','d','dd','model')
AC<-matrix(data=NA,ncol=6,nrow=length(models)*9,byrow=TRUE)
colnames(AC)<-c('slope','p','ChiSq','d','dd','model')
CS<-matrix(data=NA,ncol=6,nrow=length(models)*9,byrow=TRUE)
colnames(CS)<-c('slope','p','ChiSq','d','dd','model')
# Aminmax<-matrix(data=NA,ncol=2,nrow=length(models)*9,byrow=TRUE)
# colnames(Aminmax)<-c('min','max')
# Cminmax<-matrix(data=NA,ncol=2,nrow=length(models)*9,byrow=TRUE)
# colnames(Cminmax)<-c('min','max')
# Sminmax<-matrix(data=NA,ncol=2,nrow=length(models)*9,byrow=TRUE)
# colnames(Sminmax)<-c('min','max')
Cmean<-matrix(data=NA,ncol=1,nrow=length(models)*9,byrow=TRUE)
Smean<-matrix(data=NA,ncol=1,nrow=length(models)*9,byrow=TRUE)

count<-1
for(d in 1:length(dispersal)){
  for(s in 1:length(ddstrength)){
    for(m in 1:length(models)){
      #call and combine data:
      
      setwd(paste("~/Documents/submissions/NichePart2020/for Dryad/2020/code_editedsumm2019/neutral/",dispersal[d],'/',ddstrength[s],'/',models[m],sep=''))
      #setwd(paste('~/Documents/MATLAB/simulation2019Jan/neutral/',dispersal[d],'/',ddstrength[s],'/',models[m],sep=''))
      results_small<-read.csv(paste('results_small.csv',sep=''), header = TRUE, sep = ',')
      results_large<-read.csv(paste('results_large.csv',sep=''), header = TRUE, sep = ',')
      
      results_small<-cbind(smallpred,results_small)
      results_large<-cbind(largepred,results_large)
      
      #remove NaNs, round habitat
      results_small<-results_small[!is.nan(results_small$comphet),]
      results_small$meanC<-round(results_small$meanC,1)
      

      #Cminmax[count,1]<-min(results_small$comphet)
      #Cminmax[count,2]<-max(results_small$comphet)
      #Sminmax[count,1]<-min(results_small$S)
      #Sminmax[count,2]<-max(results_small$S)
      Cmean[count,1]<-median(results_small$comphet)
      Smean[count,1]<-median(results_small$S)
      
      #models:
      ASmod<-lmer(S~varC + (1|landscape:meanC),data=results_small,control=ctrl)
      ACmod<-lmer(comphet~varC + (1|landscape:meanC),data=results_small,control=ctrl)
      CSmod<-lmer(S~comphet + (1|landscape:meanC),data=results_small,control=ctrl)
      
      #compile:
      AS[count,1]<-summary(ASmod)$coefficients[2,1]
      AS[count,2]<-Anova(ASmod)[1,3]
      AS[count,3]<-Anova(ASmod)[1,1]
      AS[count,4]<-dispersal[d]
      AS[count,5]<-ddstrength[s]
      AS[count,6]<-models[m]
      
      AC[count,1]<-summary(ACmod)$coefficients[2,1]
      AC[count,2]<-Anova(ACmod)[1,3]
      AC[count,3]<-Anova(ACmod)[1,1]
      AC[count,4]<-dispersal[d]
      AC[count,5]<-ddstrength[s]
      AC[count,6]<-models[m]
      
      CS[count,1]<-summary(CSmod)$coefficients[2,1]
      CS[count,2]<-Anova(CSmod)[1,3]
      CS[count,3]<-Anova(CSmod)[1,1]
      CS[count,4]<-dispersal[d]
      CS[count,5]<-ddstrength[s]
      CS[count,6]<-models[m]
      count = count+1}
  }
}

#mins and maxes:
# min(Sminmax[,1])
# max(Sminmax[,2])
# min(Cminmax[,1])
# max(Cminmax[,2])
# min(Aminmax[,1])
# max(Aminmax[,2])

AS[,6][AS[,6]=='neutRdiffGneutC']<-'D'
AS[,6][AS[,6]=='neutRneutGdiffC']<-'I'
AC[,6][AC[,6]=='neutRdiffGneutC']<-'D'
AC[,6][AC[,6]=='neutRneutGdiffC']<-'I'
CS[,6][CS[,6]=='neutRdiffGneutC']<-'D'
CS[,6][CS[,6]=='neutRneutGdiffC']<-'I'

boxplot(abs(as.numeric(as.character(slope)))~model,data=AS,main="A/S")
boxplot(abs(as.numeric(as.character(slope)))~model,data=AC,main="A/C")
boxplot(abs(as.numeric(as.character(slope)))~model,data=CS,main="C/S")

allAS<-lm(I(abs(as.numeric(as.character(slope))))~model*d+model*dd,data=as.data.frame(AS))
round(Anova(allAS,type=3),2)
plot_model(allAS,type='emm',terms=c('model','d'),dot.size=6.5,
           col=c(rgb(197/255,27/255,138/255),rgb(250/255,159/255,181/255),rgb(253/255,224/255,221/255
           )),line.size=.75)+
  scale_x_discrete("type", expand=c(0,.5))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.key=element_blank(),
        axis.line = element_line(colour = "black",size=1.5),
        axis.text=element_text(size=16))

allAC<-lm(I(abs(as.numeric(as.character(slope))))~model*d+model*dd,data=as.data.frame(AC))
round(Anova(allAC,type=3),2)
plot_model(allAC,type='emm',terms=c('model','d'),dot.size=6.5,
           col=c(rgb(197/255,27/255,138/255),rgb(250/255,159/255,181/255),rgb(253/255,224/255,221/255
           )),line.size=.75)+
  scale_x_discrete("type", expand=c(0,.5))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.key=element_blank(),
        axis.line = element_line(colour = "black",size=1.5),
        axis.text=element_text(size=16))

allCS<-lm(I(abs(as.numeric(as.character(slope))))~model*d+model*dd,data=as.data.frame(CS))
round(Anova(allCS,type=3),2)
plot_model(allCS,type='emm',terms=c('model','d'),dot.size=6.5,
           col=c(rgb(197/255,27/255,138/255),rgb(250/255,159/255,181/255),rgb(253/255,224/255,221/255
           )),line.size=.75)+
  scale_x_discrete("type", expand=c(0,.5))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
       panel.background = element_blank(),
       legend.key=element_blank(),
       axis.line = element_line(colour = "black",size=1.5),
       axis.text=element_text(size=16))

#pvalue figs for supplement:
ASpvals<-lm(I(as.numeric(as.character(p)))~model*d+model*dd,data=as.data.frame(AS))
round(Anova(ASpvals,type=3),2)
plot_model(ASpvals,type='emm',terms=c('model','d'),dot.size=6.5,
           col=c(rgb(197,27,138),rgb(250,159,181),rgb(253,224,221
)),line.size=.75)+
  scale_x_discrete("type", expand=c(0,.5))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.key=element_blank(),
        axis.line = element_line(colour = "black",size=1.5),
        axis.text=element_text(size=16))+geom_hline(yintercept=0.05,size=1.5)

ACpvals<-lm(I(as.numeric(as.character(p)))~model*d+model*dd,data=as.data.frame(AC))
round(Anova(ACpvals,type=3),2)
plot_model(ACpvals,type='emm',terms=c('model','d'),dot.size=6.5,
           col=bpy.colors(n=3,alpha=1,cutoff.tails=.35),line.size=.75)+
  scale_x_discrete("type", expand=c(0,.5))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.key=element_blank(),
        axis.line = element_line(colour = "black",size=1.5),
        axis.text=element_text(size=16))+geom_hline(yintercept=0.05,size=1.5)

CSpvals<-lm(I(as.numeric(as.character(p)))~model*d+model*dd,data=as.data.frame(CS))
round(Anova(CSpvals,type=3),2)
plot_model(CSpvals,type='emm',terms=c('model','d'),dot.size=6.5,
           col=bpy.colors(n=3,alpha=1,cutoff.tails=.35),line.size=.75)+
  scale_x_discrete("type", expand=c(0,.5))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.key=element_blank(),
        axis.line = element_line(colour = "black",size=1.5),
        axis.text=element_text(size=16))+geom_hline(yintercept=.05,size=1.5)


#take outliers out:
keepI<-which(abs(as.numeric(as.character((CS[,1]))))<2.5)
Dres<-CS[CS[,6]=='D',]
Ires<-CS[keepI,]
outlierCS<-lm(I(abs(as.numeric(as.character(slope))))~model*d+model*dd+d*dd,data=as.data.frame(rbind(Dres,Ires)))
Anova(outlierCS)
plot_model(outlierCS,type='emm',terms=c('model','d'))





#here's some ways I did regression:
      
      #for plots:
      
      setwd(paste("~/Documents/submissions/NichePart2020/for Dryad/2020/code_editedsumm2019/neutral/",dispersal[d],'/',ddstrength[s],'/',models[m],sep=''))
      #setwd(paste('~/Documents/MATLAB/simulation2019Jan/neutral/',dispersal[d],'/',ddstrength[s],'/',models[m],sep=''))
      results_small<-read.csv(paste('results_small.csv',sep=''), header = TRUE, sep = ',')
      results_large<-read.csv(paste('results_large.csv',sep=''), header = TRUE, sep = ',')
      
      results_small<-cbind(smallpred,results_small)
      results_large<-cbind(largepred,results_large)
      
      results_small$meanC<-round(results_small$meanC,1)
      #linear mixed-effect models:
      noNans<-results_small[!is.nan(results_small$comphet),]
    
      theCS_small<-lmer(S~comphet +(1|landscape)+(1|meanC),data=noNans,na.action=na.omit)
      simulationOutput <- simulateResiduals(fittedModel = theAC_small, n = 750)
      plot(simulationOutput)
      Anova(theAC_small)
      
      #for violin plot:
      noNans$comphet<-as.factor(round(noNans$comphet,1))
      ggplot(noNans, aes(x=comphet, y=S)) +
        #geom_jitter(width=.4,alpha=0.25,color='blue')+
        geom_violin(alpha=0,color='red')+
        scale_x_discrete(limits=c(0,.1,.2,.3,.4),labels=c(0,.1,.2,.3,.4))

      
      noNans<-results_small[!is.nan(results_small$comphet),]
      theAC_small<-lmer(sqrt(comphet)~varC +(1|landscape)+(1|meanC),data=noNans,control=ctrl)
      
      
      
      linearplot<-plot_model(theAC_small,terms='varC',type='emm',show.data=TRUE)+aes(color='red')
      
      #since the assumptions aren't met for a gaussian distribution:
      results_small$S<-as.factor(round(results_small$S,1))
      noNans$comphet<-as.factor(round(noNans$comphet,1))
      theAS_small<-clmm(S~varC +(1|landscape)+(1|meanC),data=results_small,Hess=TRUE)
      ordinalplot<-plot_model(theAS_small,terms='varC',type='emm')
      summary(theAS_small)$coefficients
      
      theAC_small<-clmm(comphet~varC +(1|landscape)+(1|meanC),data=noNans,Hess=TRUE)
      ordinalplot<-plot_model(theAC_small,terms='varC',type='emm')
      summary(theAC_small)$coefficients
      
      
      
      
      

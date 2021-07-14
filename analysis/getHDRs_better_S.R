

#####call libraries and functions#####
########################################################################

library(nlme)
library(MuMIn)
library(ggplot2)
library(lsmeans)

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
ctrl<-lmeControl(opt='optim')#ctrl not implemented, just in case.
options(na.action=na.fail)

#####define habitat quality when species partition resources along a gradient#####
########################################################################

setwd("~/Documents/MATLAB/simulation/landscapes")
bp<-read.csv('breakpoints.csv',header=FALSE,sep=',') #where species curves 
            #intersect. i.e. defines the range of abiotic conditions for which 
            #species x performs better than any other species.
#np<-read.csv('np.csv',header=FALSE,sep=',') #niche optima for each species


#####get abiotic predictors#####
########################################################################

setwd("~/Documents/MATLAB/simulation_2018-2019")
smallpred<-read.csv('smallpred.csv',header = TRUE,sep=',')
largepred<-read.csv('largepred.csv',header = TRUE,sep = ',')
bysps<-c('annual-edge','perennial-int','annual-int','perennial-edge')
colors<-c('lightblue','black','blue','slategrey')
habitats = 0 #not working when set to 1....
plotit<-0

models<-c('neutRdiffGneutC','neutRneutGdiffC')#'neutRneutGneutC','neutRdiffGneutC','neutRneutGdiffC','neutRdiffGdiffC')
ddstrength<-c('_0_5','_1_0','_1_5')
dispersal<-c('adj','int','uni')#'adj','int','uni')
#dd<-c('dd=_1','dd=_9','dd=1_1','dd=1_9')#
#alpha<-c('alpha = _1','alpha = _9','alpha = 1_1','alpha = 1_9')#
lme_ab_coefs<-matrix(NA,ncol=8,nrow = length(models)*9,byrow=TRUE)
lme_comp_coefs<-matrix(NA,ncol=8,nrow = length(models)*9,byrow=TRUE)
cov_coefs<-matrix(NA,ncol=8,nrow = length(models)*9,byrow=TRUE)

#g_all_small<-list()
#g_all_large<-list()
#for(a in 1:length(alpha)){#length(alpha)){
if(plotit==1){
par(mfrow=c(2,2),
    mai=c(0.5,0.5,0.5,0.5))
}
for(d in 1:length(dispersal)){
for(s in 1:length(ddstrength)){
for(m in 1:length(models)){
  if(d==1&s==1&m==1){
    count=1
  }
setwd(paste('~/Documents/MATLAB/simulation2019Jan/neutral/',dispersal[d],'/',ddstrength[s],'/',models[m],sep=''))
results_small<-read.csv(paste('results_small.csv',sep=''), header = TRUE, sep = ',')
results_large<-read.csv(paste('results_large.csv',sep=''), header = TRUE, sep = ',')

results_small<-cbind(smallpred,results_small)
results_large<-cbind(largepred,results_large)


for(i in 1:20){
  results_small$hab[results_small$landscape==i]<-unlist(lapply(results_small$meanC[results_small$landscape==i],function(x) 
    if(x<=bp[i,1]){'annual-edge'}else if(x>bp[i,1]&x<=bp[i,2])
    {'perennial-int'}else if(x>bp[i,2]&x<=bp[i,3]){'annual-int'}
    else{'perennial-edge'}))}


#results_small<-results_small[results_small$landscape != 7 & results_small$landscape != 14 & results_small$landscape != 16,]
#results_large<-results_large[results_large$landscape != 7 & results_large$landscape != 14 & results_large$landscape != 16,]
if(habitats == 1){
for(sps in 1:length(bysps)){
  results<-results_small[results_small$hab == bysps[sps],]
  
  lme_ab<-lme(fit~varC,~varC|landscape,data=results,control=ctrl,na.action=na.omit)
  #diagnostics(lme_ab)
  if(plotit==1){
  plot(fit~varC,data=results,cex=.7,main=paste(models[m],',',ddstrength[s],',',dispersal[d]),
       xlab='abiotic heterogeneity',ylab='species richness (# sps/subplot)',ylim=c(0,4))
  for(i in 1:20){abline(coef(lme_ab)[i,1],coef(lme_ab)[i,2],col=colors[sps])}
  abline(summary(lme_ab)$coefficients$fixed[1],
         summary(lme_ab)$coefficients$fixed[2],col='red',lwd=3)
  }
  
  lme_comp<-lme(fit~comphet,~comphet|landscape,data=results,control=ctrl,na.action=na.omit)
  #diagnostics(lme_comp)
  if(plotit==1){
  plot(fit~comphet,data=results,cex=.7,main=paste(models[m],',',ddstrength[s],',',dispersal[d]),
       xlab='competitive heterogeneity',ylab='species richness (# sps/subplot)',ylim=c(0,4))
  for(i in 1:20){abline(coef(lme_comp)[i,1],coef(lme_comp)[i,2],col=colors[sps])}
  abline(summary(lme_comp)$coefficients$fixed[1],
         summary(lme_comp)$coefficients$fixed[2],col='red',lwd=3)
  }
  cov<-lme(sqrt(comphet)~varC,~varC|landscape,data=results,control=ctrl,na.action=na.omit)
  #diagnostics(cov)
  if(plotit==1){
  plot(sqrt(comphet)~varC,data=results,cex=.7,main=paste(models[m],',',ddstrength[s],',',dispersal[d]),
       xlab='abiotic heterogeneity',ylab='competitive heterogeneity',axes=FALSE)
  axis(1,at=c(0,.5,1,1.5,2,2.5,3),labels=c(0,.5,1,1.5,2,2.5,3),cex.axis=1)
  axis(2,at=c(0,.1,.2,.3,.4,.5,.6,.7),labels=round((c(0,.1,.2,.3,.4,.5,.6,.7))^2,2),cex.axis=1)
  box()
  for(i in 1:20){abline(coef(cov)[i,1],coef(cov)[i,2],col=colors[sps])}
  abline(summary(cov)$coefficients$fixed[1],
         summary(cov)$coefficients$fixed[2],col='red',lwd=3)}}
  
  
  
  }else{
    
lme_ab<-lme(S~varC,~varC|landscape/hab,data=results_small,control=ctrl,na.action=na.omit)
#diagnostics(lme_ab)
#pull coefficients:
ae<-grep('/annual-edge',rownames(coef(lme_ab)),fixed=TRUE)
pi<-grep('/perennial-int',rownames(coef(lme_ab)),fixed=TRUE)
ai<-grep('/annual-int',rownames(coef(lme_ab)),fixed=TRUE)
pe<-grep('/perennial-edge',rownames(coef(lme_ab)),fixed=TRUE)

ab_fixed<-round(lme_ab$coefficients$fixed[2],2)
ab_ae<-round(mean(coef(lme_ab)[c(ae),2]),4)
ab_ai<-round(mean(coef(lme_ab)[c(ai),2]),4)
ab_pe<-round(mean(coef(lme_ab)[c(pe),2]),4)
ab_pi<-round(mean(coef(lme_ab)[c(pi),2]),4)
lme_ab_coefs[count,]<-c(ab_fixed,ab_ae,ab_pi,ab_ai,ab_pe,paste(dispersal[d]),paste(ddstrength[s]),paste(models[m]))

if(plotit==1){
plot(S~varC,data=results_small,cex=.7,main=paste(models[m],',',ddstrength[s],',',dispersal[d]),
     xlab='abiotic heterogeneity',ylab='species richness (# sps/subplot)',ylim=c(0,4))
for(i in c(ae)){abline(coef(lme_ab)[i,1],coef(lme_ab)[i,2],col='lightblue')}
for(i in c(ai)){abline(coef(lme_ab)[i,1],coef(lme_ab)[i,2],col='blue')}
for(i in c(pe)){abline(coef(lme_ab)[i,1],coef(lme_ab)[i,2],col='slategrey')}
for(i in c(pi)){abline(coef(lme_ab)[i,1],coef(lme_ab)[i,2],col='black')}
abline(summary(lme_ab)$coefficients$fixed[1],
       summary(lme_ab)$coefficients$fixed[2],col='red',lwd=3)
}

lme_comp<-lme(S~comphet,~comphet|landscape/hab,data=results_small,control=ctrl,na.action=na.omit)
#diagnostics(lme_comp)
#pull coefficients
ae<-grep('/annual-edge',rownames(coef(lme_comp)),fixed=TRUE)
pi<-grep('/perennial-int',rownames(coef(lme_comp)),fixed=TRUE)
ai<-grep('/annual-int',rownames(coef(lme_comp)),fixed=TRUE)
pe<-grep('/perennial-edge',rownames(coef(lme_comp)),fixed=TRUE)

comp_fixed<-round(lme_comp$coefficients$fixed[2],2)
comp_ae<-round(mean(coef(lme_comp)[c(ae),2]),4)
comp_ai<-round(mean(coef(lme_comp)[c(ai),2]),4)
comp_pe<-round(mean(coef(lme_comp)[c(pe),2]),4)
comp_pi<-round(mean(coef(lme_comp)[c(pi),2]),4)
lme_comp_coefs[count,]<-c(comp_fixed,comp_ae,comp_pi,comp_ai,comp_pe,paste(dispersal[d]),paste(ddstrength[s]),paste(models[m]))

if(plotit==1){
plot(S~comphet,data=results_small,cex=.7,main=paste(models[m],',',ddstrength[s],',',dispersal[d]),
     xlab='competitive heterogeneity',ylab='species richness (# sps/subplot)',ylim=c(0,4))
for(i in c(ae)){abline(coef(lme_comp)[i,1],coef(lme_comp)[i,2],col='lightblue')}
for(i in c(ai)){abline(coef(lme_comp)[i,1],coef(lme_comp)[i,2],col='blue')}
for(i in c(pe)){abline(coef(lme_comp)[i,1],coef(lme_comp)[i,2],col='slategrey')}
for(i in c(pi)){abline(coef(lme_comp)[i,1],coef(lme_comp)[i,2],col='black')}
abline(summary(lme_comp)$coefficients$fixed[1],
       summary(lme_comp)$coefficients$fixed[2],col='red',lwd=3)
}

cov<-lme(sqrt(comphet)~varC,~varC|landscape/hab,data=results_small,control=ctrl,na.action=na.omit)
#diagnostics(cov)
#pull coefficients
ae<-grep('/annual-edge',rownames(coef(cov)),fixed=TRUE)
pi<-grep('/perennial-int',rownames(coef(cov)),fixed=TRUE)
ai<-grep('/annual-int',rownames(coef(cov)),fixed=TRUE)
pe<-grep('/perennial-edge',rownames(coef(cov)),fixed=TRUE)

cov_fixed<-round(cov$coefficients$fixed[2],2)
cov_ae<-round(mean(coef(cov)[c(ae),2]),4)
cov_ai<-round(mean(coef(cov)[c(ai),2]),4)
cov_pe<-round(mean(coef(cov)[c(pe),2]),4)
cov_pi<-round(mean(coef(cov)[c(pi),2]),4)
cov_coefs[count,]<-c(cov_fixed,cov_ae,cov_pi,cov_ai,cov_pe,paste(dispersal[d]),paste(ddstrength[s]),paste(models[m]))

if(plotit==1){
plot(sqrt(comphet)~varC,data=results_small,cex=.7,main=paste(models[m],',',ddstrength[s],',',dispersal[d]),
     xlab='abiotic heterogeneity',ylab='competitive heterogeneity',axes=FALSE)
axis(1,at=c(0,.5,1,1.5,2,2.5,3),labels=c(0,.5,1,1.5,2,2.5,3),cex.axis=1)
axis(2,at=c(sqrt(min(!is.na(results_small$comphet))), sqrt(max(!is.na(results_small$comphet)))),labels=round(c(min(results_small$comphet), max(results_small$comphet))^2,2),cex.axis=1)
box()
for(i in c(ae)){abline(coef(cov)[i,1],coef(cov)[i,2],col='lightblue')}
for(i in c(ai)){abline(coef(cov)[i,1],coef(cov)[i,2],col='blue')}
for(i in c(pe)){abline(coef(cov)[i,1],coef(cov)[i,2],col='slategrey')}
for(i in c(pi)){abline(coef(cov)[i,1],coef(cov)[i,2],col='black')}
abline(summary(cov)$coefficients$fixed[1],
       summary(cov)$coefficients$fixed[2],col='red',lwd=3)
}
  }
count = count+1
}
}
}
colnames(lme_ab_coefs)<-c('fixed','ae','pi','ai','pe','dispersal','ddstrength','model')
colnames(lme_comp_coefs)<-c('fixed','ae','pi','ai','pe','dispersal','ddstrength','model')
colnames(cov_coefs)<-c('fixed','ae','pi','ai','pe','dispersal','ddstrength','model')
abhet<-as.data.frame(lme_ab_coefs)
comphet<-as.data.frame(lme_comp_coefs)
abcomp_cov<-as.data.frame(cov_coefs)

#COMPARE JUST EF AND AMC MODELS:
#data<-as.numeric(rownames(abhet[as.character(abhet$model)=='neutRneutGdiffC' | as.character(abhet$model)=='neutRdiffGneutC',]))

#building models with other categorical predictors
#abhetvecF<-abs(as.numeric(as.character(abhet$fixed)))[data]
#comphetvecF<-abs(as.numeric(as.character(comphet$fixed)))[data]
#abcompvecF<-abs(as.numeric(as.character(abcomp_cov$fixed)))[data]
#dispersalvec<-as.factor(as.character(abcomp_cov$dispersal))[data]
#ddstrengthvec<-as.factor(as.character(abcomp_cov$ddstrength))[data]
#modelvec<-as.factor(as.character(abcomp_cov$model))[data]

#fix your output:
hdrdata<-matrix(abs(as.numeric(matrix(unlist(abhet[,1:5]),ncol=5,nrow=length(models)*9,byrow=FALSE))),ncol=5,nrow=length(models)*9,byrow=FALSE)
disp<-as.factor(as.character(abcomp_cov$dispersal))
dd<-as.factor(as.character(abcomp_cov$ddstrength))
mods<-as.factor(as.character(abcomp_cov$model))
scompdata<-matrix(abs(as.numeric(matrix(unlist(comphet[,1:5]),ncol=5,nrow=length(models)*9,byrow=FALSE))),ncol=5,nrow=length(models)*9,byrow=FALSE)
covdata<-matrix(abs(as.numeric(matrix(unlist(abcomp_cov[,1:5]),ncol=5,nrow=length(models)*9,byrow=FALSE))),ncol=5,nrow=length(models)*9,byrow=FALSE)
alldat<-data.frame(hdrdata,scompdata,covdata,disp,dd,mods)
colnames(alldat)<-c('fixedHDR','aeHDR','piHDR','aiHDR','peHDR','fixedCDR','aeCDR','piCDR','aiCDR','peCDR','fixedCOV','aeCOV','piCOV','aiCOV','peCOV','disp','dd','mods')

#plot normalized HDRs together for comparison.
plot(scale(fixedHDR)~jitter(fixedCOV),data=alldat,col='black',main='abiotic heterogeneity',
     xlab='strength of covariance between abiotic and 
     competitive heterogeneity (abs(slope))',
     ylab='strength of HDR (abs(slope))',xlim=c(0,.125),ylim=c(-2,5))
fixedMOD<-lm(scale(fixedHDR)~fixedCOV,data=alldat[alldat$mods=='neutRdiffGneutC',])
abline(coef(fixedMOD)[1],coef(fixedMOD)[2],col='black',lwd=3,lty=2)

points(jitter(alldat$fixedCOV),scale(alldat$fixedCDR),
       pch=16,col='black',main='competitive heterogeneity')
fixedMOD<-lm(scale(fixedCDR)~fixedCOV,data=alldat)
abline(coef(fixedMOD)[1],coef(fixedMOD)[2],col='black',lwd=3,lty=1)
legend('topright',legend = c('abiotic HDRs','competitive HDRs'),col=c('black','black'),pch=c(1,16),lty=c(2,1),lwd=c(2,2))

plot(fixedHDR~jitter(fixedCOV),data=alldat,col='springgreen',
     xlab='strength of relationship between abiotic and 
     competitive heterogeneity',
     ylab='strength of the relationship between 
     abiotic heterogeneity and species richness',pch=16)
fixedMOD<-lm(fixedHDR~fixedCOV,data=alldat)
abline(coef(fixedMOD)[1],coef(fixedMOD)[2],col='springgreen',lwd=3,lty=1)

plot(fixedCDR~jitter(fixedCOV),data=alldat,col='lavenderblush4',
     xlab='strength of relationship between abiotic and 
     competitive heterogeneity',
     ylab='strength of the relationship between 
     competitive heterogeneity and species richness',pch=16)
fixedMOD<-lm(fixedCDR~fixedCOV,data=alldat)
abline(coef(fixedMOD)[1],coef(fixedMOD)[2],col='lavenderblush4',lwd=3,lty=1)
#legend('topright',legend = c('abiotic HDRs','competitive HDRs'),col=c('black','black'),pch=c(1,16),lty=c(2,1),lwd=c(2,2))

colors<-palette(c("blue","red"))

plot(fixedCDR~jitter(fixedCOV),data=alldat[alldat$mods=='neutRdiffGneutC',],col='blue',main='comp heterogeneity',
     xlab='strength of covariance between abiotic and 
     competitive heterogeneity (abs(slope))',
     ylab='strength of HDR (abs(slope))',xlim=c(0,.125),ylim=c(0,12),pch=as.numeric(alldat$disp[alldat$mods=='neutRdiffGneutC']))
points(jitter(alldat$fixedCOV[alldat$mods=='neutRneutGdiffC']),alldat$fixedCDR[alldat$mods=='neutRneutGdiffC'],
       pch=as.numeric(alldat$disp[alldat$mods=='neutRneutGdiffC']),col='red')
legend('bottomright',legend=c('EF','AMC'),pch=c(15,15),pt.cex=2,col=c('blue','red'),bty='n')
fixedMOD<-lm(fixedCDR~fixedCOV,data=alldat[alldat$disp=='adj',])
abline(coef(fixedMOD)[1],coef(fixedMOD)[2],col='black',lwd=3,lty=1)
fixedMOD<-lm(fixedCDR~fixedCOV,data=alldat[alldat$disp=='int',])
abline(coef(fixedMOD)[1],coef(fixedMOD)[2],col='black',lwd=3,lty=2)
fixedMOD<-lm(fixedCDR~fixedCOV,data=alldat[alldat$disp=='uni',])
abline(coef(fixedMOD)[1],coef(fixedMOD)[2],col='black',lwd=3,lty=3)
legend('topright',legend=c('adjacent','intermediate','universal'),lty=c(1,2,3),lwd=c(2,2,2),bty='n',pch=c(1,2,3))

plot(fixedHDR~jitter(fixedCOV),data=alldat,col=mods,main='abiotic heterogeneity',
     xlab='strength of covariance between abiotic and 
     competitive heterogeneity (abs(slope))',
     ylab='strength of HDR (abs(slope))',xlim=c(0,.125),pch=16)
#points(jitter(alldat$fixedCOV[alldat$mods=='neutRneutGdiffC']),alldat$fixedHDR[alldat$mods=='neutRneutGdiffC'],
#       pch=16,col='red')
legend('bottomright',legend=c('EF','AMC'),pch=c(15,15),pt.cex=2,col=c('blue','red'),bty='n')
fixedMOD<-lm(fixedHDR~fixedCOV,data=alldat)
abline(coef(fixedMOD)[1],coef(fixedMOD)[2],col='black',lwd=3,lty=1)

fixedMOD<-lm(fixedHDR~fixedCOV,data=alldat[alldat$disp=='adj',])
abline(coef(fixedMOD)[1],coef(fixedMOD)[2],col='black',lwd=3,lty=1)
fixedMOD<-lm(fixedHDR~fixedCOV,data=alldat[alldat$disp=='int',])
abline(coef(fixedMOD)[1],coef(fixedMOD)[2],col='black',lwd=3,lty=2)
fixedMOD<-lm(fixedHDR~fixedCOV,data=alldat[alldat$disp=='uni',])
abline(coef(fixedMOD)[1],coef(fixedMOD)[2],col='black',lwd=3,lty=3)
legend('bottomright',legend=c('adjacent','intermediate','universal'),lty=c(1,2,3),lwd=c(3,3,3),bty='n')


plot(scale(fixedCDR)~fixedCOV,data=alldat[alldat$disp=='adj',],col=mods,main='abiotic heterogeneity',
     xlab='strength of covariance between abiotic and 
     competitive heterogeneity (abs(slope))',
     ylab='strength of CDR (abs(slope))',xlim=c(0,.125),ylim=c(-2,5),pch=16)
fixedMOD<-lm(scale(fixedCDR)~fixedCOV,data=alldat[alldat$disp=='adj',])
abline(coef(fixedMOD)[1],coef(fixedMOD)[2],col='black',lwd=3,lty=1)

points(jitter(alldat$fixedCOV[alldat$disp=='int']),scale(alldat$fixedCDR[alldat$disp=='int']),
       pch=16,col=mods)
fixedMOD<-lm(scale(fixedCDR)~fixedCOV,data=alldat[alldat$disp=='int',])
abline(coef(fixedMOD)[1],coef(fixedMOD)[2],col='black',lwd=3,lty=2)

points(jitter(alldat$fixedCOV[alldat$disp=='uni']),scale(alldat$fixedCDR[alldat$disp=='uni']),
       pch=16,col=mods)
fixedMOD<-lm(scale(fixedCDR)~fixedCOV,data=alldat[alldat$disp=='uni',])
abline(coef(fixedMOD)[1],coef(fixedMOD)[2],col='black',lwd=3,lty=3)


plot(scale(fixedCDR)~fixedCOV,data=alldat[alldat$disp=='int',],col=mods,main='abiotic heterogeneity',
     xlab='strength of covariance between abiotic and 
     competitive heterogeneity (abs(slope))',
     ylab='strength of CDR (abs(slope))',xlim=c(0,.125),ylim=c(-2,5),pch=1)

plot(scale(fixedCDR)~fixedCOV,data=alldat[alldat$disp=='uni',],col=mods,main='abiotic heterogeneity',
     xlab='strength of covariance between abiotic and 
     competitive heterogeneity (abs(slope))',
     ylab='strength of CDR (abs(slope))',xlim=c(0,.125),ylim=c(-2,5))


#HDRS
#fixed
plot(fixedHDR~fixedCOV,data=alldat,
     pch=16,col='black',main='abiotic heterogeneity',
     xlab='strength of covariance between abiotic and 
     competitive heterogeneity (abs(slope))',
     ylab='strength of HDR (abs(slope))',ylim=c(0,1),xlim=c(0,.125))
fixedMOD<-lm(fixedHDR~fixedCOV,data=alldat)
abline(coef(fixedMOD)[1],coef(fixedMOD)[2],col='red',lwd=3)


#ae
points(alldat$aeCOV,alldat$aeHDR,col='lightblue',pch=16)
aeMOD<-lm(aeHDR~aeCOV,data=alldat)
abline(coef(aeMOD)[1],coef(aeMOD)[2],col='lightblue')

#pi
points(alldat$piCOV,alldat$piHDR,col='black',pch=16)
piMOD<-lm(piHDR~piCOV,data=alldat)
abline(coef(piMOD)[1],coef(piMOD)[2],col='black')

#ai
points(alldat$aiCOV,alldat$aiHDR,col='blue',pch=16)
aiMOD<-lm(aiHDR~aiCOV,data=alldat)
abline(coef(aiMOD)[1],coef(aiMOD)[2],col='blue')

#pe
points(alldat$peCOV,alldat$peHDR,col='slategrey',pch=16)
peMOD<-lm(peHDR~peCOV,data=alldat)
abline(coef(peMOD)[1],coef(peMOD)[2],col='slategrey')
#ALL OF THESE ARE SIGNIFICANT


#CDRS
#fixed
plot(fixedCDR~fixedCOV,data=alldat,
     pch=16,col='black',main='competitive heterogeneity',
     xlab='strength of covariance between abiotic and 
     competitive heterogeneity (abs(slope))',
     ylab='strength of CDR (abs(slope))',ylim=c(0,15),xlim=c(0,.125))
fixedMOD<-lm(fixedCDR~fixedCOV,data=alldat)
abline(coef(fixedMOD)[1],coef(fixedMOD)[2],col='red',lwd=3)


#ae
points(alldat$aeCOV,alldat$aeCDR,col='lightblue',pch=16)
aeMOD<-lm(aeCDR~aeCOV,data=alldat)
abline(coef(aeMOD)[1],coef(aeMOD)[2],col='lightblue')

#pi
points(alldat$piCOV,alldat$piCDR,col='black',pch=16)
piMOD<-lm(piCDR~piCOV,data=alldat)
abline(coef(piMOD)[1],coef(piMOD)[2],col='black')

#ai
points(alldat$aiCOV,alldat$aiCDR,col='blue',pch=16)
aiMOD<-lm(aiCDR~aiCOV,data=alldat)
abline(coef(aiMOD)[1],coef(aiMOD)[2],col='blue')

#pe
points(alldat$peCOV,alldat$peCDR,col='slategrey',pch=16)
peMOD<-lm(peCDR~peCOV,data=alldat)
abline(coef(peMOD)[1],coef(peMOD)[2],col='slategrey')
#ANNUALS HAVE A STRONGER NEG SLOPE HERE... which means as the covariance between 
#abhet and comphet increases, the strength of the S~comphet relationship decreases, but
#only in annual species.
#BUT THEYRE NOT SIGNIFICANT SO....

###THIS IS ALL YOU NEED TO LOOK AT!

#fixed
plot(scale(fixedCOV)~disp,data=alldat,col='gold')#,ylim=c(0,.08))
sig<-aov(fixedCDR~disp,data=alldat)
pairs(lsmeans(sig,~disp))
legend('bottomright',legend=c('n.s.'),bty='n')
#NS
plot(fixedCOV~dd,data=alldat,col='lightpink')#,ylim=c(0,.08))
sig<-aov(fixedCOV~dd,data=alldat)
pairs(lsmeans(sig,~dd))
legend('topright',legend=c('n.s.'),bty='n')
#NS
plot(scale(fixedCOV)~mods,data=alldat,border=c('blue','red'),pch=16,ylim=c(-1,3))
sig<-aov(fixedCOV~mods,data=alldat)
pairs(lsmeans(sig,~mods))
#SIGNIFICANT
legend('topright',legend=c('p < 0.0001'),bty='n')

#Main fixed effects
plot(scale(fixedHDR)~disp,data=alldat,col='gold')
sig<-aov(fixedHDR~disp,data=alldat)
pairs(lsmeans(sig,~disp))
#NS
plot(scale(fixedHDR)~dd,data=alldat,col='lightpink')
sig<-aov(fixedHDR~dd,data=alldat)
pairs(lsmeans(sig,~dd))
#NS
plot(scale(fixedHDR)~mods,data=alldat,col=c('blue','red'),pch=16,add=TRUE)
sig<-aov(fixedHDR~mods,data=alldat)
pairs(lsmeans(sig,~mods))
legend('topright',legend=c('p < 0.0001'),bty='n')
#SIGNIFICANT

plot(scale(fixedCDR)~disp,data=alldat,col='gold')
sig<-aov(fixedCDR~disp,data=alldat)
pairs(lsmeans(sig,~disp))
#NS
plot(fixedCDR~dd,data=alldat,col='lightpink')
sig<-aov(fixedCDR~dd,data=alldat)
pairs(lsmeans(sig,~dd))
#NS
plot(fixedCDR~mods,data=alldat,col=c('blue','red'),pch=16)
sig<-aov(fixedCDR~mods,data=alldat)
pairs(lsmeans(sig,~mods))
#NS
legend('topleft',legend=c('p = 0.49'),bty='n')


#fixed interaction effects:

boxplot(fixedCOV~disp,data=alldat[alldat$mods=='neutRdiffGneutC',],ylim=c(0,.1), col='blue')
boxplot(fixedCOV~disp,data=alldat[alldat$mods=='neutRneutGdiffC',],add=TRUE,col='red')

boxplot(fixedCOV~dd,data=alldat[alldat$mods=='neutRdiffGneutC',],ylim=c(0,.1),col='blue')
boxplot(fixedCOV~dd,data=alldat[alldat$mods=='neutRneutGdiffC',],add=TRUE,col='red')

boxplot(fixedHDR~disp,data=alldat[alldat$mods=='neutRdiffGneutC',], col='blue')
boxplot(fixedHDR~disp,data=alldat[alldat$mods=='neutRneutGdiffC',],add=TRUE,col='red')
intmod<-aov(scale(fixedHDR)~as.factor(disp)*as.factor(mods),data=alldat)
summary(intmod)
boxplot(fixedHDR~dd,data=alldat[alldat$mods=='neutRdiffGneutC',],ylim=c(0,.5),col='blue')
boxplot(fixedHDR~dd,data=alldat[alldat$mods=='neutRneutGdiffC',],add=TRUE,col='red')
intmod<-aov(fixedHDR~as.factor(dd)*as.factor(mods),data=alldat)
summary(intmod)


boxplot(fixedCDR~disp,data=alldat[alldat$mods=='neutRdiffGneutC',],ylim=c(0,12), col='blue',main='dispersal x assembly mechanism')
legend('bottomright',legend=c('EF','AMC'),pch=c(15,15),pt.cex=2,col=c('blue','red'),bty='n')
boxplot(fixedCDR~disp,data=alldat[alldat$mods=='neutRneutGdiffC',],add=TRUE,col='red')
legend('topright',legend=c('p=0.0015'),bty='n')
intmod<-aov(fixedCDR~as.factor(disp)*as.factor(mods),data=alldat)
summary(intmod)

boxplot(scale(fixedCDR)~dd,data=alldat[alldat$mods=='neutRdiffGneutC',],ylim=c(-2,2),col='blue')
boxplot(scale(fixedCDR)~dd,data=alldat[alldat$mods=='neutRneutGdiffC',],add=TRUE,col='red')

#ae
plot(aeCOV~disp,data=alldat,col='gold')
sig<-aov(aeCOV~disp,data=alldat)
pairs(lsmeans(sig,~disp))
plot(aeCOV~dd,data=alldat,col='lightpink')
sig<-aov(aeCOV~dd,data=alldat)
pairs(lsmeans(sig,~dd))
plot(aeCOV~mods,data=alldat,col=c('blue','red'))
sig<-aov(aeCOV~mods,data=alldat)
pairs(lsmeans(sig,~mods))

#pi
plot(piCOV~disp,data=alldat)
sig<-aov(piCOV~disp,data=alldat)
pairs(lsmeans(sig,~disp))
plot(piCOV~dd,data=alldat)
sig<-aov(piCOV~dd,data=alldat)
pairs(lsmeans(sig,~dd))
plot(piCOV~mods,data=alldat,col=c('blue','red'))
sig<-aov(piCOV~mods,data=alldat)
pairs(lsmeans(sig,~mods))

#ai
plot(aiCOV~disp,data=alldat)
sig<-aov(aiCOV~disp,data=alldat)
pairs(lsmeans(sig,~disp))
plot(aiCOV~dd,data=alldat)
sig<-aov(aiCOV~dd,data=alldat)
pairs(lsmeans(sig,~dd))
plot(aiCOV~mods,data=alldat)
sig<-aov(aiCOV~mods,data=alldat)
pairs(lsmeans(sig,~mods))

#pe
plot(peCOV~disp,data=alldat)
sig<-aov(aiCOV~disp,data=alldat)
pairs(lsmeans(sig,~disp))
plot(peCOV~dd,data=alldat)
sig<-aov(peCOV~dd,data=alldat)
pairs(lsmeans(sig,~dd))
plot(peCOV~mods,data=alldat)
sig<-aov(peCOV~mods,data=alldat)
pairs(lsmeans(sig,~mods))
#ONLY MOD IS SIGNIFICANT

#interaction effects:
#FIXED
summary(lm(fixedCOV~mods*disp,data=alldat)) #YES INTERACTION EFFECT IN UNI DISP:
sig<-aov(fixedCOV~mods*disp,data=alldat)
summary(sig)
#DISP AND MODS:DISP ARE SIG
plot(fixedCOV~disp,data=alldat)
plot(fixedCOV~disp,data=alldat[alldat$mods=='neutRneutGdiffC',])
plot(fixedCOV~disp,data=alldat[alldat$mods=='neutRdiffGneutC',])

summary(lm(fixedCOV~mods*dd,data=alldat)) #no interaction effects
sig<-aov(fixedCOV~mods*dd,data=alldat)
summary(sig)
#NS

#ae
summary(lm(aeCOV~mods*disp,data=alldat))
sig<-aov(aeCOV~mods*disp,data=alldat)
summary(sig)
#NS
summary(lm(aeCOV~mods*dd,data=alldat))
sig<-aov(aeCOV~mods*dd,data=alldat)
summary(sig)
#NS

#pi
summary(lm(piCOV~mods*disp,data=alldat))
#sig here
sig<-aov(piCOV~mods*disp,data=alldat)
summary(sig)
#but not here.
plot(piCOV~disp,data=alldat[alldat$mods=='neutRneutGdiffC',])
plot(piCOV~disp,data=alldat[alldat$mods=='neutRdiffGneutC',])

summary(lm(piCOV~mods*dd,data=alldat))
sig<-aov(piCOV~mods*dd,data=alldat)
summary(sig)
#NS

#ai
summary(lm(aiCOV~mods*disp,data=alldat))#EVERYTHING IS SIGNIFICANT
sig<-aov(aiCOV~mods*disp,data=alldat)
summary(sig)
plot(aiCOV~disp,data=alldat)
plot(aiCOV~disp,data=alldat[alldat$mods=='neutRneutGdiffC',])
plot(aiCOV~disp,data=alldat[alldat$mods=='neutRdiffGneutC',])

summary(lm(aiCOV~mods*dd,data=alldat))
sig<-aov(aiCOV~mods*dd,data=alldat)
summary(sig)
#NS

#pe
summary(lm(peCOV~mods*disp,data=alldat))
sig<-aov(peCOV~mods*disp,data=alldat)#unid is sigdiff pos
summary(sig)
plot(peCOV~disp,data=alldat)

summary(lm(peCOV~mods*dd,data=alldat))
sig<-aov(peCOV~mods*dd,data=alldat)#unid is sigdiff pos
summary(sig)

########################################################################################

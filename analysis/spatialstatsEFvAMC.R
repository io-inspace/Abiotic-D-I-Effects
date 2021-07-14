library(scales)
library(ggplot2)

load("/Users/s_catella/Documents/submissions/NichePart2020/for Dryad/2020/code_editedsumm2019/neutral/results/ISAR_outputmatrix.RData")
load("/Users/s_catella/Documents/submissions/NichePart2020/for Dryad/2020/code_editedsumm2019/neutral/results/Pc_outputmatrix.RData")

ISARdat<-as.data.frame(ISAR_outputmatrix)
Pdat<-as.data.frame(Pc_outputmatrix)
unique(ISARdat$model)
unique(Pdat$model)

#OK THE MAIN FIGURE IS GOING TO BE STATS FOR THE EXAMPLE POINT PATTERNS SHOWN IN THE METHODS FIG.

EF<-ISARdat[ISARdat$model == 'neutRdiffGneutC' & ISARdat$dispersal == 'int' & ISARdat$ddstrength == '_1_0',]
EFnums<-as.vector(rep(0,720))
for(i in 1:20){
  EFnums<-cbind(EFnums,as.numeric(as.character(EF[,i+1])))
}

AMC<-ISARdat[ISARdat$model == 'neutRneutGdiffC'& ISARdat$dispersal == 'int' & ISARdat$ddstrength == '_1_0',]
AMCnums<-as.vector(rep(0,720))
for(i in 1:20){
  AMCnums<-cbind(AMCnums,as.numeric(as.character(AMC[,i+1])))
}

N<-ISARdat[ISARdat$model == 'neutRneutGneutC'& ISARdat$dispersal == 'int' & ISARdat$ddstrength == '_1_0',]
Nnums<-as.vector(rep(0,720))
for(i in 1:20){
  Nnums<-cbind(Nnums,as.numeric(as.character(N[,i+1])))
}

#you want to normalize your radius to between 0 and 1:
rad<-seq(0,10,by=.5) #round(sapply(X=seq(0,10,.5),function(x) (x-0)/(10-0)),3)
par(mfrow=c(1,2),mar=c(2,2.5,2,2))
b1<-boxplot(EFnums,col=alpha(rgb(215/255,48/255,31/255),.5),cex.axis=1.25,cex=.5,whisklty=1,names=as.factor(rad),ylim=c(0,1),ylab=NA,xlab=NA,staplecol='red',medcol='red',outpch=16,outcol='red',whiskcol='red')
b2<-boxplot(Nnums,col=alpha('black',.5),whisklty=1,cex=.5,names=as.factor(rad),axes=FALSE,staplecol='black',medcol='black',outpch=16,outcol='black',whiskcol='black',add=TRUE)
b3<-boxplot(AMCnums,col=alpha(rgb(253/255,204/255,138/255),.5),whisklty=1,cex=.5,names=as.factor(rad),axes=FALSE,staplecol='orange',medcol='orange',outpch=16,outcol='orange',whiskcol='orange',add=TRUE)

par(mfrow=c(1,1))
datapoints<-data.frame(c(b1$stats[3,],b2$stats[3,],b3$stats[3,]),
                           c(b1$stats[1,],b2$stats[1,],b3$stats[1,]),
                           c(b1$stats[5,],b2$stats[5,],b3$stats[5,]),
                           rep(rad,3),c(rep('Direct',21),rep('Neutral',21),
                                        rep('Indirect',21)))
colnames(datapoints)<-c('median','lower','upper','radius','model')
ggplot(data=datapoints, aes(x=radius, y=median, color=model)) + 
  geom_line(size=1) + geom_ribbon(aes(ymin=lower, ymax=upper), 
                                           linetype=2, alpha=0.1)+ylim(0,1)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
         panel.background = element_blank(),
         axis.line = element_line(colour = "black"),
         axis.text=element_text(size=12))+scale_color_manual(values = c(rgb(215/255,48/255,31/255),rgb(253/255,204/255,138/255),'black'))


EF<-Pdat[Pdat$model == 'neutRdiffGneutC' & Pdat$dispersal == 'int' & Pdat$ddstrength == '_1_0',]
EFnums<-as.vector(rep(0,720))
for(i in 1:20){
  EFnums<-cbind(EFnums,as.numeric(as.character(EF[,i+1])))
}

AMC<-Pdat[Pdat$model == 'neutRneutGdiffC'& Pdat$dispersal == 'int' & Pdat$ddstrength == '_1_0',]
AMCnums<-as.vector(rep(0,720))
for(i in 1:20){
  AMCnums<-cbind(AMCnums,as.numeric(as.character(AMC[,i+1])))
}

N<-Pdat[Pdat$model == 'neutRneutGneutC'& Pdat$dispersal == 'int' & Pdat$ddstrength == '_1_0',]
Nnums<-as.vector(rep(0,720))
for(i in 1:20){
  Nnums<-cbind(Nnums,as.numeric(as.character(N[,i+1])))
}

b1<-boxplot(EFnums,col=alpha(rgb(215/255,48/255,31/255),.5),cex.axis=1.25,cex=.5,whisklty=1,names=as.factor(rad),ylim=c(0,1),ylab=NA,xlab=NA,staplecol='red',medcol='red',outpch=16,outcol='red',whiskcol='red',print=FALSE)
b2<-boxplot(Nnums,col=alpha('black',.5),whisklty=1,cex=.5,names=as.factor(rad),axes=FALSE,staplecol='black',medcol='black',outpch=16,outcol='black',whiskcol='black',add=TRUE,print=FALSE)
b3<-boxplot(AMCnums,col=alpha(rgb(253/255,204/255,138/255),.5),whisklty=1,cex=.5,names=as.factor(rad),axes=FALSE,staplecol='orange',medcol='orange',outpch=16,outcol='orange',whiskcol='orange',add=TRUE,print=FALSE)

datapoints<-data.frame(c(b1$stats[3,],b2$stats[3,],b3$stats[3,]),
                       c(b1$stats[1,],b2$stats[1,],b3$stats[1,]),
                       c(b1$stats[5,],b2$stats[5,],b3$stats[5,]),
                       rep(rad,3),c(rep('Direct',21),rep('Neutral',21),
                                    rep('Indirect',21)))
colnames(datapoints)<-c('median','lower','upper','radius','model')
par(mfrow=c(1,1))
ggplot(data=datapoints, aes(x=radius, y=median, color=model)) + 
  geom_line(size=1) + geom_ribbon(aes(ymin=lower, ymax=upper), 
                                  linetype=2, alpha=0.1)+ylim(0,1)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size=12))+scale_color_manual(values = c(rgb(215/255,48/255,31/255),rgb(253/255,204/255,138/255),'black'))


+
  legend('topright',legend=c('Direct','Indirect','Neutral'),col=c('red','orange','black'),pch=15,bty='n',title="Abiotic effect", cex=1.25)


#This is code to make pretty supplemental figures across parameter space:

#pos1

EF<-ISARdat[ISARdat$model == 'neutRdiffGneutC' & ISARdat$dispersal == 'adj' & ISARdat$ddstrength == '_0_5',]
EFnums<-as.vector(rep(0,720))
for(i in 1:20){
  EFnums<-cbind(EFnums,as.numeric(as.character(EF[,i+1])))
}

AMC<-ISARdat[ISARdat$model == 'neutRneutGdiffC'& ISARdat$dispersal == 'adj' & ISARdat$ddstrength == '_0_5',]
AMCnums<-as.vector(rep(0,720))
for(i in 1:20){
 AMCnums<-cbind(AMCnums,as.numeric(as.character(AMC[,i+1])))
}

N<-ISARdat[ISARdat$model == 'neutRneutGneutC'& ISARdat$dispersal == 'adj' & ISARdat$ddstrength == '_0_5',]
Nnums<-as.vector(rep(0,720))
for(i in 1:20){
  Nnums<-cbind(Nnums,as.numeric(as.character(N[,i+1])))
}

#you want to normalize your radius to between 0 and 1:
rad<-seq(0,10,by=.5) #round(sapply(X=seq(0,10,.5),function(x) (x-0)/(10-0)),3)
par(mfrow=c(3,3),mar=c(2,2,1,1))
boxplot(EFnums,col=alpha(rgb(215/255,48/255,31/255),.5),cex.axis=1.1,cex=.5,whisklty=1,names=as.factor(rad),ylim=c(0,1),ylab=NA,xlab=NA,staplecol='red',medcol='red',outpch=16,outcol='red',whiskcol='red',xaxt='n')
boxplot(Nnums,col=alpha('black',.5),whisklty=1,cex=.5,names=as.factor(rad),axes=FALSE,staplecol='black',medcol='black',outpch=16,outcol='black',whiskcol='black',add=TRUE)
boxplot(AMCnums,col=alpha(rgb(253/255,204/255,138/255),.5),whisklty=1,cex=.5,names=as.factor(rad),axes=FALSE,staplecol='orange',medcol='orange',outpch=16,outcol='orange',whiskcol='orange',add=TRUE)

#pos2

EF<-ISARdat[ISARdat$model == 'neutRdiffGneutC' & ISARdat$dispersal == 'adj' & ISARdat$ddstrength == '_1_0',]
EFnums<-as.vector(rep(0,720))
for(i in 1:20){
  EFnums<-cbind(EFnums,as.numeric(as.character(EF[,i+1])))
}

AMC<-ISARdat[ISARdat$model == 'neutRneutGdiffC'& ISARdat$dispersal == 'adj' & ISARdat$ddstrength == '_1_0',]
AMCnums<-as.vector(rep(0,720))
for(i in 1:20){
  AMCnums<-cbind(AMCnums,as.numeric(as.character(AMC[,i+1])))
}

N<-ISARdat[ISARdat$model == 'neutRneutGneutC'& ISARdat$dispersal == 'adj' & ISARdat$ddstrength == '_1_0',]
Nnums<-as.vector(rep(0,720))
for(i in 1:20){
  Nnums<-cbind(Nnums,as.numeric(as.character(N[,i+1])))
}

boxplot(EFnums,col=alpha(rgb(215/255,48/255,31/255),.5),cex.axis=1.1,cex=.5,
        whisklty=1,names=as.factor(rad),ylim=c(0,1),ylab=NA,xlab=NA,
        staplecol='red',medcol='red',outpch=16,outcol='red',
        whiskcol='red',xaxt='n',yaxt='n')
boxplot(Nnums,col=alpha('black',.5),whisklty=1,cex=.5,
        names=as.factor(rad),axes=FALSE,staplecol='black',
        medcol='black',outpch=16,outcol='black',whiskcol='black',
        add=TRUE)
boxplot(AMCnums,col=alpha(rgb(253/255,204/255,138/255),.5),whisklty=1,cex=.5,
        names=as.factor(rad),axes=FALSE,staplecol='orange',
        medcol='orange',outpch=16,outcol='orange',whiskcol='orange',add=TRUE)


#pos3

EF<-ISARdat[ISARdat$model == 'neutRdiffGneutC' & ISARdat$dispersal == 'adj' & ISARdat$ddstrength == '_1_5',]
EFnums<-as.vector(rep(0,720))
for(i in 1:20){
  EFnums<-cbind(EFnums,as.numeric(as.character(EF[,i+1])))
}

AMC<-ISARdat[ISARdat$model == 'neutRneutGdiffC'& ISARdat$dispersal == 'adj' & ISARdat$ddstrength == '_1_5',]
AMCnums<-as.vector(rep(0,720))
for(i in 1:20){
  AMCnums<-cbind(AMCnums,as.numeric(as.character(AMC[,i+1])))
}

N<-ISARdat[ISARdat$model == 'neutRneutGneutC'& ISARdat$dispersal == 'adj' & ISARdat$ddstrength == '_1_5',]
Nnums<-as.vector(rep(0,720))
for(i in 1:20){
  Nnums<-cbind(Nnums,as.numeric(as.character(N[,i+1])))
}

boxplot(EFnums,col=alpha(rgb(215/255,48/255,31/255),.5),cex.axis=1.1,cex=.5,whisklty=1,
        names=as.factor(rad),ylim=c(0,1),ylab=NA,xlab=NA,
        staplecol='red',medcol='red',outpch=16,outcol='red',
        whiskcol='red',xaxt='n',yaxt='n')
boxplot(Nnums,col=alpha('black',.5),whisklty=1,cex=.5,
                                                 names=as.factor(rad),axes=FALSE,staplecol='black',
                                                 medcol='black',outpch=16,outcol='black',whiskcol='black',add=TRUE)

boxplot(AMCnums,col=alpha(rgb(253/255,204/255,138/255),.5),whisklty=1,cex=.5,
        names=as.factor(rad),axes=FALSE,staplecol='orange',medcol='orange',
        outpch=16,outcol='orange',whiskcol='orange',add=TRUE)

#pos4

EF<-ISARdat[ISARdat$model == 'neutRdiffGneutC' & ISARdat$dispersal == 'int' & ISARdat$ddstrength == '_0_5',]
EFnums<-as.vector(rep(0,720))
for(i in 1:20){
  EFnums<-cbind(EFnums,as.numeric(as.character(EF[,i+1])))
}

AMC<-ISARdat[ISARdat$model == 'neutRneutGdiffC'& ISARdat$dispersal == 'int' & ISARdat$ddstrength == '_0_5',]
AMCnums<-as.vector(rep(0,720))
for(i in 1:20){
  AMCnums<-cbind(AMCnums,as.numeric(as.character(AMC[,i+1])))
}

N<-ISARdat[ISARdat$model == 'neutRneutGneutC'& ISARdat$dispersal == 'int' & ISARdat$ddstrength == '_0_5',]
Nnums<-as.vector(rep(0,720))
for(i in 1:20){
  Nnums<-cbind(Nnums,as.numeric(as.character(N[,i+1])))
}

boxplot(EFnums,col=alpha(rgb(215/255,48/255,31/255),.5),cex.axis=1.1,cex=.5,
        whisklty=1,names=as.factor(rad),ylim=c(0,1),ylab=NA,xlab=NA,
        staplecol='red',medcol='red',outpch=16,outcol='red',
        whiskcol='red',xaxt='n')
boxplot(Nnums,col=alpha('black',.5),whisklty=1,cex=.5,
        names=as.factor(rad),axes=FALSE,staplecol='black',
        medcol='black',outpch=16,outcol='black',whiskcol='black',
        add=TRUE)
boxplot(AMCnums,col=alpha(rgb(253/255,204/255,138/255),.5),whisklty=1,cex=.5,
        names=as.factor(rad),axes=FALSE,staplecol='orange',
        medcol='orange',outpch=16,outcol='orange',whiskcol='orange',add=TRUE)

#pos5

EF<-ISARdat[ISARdat$model == 'neutRdiffGneutC' & ISARdat$dispersal == 'int' & ISARdat$ddstrength == '_1_0',]
EFnums<-as.vector(rep(0,720))
for(i in 1:20){
EFnums<-cbind(EFnums,as.numeric(as.character(EF[,i+1])))
}

AMC<-ISARdat[ISARdat$model == 'neutRneutGdiffC'& ISARdat$dispersal == 'int' & ISARdat$ddstrength == '_1_0',]
AMCnums<-as.vector(rep(0,720))
for(i in 1:20){
  AMCnums<-cbind(AMCnums,as.numeric(as.character(AMC[,i+1])))
}

N<-ISARdat[ISARdat$model == 'neutRneutGneutC'& ISARdat$dispersal == 'int' & ISARdat$ddstrength == '_1_0',]
Nnums<-as.vector(rep(0,720))
for(i in 1:20){
  Nnums<-cbind(Nnums,as.numeric(as.character(N[,i+1])))
}


# EFAMC<-ISARdat[ISARdat$model == 'neutRdiffGdiffC'& ISARdat$dispersal == 'int' & ISARdat$ddstrength == '_1_0',]
# EFAMCnums<-as.vector(rep(0,720))
# for(i in 1:20){
#   EFAMCnums<-cbind(EFAMCnums,as.numeric(as.character(EFAMC[,i+1])))
# }
# 
# par(cex.axis=1.25,mar=c(2,2.5,1,1),cex=.75)
# #boxplot(Nnums,col='grey',names=as.factor(rad),staplecol='black',medcol='black',outpch=16,outcol='black',whiskcol='black',add=TRUE)
boxplot(EFnums,col=alpha(rgb(215/255,48/255,31/255),.5),cex.axis=1.1,cex=.5,
        whisklty=1,names=as.factor(rad),ylim=c(0,1),ylab=NA,xlab=NA,
        staplecol='red',medcol='red',outpch=16,outcol='red',
        whiskcol='red',xaxt='n',yaxt='n')
boxplot(Nnums,col=alpha('black',.5),whisklty=1,cex=.5,
        names=as.factor(rad),axes=FALSE,staplecol='black',
        medcol='black',outpch=16,outcol='black',whiskcol='black',
        add=TRUE)
boxplot(AMCnums,col=alpha(rgb(253/255,204/255,138/255),.5),whisklty=1,cex=.5,
        names=as.factor(rad),axes=FALSE,staplecol='orange',
        medcol='orange',outpch=16,outcol='orange',whiskcol='orange',add=TRUE)

legend('bottomright',legend=c('D','I','N'),col=c('red','orange','black'),
       pch=15,bty='n',title="Abiotic Effect")

#boxplot(EFAMCnums,col=alpha('purple',.2),names=as.factor(rad),ylim=c(0,1),ylab=NA,xlab=NA,staplecol='purple',medcol='purple',outpch=16,outcol='purple',whiskcol='purple',add=TRUE)

#pos6

EF<-ISARdat[ISARdat$model == 'neutRdiffGneutC' & ISARdat$dispersal == 'int' & ISARdat$ddstrength == '_1_5',]
EFnums<-as.vector(rep(0,720))
for(i in 1:20){
  EFnums<-cbind(EFnums,as.numeric(as.character(EF[,i+1])))
}

AMC<-ISARdat[ISARdat$model == 'neutRneutGdiffC'& ISARdat$dispersal == 'int' & ISARdat$ddstrength == '_1_5',]
AMCnums<-as.vector(rep(0,720))
for(i in 1:20){
  AMCnums<-cbind(AMCnums,as.numeric(as.character(AMC[,i+1])))
}

N<-ISARdat[ISARdat$model == 'neutRneutGneutC'& ISARdat$dispersal == 'int' & ISARdat$ddstrength == '_1_5',]
Nnums<-as.vector(rep(0,720))
for(i in 1:20){
  Nnums<-cbind(Nnums,as.numeric(as.character(N[,i+1])))
}

boxplot(EFnums,col=alpha(rgb(215/255,48/255,31/255),.5),cex.axis=1.1,cex=.5,
        whisklty=1,names=as.factor(rad),ylim=c(0,1),ylab=NA,xlab=NA,
        staplecol='red',medcol='red',outpch=16,outcol='red',
        whiskcol='red',xaxt='n',yaxt='n')
boxplot(Nnums,col=alpha('black',.5),whisklty=1,cex=.5,
        names=as.factor(rad),axes=FALSE,staplecol='black',
        medcol='black',outpch=16,outcol='black',whiskcol='black',
        add=TRUE)

boxplot(AMCnums,col=alpha(rgb(253/255,204/255,138/255),.5),whisklty=1,cex=.5,
        names=as.factor(rad),axes=FALSE,staplecol='orange',
        medcol='orange',outpch=16,outcol='orange',whiskcol='orange',add=TRUE)

#pos7

EF<-ISARdat[ISARdat$model == 'neutRdiffGneutC' & ISARdat$dispersal == 'uni' & ISARdat$ddstrength == '_0_5',]
EFnums<-as.vector(rep(0,720))
for(i in 1:20){
  EFnums<-cbind(EFnums,as.numeric(as.character(EF[,i+1])))
}

AMC<-ISARdat[ISARdat$model == 'neutRneutGdiffC'& ISARdat$dispersal == 'uni' & ISARdat$ddstrength == '_0_5',]
AMCnums<-as.vector(rep(0,720))
for(i in 1:20){
  AMCnums<-cbind(AMCnums,as.numeric(as.character(AMC[,i+1])))
}

N<-ISARdat[ISARdat$model == 'neutRneutGneutC'& ISARdat$dispersal == 'uni' & ISARdat$ddstrength == '_0_5',]
Nnums<-as.vector(rep(0,720))
for(i in 1:20){
  Nnums<-cbind(Nnums,as.numeric(as.character(N[,i+1])))
}

boxplot(EFnums,col=alpha(rgb(215/255,48/255,31/255),.5),cex.axis=1.1,cex=.5,
        whisklty=1,names=as.factor(rad),ylim=c(0,1),ylab=NA,xlab=NA,
        staplecol='red',medcol='red',outpch=16,outcol='red',
        whiskcol='red')
boxplot(Nnums,col=alpha('black',.5),whisklty=1,cex=.5,
        names=as.factor(rad),axes=FALSE,staplecol='black',
        medcol='black',outpch=16,outcol='black',whiskcol='black',
        add=TRUE)
boxplot(AMCnums,col=alpha(rgb(253/255,204/255,138/255),.5),whisklty=1,cex=.5,
        names=as.factor(rad),axes=FALSE,staplecol='orange',
        medcol='orange',outpch=16,outcol='orange',whiskcol='orange',add=TRUE)
#pos8

EF<-ISARdat[ISARdat$model == 'neutRdiffGneutC' & ISARdat$dispersal == 'uni' & ISARdat$ddstrength == '_1_0',]
EFnums<-as.vector(rep(0,720))
for(i in 1:20){
  EFnums<-cbind(EFnums,as.numeric(as.character(EF[,i+1])))
}

AMC<-ISARdat[ISARdat$model == 'neutRneutGdiffC'& ISARdat$dispersal == 'uni' & ISARdat$ddstrength == '_1_0',]
AMCnums<-as.vector(rep(0,720))
for(i in 1:20){
  AMCnums<-cbind(AMCnums,as.numeric(as.character(AMC[,i+1])))
}

N<-ISARdat[ISARdat$model == 'neutRneutGneutC'& ISARdat$dispersal == 'uni' & ISARdat$ddstrength == '_1_0',]
Nnums<-as.vector(rep(0,720))
for(i in 1:20){
  Nnums<-cbind(Nnums,as.numeric(as.character(N[,i+1])))
}

boxplot(EFnums,col=alpha(rgb(215/255,48/255,31/255),.5),cex.axis=1.1,cex=.5,
        whisklty=1,names=as.factor(rad),ylim=c(0,1),ylab=NA,xlab=NA,
        staplecol='red',medcol='red',outpch=16,outcol='red',
        whiskcol='red',yaxt='n')
boxplot(Nnums,col=alpha('black',.5),whisklty=1,cex=.5,
        names=as.factor(rad),axes=FALSE,staplecol='black',
        medcol='black',outpch=16,outcol='black',whiskcol='black',
        add=TRUE)
boxplot(AMCnums,col=alpha(rgb(253/255,204/255,138/255),.5),whisklty=1,cex=.5,
        names=as.factor(rad),axes=FALSE,staplecol='orange',
        medcol='orange',outpch=16,outcol='orange',whiskcol='orange',add=TRUE)


#pos9

EF<-ISARdat[ISARdat$model == 'neutRdiffGneutC' & ISARdat$dispersal == 'uni' & ISARdat$ddstrength == '_1_5',]
EFnums<-as.vector(rep(0,720))
for(i in 1:20){
  EFnums<-cbind(EFnums,as.numeric(as.character(EF[,i+1])))
}

AMC<-ISARdat[ISARdat$model == 'neutRneutGdiffC'& ISARdat$dispersal == 'uni' & ISARdat$ddstrength == '_1_5',]
AMCnums<-as.vector(rep(0,720))
for(i in 1:20){
  AMCnums<-cbind(AMCnums,as.numeric(as.character(AMC[,i+1])))
}

N<-ISARdat[ISARdat$model == 'neutRneutGneutC'& ISARdat$dispersal == 'uni' & ISARdat$ddstrength == '_1_5',]
Nnums<-as.vector(rep(0,720))
for(i in 1:20){
  Nnums<-cbind(Nnums,as.numeric(as.character(N[,i+1])))
}
boxplot(EFnums,col=alpha(rgb(215/255,48/255,31/255),.5),cex.axis=1.1,cex=.5,
        whisklty=1,names=as.factor(rad),ylim=c(0,1),ylab=NA,xlab=NA,
        staplecol='red',medcol='red',outpch=16,outcol='red',
        whiskcol='red',yaxt='n')
boxplot(Nnums,col=alpha('black',.5),whisklty=1,cex=.5,
                                        names=as.factor(rad),axes=FALSE,staplecol='black',
                                        medcol='black',outpch=16,outcol='black',whiskcol='black',
                                        add=TRUE)

boxplot(AMCnums,col=alpha(rgb(253/255,204/255,138/255),.5),whisklty=1,cex=.5,
        names=as.factor(rad),axes=FALSE,staplecol='orange',
        medcol='orange',outpch=16,outcol='orange',whiskcol='orange',add=TRUE)

#Pdat:
#This is code to make pretty supplemental figures across parameter space:

#pos1

EF<-Pdat[Pdat$model == 'neutRdiffGneutC' & Pdat$dispersal == 'adj' & Pdat$ddstrength == '_0_5',]
EFnums<-as.vector(rep(0,720))
for(i in 1:20){
  EFnums<-cbind(EFnums,as.numeric(as.character(EF[,i+1])))
}

AMC<-Pdat[Pdat$model == 'neutRneutGdiffC'& Pdat$dispersal == 'adj' & Pdat$ddstrength == '_0_5',]
AMCnums<-as.vector(rep(0,720))
for(i in 1:20){
  AMCnums<-cbind(AMCnums,as.numeric(as.character(AMC[,i+1])))
}

N<-Pdat[Pdat$model == 'neutRneutGneutC'& Pdat$dispersal == 'adj' & Pdat$ddstrength == '_0_5',]
Nnums<-as.vector(rep(0,720))
for(i in 1:20){
  Nnums<-cbind(Nnums,as.numeric(as.character(N[,i+1])))
}

#you want to normalize your radius to between 0 and 1:
rad<-seq(0,10,by=.5) #round(sapply(X=seq(0,10,.5),function(x) (x-0)/(10-0)),3)
par(mfrow=c(3,3),mar=c(2,2,1,1))
boxplot(EFnums,col=alpha(rgb(215/255,48/255,31/255),.5),cex.axis=1.1,cex=.5,whisklty=1,names=as.factor(rad),ylim=c(0,1),ylab=NA,xlab=NA,staplecol='red',medcol='red',outpch=16,outcol='red',whiskcol='red',xaxt='n')
boxplot(Nnums,col=alpha('black',.5),whisklty=1,cex=.5,names=as.factor(rad),axes=FALSE,staplecol='black',medcol='black',outpch=16,outcol='black',whiskcol='black',add=TRUE)
boxplot(AMCnums,col=alpha(rgb(253/255,204/255,138/255),.5),whisklty=1,cex=.5,names=as.factor(rad),axes=FALSE,staplecol='orange',medcol='orange',outpch=16,outcol='orange',whiskcol='orange',add=TRUE)

#pos2

EF<-Pdat[Pdat$model == 'neutRdiffGneutC' & Pdat$dispersal == 'adj' & Pdat$ddstrength == '_1_0',]
EFnums<-as.vector(rep(0,720))
for(i in 1:20){
  EFnums<-cbind(EFnums,as.numeric(as.character(EF[,i+1])))
}

AMC<-Pdat[Pdat$model == 'neutRneutGdiffC'& Pdat$dispersal == 'adj' & Pdat$ddstrength == '_1_0',]
AMCnums<-as.vector(rep(0,720))
for(i in 1:20){
  AMCnums<-cbind(AMCnums,as.numeric(as.character(AMC[,i+1])))
}

N<-Pdat[Pdat$model == 'neutRneutGneutC'& Pdat$dispersal == 'adj' & Pdat$ddstrength == '_1_0',]
Nnums<-as.vector(rep(0,720))
for(i in 1:20){
  Nnums<-cbind(Nnums,as.numeric(as.character(N[,i+1])))
}

boxplot(EFnums,col=alpha(rgb(215/255,48/255,31/255),.5),cex.axis=1.1,cex=.5,
        whisklty=1,names=as.factor(rad),ylim=c(0,1),ylab=NA,xlab=NA,
        staplecol='red',medcol='red',outpch=16,outcol='red',
        whiskcol='red',xaxt='n',yaxt='n')
boxplot(Nnums,col=alpha('black',.5),whisklty=1,cex=.5,
        names=as.factor(rad),axes=FALSE,staplecol='black',
        medcol='black',outpch=16,outcol='black',whiskcol='black',
        add=TRUE)
boxplot(AMCnums,col=alpha(rgb(253/255,204/255,138/255),.5),whisklty=1,cex=.5,
        names=as.factor(rad),axes=FALSE,staplecol='orange',
        medcol='orange',outpch=16,outcol='orange',whiskcol='orange',add=TRUE)


#pos3

EF<-Pdat[Pdat$model == 'neutRdiffGneutC' & Pdat$dispersal == 'adj' & Pdat$ddstrength == '_1_5',]
EFnums<-as.vector(rep(0,720))
for(i in 1:20){
  EFnums<-cbind(EFnums,as.numeric(as.character(EF[,i+1])))
}

AMC<-Pdat[Pdat$model == 'neutRneutGdiffC'& Pdat$dispersal == 'adj' & Pdat$ddstrength == '_1_5',]
AMCnums<-as.vector(rep(0,720))
for(i in 1:20){
  AMCnums<-cbind(AMCnums,as.numeric(as.character(AMC[,i+1])))
}

N<-Pdat[Pdat$model == 'neutRneutGneutC'& Pdat$dispersal == 'adj' & Pdat$ddstrength == '_1_5',]
Nnums<-as.vector(rep(0,720))
for(i in 1:20){
  Nnums<-cbind(Nnums,as.numeric(as.character(N[,i+1])))
}

boxplot(EFnums,col=alpha(rgb(215/255,48/255,31/255),.5),cex.axis=1.1,cex=.5,whisklty=1,
        names=as.factor(rad),ylim=c(0,1),ylab=NA,xlab=NA,
        staplecol='red',medcol='red',outpch=16,outcol='red',
        whiskcol='red',xaxt='n',yaxt='n')
boxplot(Nnums,col=alpha('black',.5),whisklty=1,cex=.5,
        names=as.factor(rad),axes=FALSE,staplecol='black',
        medcol='black',outpch=16,outcol='black',whiskcol='black',add=TRUE)
boxplot(AMCnums,col=alpha(rgb(253/255,204/255,138/255),.5),whisklty=1,cex=.5,
        names=as.factor(rad),axes=FALSE,staplecol='orange',medcol='orange',
        outpch=16,outcol='orange',whiskcol='orange',add=TRUE)

#pos4

EF<-Pdat[Pdat$model == 'neutRdiffGneutC' & Pdat$dispersal == 'int' & Pdat$ddstrength == '_0_5',]
EFnums<-as.vector(rep(0,720))
for(i in 1:20){
  EFnums<-cbind(EFnums,as.numeric(as.character(EF[,i+1])))
}

AMC<-Pdat[Pdat$model == 'neutRneutGdiffC'& Pdat$dispersal == 'int' & Pdat$ddstrength == '_0_5',]
AMCnums<-as.vector(rep(0,720))
for(i in 1:20){
  AMCnums<-cbind(AMCnums,as.numeric(as.character(AMC[,i+1])))
}

N<-Pdat[Pdat$model == 'neutRneutGneutC'& Pdat$dispersal == 'int' & Pdat$ddstrength == '_0_5',]
Nnums<-as.vector(rep(0,720))
for(i in 1:20){
  Nnums<-cbind(Nnums,as.numeric(as.character(N[,i+1])))
}

boxplot(EFnums,col=alpha(rgb(215/255,48/255,31/255),.5),cex.axis=1.1,cex=.5,
        whisklty=1,names=as.factor(rad),ylim=c(0,1),ylab=NA,xlab=NA,
        staplecol='red',medcol='red',outpch=16,outcol='red',
        whiskcol='red',xaxt='n')
boxplot(Nnums,col=alpha('black',.5),whisklty=1,cex=.5,
        names=as.factor(rad),axes=FALSE,staplecol='black',
        medcol='black',outpch=16,outcol='black',whiskcol='black',
        add=TRUE)
boxplot(AMCnums,col=alpha(rgb(253/255,204/255,138/255),.5),whisklty=1,cex=.5,
        names=as.factor(rad),axes=FALSE,staplecol='orange',
        medcol='orange',outpch=16,outcol='orange',whiskcol='orange',add=TRUE)

#pos5

EF<-Pdat[Pdat$model == 'neutRdiffGneutC' & Pdat$dispersal == 'int' & Pdat$ddstrength == '_1_0',]
EFnums<-as.vector(rep(0,720))
for(i in 1:20){
  EFnums<-cbind(EFnums,as.numeric(as.character(EF[,i+1])))
}

AMC<-Pdat[Pdat$model == 'neutRneutGdiffC'& Pdat$dispersal == 'int' & Pdat$ddstrength == '_1_0',]
AMCnums<-as.vector(rep(0,720))
for(i in 1:20){
  AMCnums<-cbind(AMCnums,as.numeric(as.character(AMC[,i+1])))
}

N<-Pdat[Pdat$model == 'neutRneutGneutC'& Pdat$dispersal == 'int' & Pdat$ddstrength == '_1_0',]
Nnums<-as.vector(rep(0,720))
for(i in 1:20){
  Nnums<-cbind(Nnums,as.numeric(as.character(N[,i+1])))
}


# EFAMC<-Pdat[Pdat$model == 'neutRdiffGdiffC'& Pdat$dispersal == 'int' & Pdat$ddstrength == '_1_0',]
# EFAMCnums<-as.vector(rep(0,720))
# for(i in 1:20){
#   EFAMCnums<-cbind(EFAMCnums,as.numeric(as.character(EFAMC[,i+1])))
# }
# 
# par(cex.axis=1.25,mar=c(2,2.5,1,1),cex=.75)
# #boxplot(Nnums,col='grey',names=as.factor(rad),staplecol='black',medcol='black',outpch=16,outcol='black',whiskcol='black',add=TRUE)
boxplot(EFnums,col=alpha(rgb(215/255,48/255,31/255),.5),cex.axis=1.1,cex=.5,
        whisklty=1,names=as.factor(rad),ylim=c(0,1),ylab=NA,xlab=NA,
        staplecol='red',medcol='red',outpch=16,outcol='red',
        whiskcol='red',xaxt='n',yaxt='n')
boxplot(Nnums,col=alpha('black',.5),whisklty=1,cex=.5,
        names=as.factor(rad),axes=FALSE,staplecol='black',
        medcol='black',outpch=16,outcol='black',whiskcol='black',
        add=TRUE)
boxplot(AMCnums,col=alpha(rgb(253/255,204/255,138/255),.5),whisklty=1,cex=.5,
        names=as.factor(rad),axes=FALSE,staplecol='orange',
        medcol='orange',outpch=16,outcol='orange',whiskcol='orange',add=TRUE)

legend('topright',legend=c('D','I','N'),col=c('red','orange','black'),
       pch=15,bty='n',title="Abiotic Effect")

#boxplot(EFAMCnums,col=alpha('purple',.2),names=as.factor(rad),ylim=c(0,1),ylab=NA,xlab=NA,staplecol='purple',medcol='purple',outpch=16,outcol='purple',whiskcol='purple',add=TRUE)

#pos6

EF<-Pdat[Pdat$model == 'neutRdiffGneutC' & Pdat$dispersal == 'int' & Pdat$ddstrength == '_1_5',]
EFnums<-as.vector(rep(0,720))
for(i in 1:20){
  EFnums<-cbind(EFnums,as.numeric(as.character(EF[,i+1])))
}

AMC<-Pdat[Pdat$model == 'neutRneutGdiffC'& Pdat$dispersal == 'int' & Pdat$ddstrength == '_1_5',]
AMCnums<-as.vector(rep(0,720))
for(i in 1:20){
  AMCnums<-cbind(AMCnums,as.numeric(as.character(AMC[,i+1])))
}

N<-Pdat[Pdat$model == 'neutRneutGneutC'& Pdat$dispersal == 'int' & Pdat$ddstrength == '_1_5',]
Nnums<-as.vector(rep(0,720))
for(i in 1:20){
  Nnums<-cbind(Nnums,as.numeric(as.character(N[,i+1])))
}

boxplot(EFnums,col=alpha(rgb(215/255,48/255,31/255),.5),cex.axis=1.1,cex=.5,
        whisklty=1,names=as.factor(rad),ylim=c(0,1),ylab=NA,xlab=NA,
        staplecol='red',medcol='red',outpch=16,outcol='red',
        whiskcol='red',xaxt='n',yaxt='n')
boxplot(Nnums,col=alpha('black',.5),whisklty=1,cex=.5,
        names=as.factor(rad),axes=FALSE,staplecol='black',
        medcol='black',outpch=16,outcol='black',whiskcol='black',
        add=TRUE)
boxplot(AMCnums,col=alpha(rgb(253/255,204/255,138/255),.5),whisklty=1,cex=.5,
        names=as.factor(rad),axes=FALSE,staplecol='orange',
        medcol='orange',outpch=16,outcol='orange',whiskcol='orange',add=TRUE)

#pos7

EF<-Pdat[Pdat$model == 'neutRdiffGneutC' & Pdat$dispersal == 'uni' & Pdat$ddstrength == '_0_5',]
EFnums<-as.vector(rep(0,720))
for(i in 1:20){
  EFnums<-cbind(EFnums,as.numeric(as.character(EF[,i+1])))
}

AMC<-Pdat[Pdat$model == 'neutRneutGdiffC'& Pdat$dispersal == 'uni' & Pdat$ddstrength == '_0_5',]
AMCnums<-as.vector(rep(0,720))
for(i in 1:20){
  AMCnums<-cbind(AMCnums,as.numeric(as.character(AMC[,i+1])))
}

N<-Pdat[Pdat$model == 'neutRneutGneutC'& Pdat$dispersal == 'uni' & Pdat$ddstrength == '_0_5',]
Nnums<-as.vector(rep(0,720))
for(i in 1:20){
  Nnums<-cbind(Nnums,as.numeric(as.character(N[,i+1])))
}

boxplot(EFnums,col=alpha(rgb(215/255,48/255,31/255),.5),cex.axis=1.1,cex=.5,
        whisklty=1,names=as.factor(rad),ylim=c(0,1),ylab=NA,xlab=NA,
        staplecol='red',medcol='red',outpch=16,outcol='red',
        whiskcol='red')
boxplot(Nnums,col=alpha('black',.5),whisklty=1,cex=.5,
        names=as.factor(rad),axes=FALSE,staplecol='black',
        medcol='black',outpch=16,outcol='black',whiskcol='black',
        add=TRUE)
boxplot(AMCnums,col=alpha(rgb(253/255,204/255,138/255),.5),whisklty=1,cex=.5,
        names=as.factor(rad),axes=FALSE,staplecol='orange',
        medcol='orange',outpch=16,outcol='orange',whiskcol='orange',add=TRUE)

#pos8

EF<-Pdat[Pdat$model == 'neutRdiffGneutC' & Pdat$dispersal == 'uni' & Pdat$ddstrength == '_1_0',]
EFnums<-as.vector(rep(0,720))
for(i in 1:20){
  EFnums<-cbind(EFnums,as.numeric(as.character(EF[,i+1])))
}

AMC<-Pdat[Pdat$model == 'neutRneutGdiffC'& Pdat$dispersal == 'uni' & Pdat$ddstrength == '_1_0',]
AMCnums<-as.vector(rep(0,720))
for(i in 1:20){
  AMCnums<-cbind(AMCnums,as.numeric(as.character(AMC[,i+1])))
}

N<-Pdat[Pdat$model == 'neutRneutGneutC'& Pdat$dispersal == 'uni' & Pdat$ddstrength == '_1_0',]
Nnums<-as.vector(rep(0,720))
for(i in 1:20){
  Nnums<-cbind(Nnums,as.numeric(as.character(N[,i+1])))
}

boxplot(EFnums,col=alpha(rgb(215/255,48/255,31/255),.5),cex.axis=1.1,cex=.5,
        whisklty=1,names=as.factor(rad),ylim=c(0,1),ylab=NA,xlab=NA,
        staplecol='red',medcol='red',outpch=16,outcol='red',
        whiskcol='red',yaxt='n')
boxplot(Nnums,col=alpha('black',.5),whisklty=1,cex=.5,
        names=as.factor(rad),axes=FALSE,staplecol='black',
        medcol='black',outpch=16,outcol='black',whiskcol='black',
        add=TRUE)
boxplot(AMCnums,col=alpha(rgb(253/255,204/255,138/255),.5),whisklty=1,cex=.5,
        names=as.factor(rad),axes=FALSE,staplecol='orange',
        medcol='orange',outpch=16,outcol='orange',whiskcol='orange',add=TRUE)


#pos9

EF<-Pdat[Pdat$model == 'neutRdiffGneutC' & Pdat$dispersal == 'uni' & Pdat$ddstrength == '_1_5',]
EFnums<-as.vector(rep(0,720))
for(i in 1:20){
  EFnums<-cbind(EFnums,as.numeric(as.character(EF[,i+1])))
}

AMC<-Pdat[Pdat$model == 'neutRneutGdiffC'& Pdat$dispersal == 'uni' & Pdat$ddstrength == '_1_5',]
AMCnums<-as.vector(rep(0,720))
for(i in 1:20){
  AMCnums<-cbind(AMCnums,as.numeric(as.character(AMC[,i+1])))
}

N<-Pdat[Pdat$model == 'neutRneutGneutC'& Pdat$dispersal == 'uni' & Pdat$ddstrength == '_1_5',]
Nnums<-as.vector(rep(0,720))
for(i in 1:20){
  Nnums<-cbind(Nnums,as.numeric(as.character(N[,i+1])))
}
boxplot(EFnums,col=alpha(rgb(215/255,48/255,31/255),.5),cex.axis=1.1,cex=.5,
        whisklty=1,names=as.factor(rad),ylim=c(0,1),ylab=NA,xlab=NA,
        staplecol='red',medcol='red',outpch=16,outcol='red',
        whiskcol='red',yaxt='n')
boxplot(Nnums,col=alpha('black',.5),whisklty=1,cex=.5,
        names=as.factor(rad),axes=FALSE,staplecol='black',
        medcol='black',outpch=16,outcol='black',whiskcol='black',
        add=TRUE)
boxplot(AMCnums,col=alpha(rgb(253/255,204/255,138/255),.5),whisklty=1,cex=.5,
        names=as.factor(rad),axes=FALSE,staplecol='orange',
        medcol='orange',outpch=16,outcol='orange',whiskcol='orange',add=TRUE)



#S####
library(ggplot2)
#bysps<-c('annual-edge','perennial-int','annual-int','perennial-edge')
#for(sp in 1:length(bysps)){#:length(bysps)){
#get abiotic predictors
setwd("~/Documents/MATLAB/simulation_2018-2019")
smallpred<-read.csv('smallpred.csv',header = TRUE,sep=',')
largepred<-read.csv('largepred.csv',header = TRUE,sep = ',')
#get habitat
setwd("~/Documents/MATLAB/simulation/landscapes")
bp<-read.csv('breakpoints.csv',header=FALSE,sep=',')
np<-read.csv('np.csv',header=FALSE,sep=',')

model<-c('neutRdiffGneutC','neutRneutGdiffC','neutRdiffGdiffC','neutRdiffG#diffC')
ddstrength<-c('_0_5','_1_0','_1_5')#,'alpha = _9','alpha = 1_1')#
dispersal<-c('adj','int','uni')
colors<-c('blue','red','purple','orange')
#cont<-c('epsilonE','epsilonC','epsilonEC','epsilonE#C')
mod<-matrix(c(1:16),nrow=4,ncol=4)
totalcontrib<-1

all_micro<-list()
all_small<-list()
all_large<-list()
g_pos_micro<-list()
g_neg_micro<-list()
g_pos_small<-list()
g_neg_small<-list()
g_pos_large<-list()
g_neg_large<-list()
g_all_micro<-list()
g_all_small<-list()
g_all_large<-list()
for(d in 1:length(dispersal)){
  for(s in 1:length(ddstrength)){#length(alpha)){
    
    setwd(paste('~/Documents/submissions/NichePart2020/uncorrelatedforspatialstats_editedsumm2019/neutral/',dispersal[d],'/',ddstrength[s],'/','neutRneutGneutC',sep=''))
    #setwd(paste('~/Documents/MATLAB/simulation2019Jan/neutral/',dispersal[d],'/',ddstrength[s],'/','neutRneutGneutC',sep=''))
    results_micro<-read.csv(paste('results_micro.csv',sep=''),header = TRUE, sep=',')
    results_small<-read.csv(paste('results_small.csv',sep=''), header = TRUE, sep = ',')
    results_large<-read.csv(paste('results_large.csv',sep=''), header = TRUE, sep = ',')
    #results_small<-results_small[results_small$landscape != 7 & results_small$landscape != 14 & results_small$landscape != 16,]
    #results_large<-results_large[results_large$landscape != 7 & results_large$landscape != 14 & results_large$landscape != 16,]
    results_small<-cbind(smallpred,results_small)
    results_large<-cbind(largepred,results_large)
    
    #for(l in unique(results_small$landscape)){
    #  results_small$hab[results_small$landscape==l]<-unlist(lapply(results_small$meanC[results_small$landscape==l],function(x) 
    #    if(x<=bp[l,1]){'annual-edge'}else if(x>bp[l,1]&x<=bp[l,2])
    #    {'perennial-int'}else if(x>bp[l,2]&x<=bp[l,3]){'annual-int'}
    #    else{'perennial-edge'}))}
    
    #results_small<-results_small[results_small$hab==bysps[sp],]
    v_data_micro<-data.frame(results_micro$landscape,results_micro$S)
    colnames(v_data_micro)<-c('landscape','v0')
    
    v_data_small<-data.frame(results_small$landscape,results_small$S)
    colnames(v_data_small)<-c('landscape','v0')
    
    v_data_large<-data.frame(results_large$landscape,results_large$S)
    colnames(v_data_large)<-c('landscape','v0')
    
    #plot(v_data_small$v0~results_small$varC,xlim=c(0,5),ylim=c(-4,4),cex=.5,pch=16,
    #col='black',main=paste('small subunits,',ddstrength[s],',',dispersal[d],sep=' '))
    
    micro_data<-matrix(NA,nrow=dim(v_data_micro)[1],ncol=4,byrow=FALSE)
    small_data<-matrix(NA,nrow=dim(v_data_small)[1],ncol=4,byrow=FALSE)
    large_data<-matrix(NA,nrow=dim(v_data_large)[1],ncol=4,byrow=FALSE)
    
    for(m in 1:length(model)){
      setwd(paste('~/Documents/submissions/NichePart2020/uncorrelatedforspatialstats_editedsumm2019/neutral/',dispersal[d],'/',ddstrength[s],'/',model[m],sep=''))
      results_micro<-read.csv(paste('results_micro.csv',sep=''),header = TRUE, sep = ',')
      results_small<-read.csv(paste('results_small.csv',sep=''), header = TRUE, sep = ',')
      results_large<-read.csv(paste('results_large.csv',sep=''), header = TRUE, sep = ',')
      #results_small<-results_small[results_small$landscape != 7 & results_small$landscape != 14 & results_small$landscape != 16,]
      #results_large<-results_large[results_large$landscape != 7 & results_large$landscape != 14 & results_large$landscape != 16,]
      #results_small<-cbind(smallpred,results_small)
      #results_large<-cbind(largepred,results_large)
      
      #for(l in unique(results_small$landscape)){
      #  results_small$hab[results_small$landscape==l]<-unlist(lapply(results_small$meanC[results_small$landscape==l],function(x) 
      #    if(x<=bp[l,1]){'annual-edge'}else if(x>bp[l,1]&x<=bp[l,2])
      #    {'perennial-int'}else if(x>bp[l,2]&x<=bp[l,3]){'annual-int'}
      #    else{'perennial-edge'}))}
      
      #results_small<-results_small[results_small$hab==bysps[sp],]
      v_micro<-results_micro$S
      if (m == 1 | m == 2){
      epsilon_micro<-(v_micro - v_data_micro$v0)}else{
      epsilon_micro<-v_micro - (v_data_micro$v0 + micro_data[,1] + micro_data[,2])
      }
      
      v_small<-results_small$S
      if(m == 1 | m == 2){
      epsilon_small<-(v_small - v_data_small$v0)}else{
      epsilon_small<-v_small - (v_data_small$v0 + small_data[,1] + small_data[,2])  
      }
      #v_data_small$epsilonE<-(v_data_small$v0-v_data_small$vE)
      
      v_large<-results_large$S
      if(m == 1 | m == 2){
      epsilon_large<-(v_large-v_data_large$v0)}else{
        epsilon_large<-v_large - (v_data_large$v0 + large_data[,1] + large_data[,2])
      }
      #points(results_small$varC,epsilon,pch=16,cex=.5,col=colors[m])
      
      micro_data[,m]<-epsilon_micro
      small_data[,m]<-epsilon_small
      large_data[,m]<-epsilon_large
    }
    micro_corr<-micro_data[,3] - micro_data[,4]
    small_corr<-small_data[,3] - small_data[,4]
    large_corr<-large_data[,3] - large_data[,4]
    
    micro_data<-cbind(micro_data,micro_corr)
    small_data<-cbind(small_data,small_corr)
    large_data<-cbind(large_data,large_corr)
    
    colnames(micro_data)<-c('E','C','ECtotal','ECind','ECcorr')
    colnames(small_data)<-c('E','C','ECtotal','ECind','ECcorr')
    colnames(large_data)<-c('E','C','ECtotal','ECind','ECcorr')
    
    
    all_micro[[mod[d,s]]]<-data.frame(v_data_micro,micro_data)
    all_small[[mod[d,s]]]<-data.frame(v_data_small,small_data)
    all_large[[mod[d,s]]]<-data.frame(v_data_large,large_data)


    if(totalcontrib == 1){
      if(abscontrib == 1){
        allmeans_micro<-apply(all_micro[[mod[d,s]]][,c(2:4,6:7)],2,function(x) mean(abs(x)))
        allmeans_small<-apply(all_small[[mod[d,s]]][,c(2:4,6:7)],2,function(x) mean(abs(x)))
        allmeans_large<-apply(all_large[[mod[d,s]]][,c(2:4,6:7)],2,function(x) mean(abs(x)))}else{
          allmeans_micro<-apply(all_micro[[mod[d,s]]][,c(2:4,6:7)],2,function(x) mean((x)))
          allmeans_small<-apply(all_small[[mod[d,s]]][,c(2:4,6:7)],2,function(x) mean((x)))
          allmeans_large<-apply(all_large[[mod[d,s]]][,c(2:4,6:7)],2,function(x) mean((x)))}}else{
            #SPLIT INTO NET ABOVE 0 AND NET BELOW 0
            if(abscontrib == 1){
              negmean_micro<-apply(all_micro[[mod[d,s]]][,c(2:4,6:7)],2,function(x) mean(abs(x[x<0])))
              posmean_micro<-apply(all_micro[[mod[d,s]]][,c(2:4,6:7)],2,function(x) mean(abs(x[x>0])))
              
              negmean_small<-apply(all_small[[mod[d,s]]][,c(2:4,6:7)],2,function(x) mean(abs(x[x<0])))
              posmean_small<-apply(all_small[[mod[d,s]]][,c(2:4,6:7)],2,function(x) mean(abs(x[x>0])))
              
              negmean_large<-apply(all_large[[mod[d,s]]][,c(2:4,6:7)],2,function(x) mean(abs(x[x<0])))
              posmean_large<-apply(all_large[[mod[d,s]]][,c(2:4,6:7)],2,function(x) mean(abs(x[x>0])))}else{
                
                negmean_micro<-apply(all_micro[[mod[d,s]]][,c(2:4,6:7)],2,function(x) mean((x[x<0])))
                posmean_micro<-apply(all_micro[[mod[d,s]]][,c(2:4,6:7)],2,function(x) mean((x[x>0])))
                
                negmean_small<-apply(all_small[[mod[d,s]]][,c(2:4,6:7)],2,function(x) mean((x[x<0])))
                posmean_small<-apply(all_small[[mod[d,s]]][,c(2:4,6:7)],2,function(x) mean((x[x>0])))
                
                negmean_large<-apply(all_large[[mod[d,s]]][,c(2:4,6:7)],2,function(x) mean((x[x<0])))
                posmean_large<-apply(all_large[[mod[d,s]]][,c(2:4,6:7)],2,function(x) mean((x[x>0])))} 
          }
    
    if(totalcontrib == 1){
      allmeans_micro<-as.data.frame(allmeans_micro)
      allmeans_small<-as.data.frame(allmeans_small)
      allmeans_large<-as.data.frame(allmeans_large)}else{
        negmeans_micro<-as.data.frame(negmean_micro)
        posmeans_micro<-as.data.frame(posmean_micro)
        negmeans_small<-as.data.frame(negmean_small)
        posmeans_small<-as.data.frame(posmean_small)
        negmeans_large<-as.data.frame(negmean_large)
        posmeans_large<-as.data.frame(posmean_large)
      }
    
    if(totalcontrib==1){
      group<-rep(paste(dispersal[d],ddstrength[s],sep=''),5)
      type<-rbind('0','E','C','ECind','ECcorr')
      g_all_micro[[mod[d,s]]]<-as.data.frame(cbind(group,type,allmeans_micro))
      g_all_small[[mod[d,s]]]<-as.data.frame(cbind(group,type,allmeans_small))
      g_all_large[[mod[d,s]]]<-as.data.frame(cbind(group,type,allmeans_large))}else{
        group<-rep(paste(dispersal[d],ddstrength[s],sep=''),5)
        type<-rbind('0','E','C','E_C','EC')
        g_neg_micro[[mod[d,s]]]<-as.data.frame(cbind(group,type,negmeans_micro))
        g_pos_micro[[mod[d,s]]]<-as.data.frame(cbind(group,type,posmeans_micro))
        
        g_neg_small[[mod[d,s]]]<-as.data.frame(cbind(group,type,negmeans_small))
        g_pos_small[[mod[d,s]]]<-as.data.frame(cbind(group,type,posmeans_small))
        
        g_neg_large[[mod[d,s]]]<-as.data.frame(cbind(group,type,negmeans_large))
        g_pos_large[[mod[d,s]]]<-as.data.frame(cbind(group,type,posmeans_large))
      }
    
  }
}


if(totalcontrib==1){
  g_all_micro<-do.call(rbind,g_all_micro)
  #plot with neutral values?
  #micro<-ggplot(data = g_all_micro, aes(x = group, y = allmeans_micro, fill = factor(type,levels=c("ECcorr", "ECind", "C", "E","0")))) + 
  #  geom_bar(stat = "identity")+scale_fill_manual("legend", values = c("0" = "black", "E" = "blue", "C" = "red", "ECind" = "orange", "ECcorr" = "purple"))+ggtitle('microsites') + ylim(-5,5)+ scale_x_discrete(limits=c('adj_0_5','adj_1_0','adj_1_5','int_0_5','int_1_0','int_1_5','uni_0_5','uni_1_0','uni_1_5'))
  #print(micro)

    g_all_micro_alt<-g_all_micro[g_all_micro$type!='0',]
  micro<-ggplot(data = g_all_micro_alt, aes(x = group, y = allmeans_micro, fill = factor(type,levels=c("ECcorr", "ECind", "C", "E")))) + 
    geom_bar(stat = "identity")+scale_fill_manual("legend", values = c("E" = "blue", "C" = "red", "ECind" = "orange", "ECcorr" = "purple"))+ggtitle('microsites') + ylim(-2.25,1)+ scale_x_discrete(limits=c('adj_0_5','adj_1_0','adj_1_5','int_0_5','int_1_0','int_1_5','uni_0_5','uni_1_0','uni_1_5'))+ylab('deltaS')
  print(micro)
  
  
  g_all_small<-do.call(rbind,g_all_small)
  #small<-ggplot(data = g_all_small, aes(x = group, y = allmeans_small, fill = factor(type,levels=c("ECcorr","ECind", "C", "E","0")))) + 
  #  geom_bar(stat = "identity")+scale_fill_manual("legend", values = c("0" = "black", "E" = "blue", "C" = "red","ECind" = "orange","ECcorr" = "purple"))+ggtitle('small subunits') + ylim(-5,5)+ scale_x_discrete(limits=c('adj_0_5','adj_1_0','adj_1_5','int_0_5','int_1_0','int_1_5','uni_0_5','uni_1_0','uni_1_5'))
  #print(small)
  
  g_all_small_alt<-g_all_small[g_all_small$type!='0',]
  small<-ggplot(data = g_all_small_alt, aes(x = group, y = allmeans_small, fill = factor(type,levels=c("ECcorr", "ECind", "C", "E")))) + 
    geom_bar(stat = "identity")+scale_fill_manual("legend", values = c("E" = "blue", "C" = "red", "ECind" = "orange", "ECcorr" = "purple"))+ggtitle('small subunits') + ylim(-2.25,1)+ scale_x_discrete(limits=c('adj_0_5','adj_1_0','adj_1_5','int_0_5','int_1_0','int_1_5','uni_0_5','uni_1_0','uni_1_5'))
  print(small)
  
  
  g_all_large<-do.call(rbind,g_all_large)
  #large<-ggplot(data = g_all_large, aes(x = group, y = allmeans_large, fill = factor(type,levels=c("ECcorr","ECind", "C", "E","0")))) + 
  #  geom_bar(stat = "identity")+scale_fill_manual("legend", values = c("0" = "black", "E" = "blue", "C" = "red","ECind" = "orange","ECcorr" = "purple"))+ggtitle('large subunits') + ylim(-5,5)+ scale_x_discrete(limits=c('adj_0_5','adj_1_0','adj_1_5','int_0_5','int_1_0','int_1_5','uni_0_5','uni_1_0','uni_1_5'))
  #print(large)
  
  g_all_large_alt<-g_all_large[g_all_large$type!='0',]
  large<-ggplot(data = g_all_large_alt, aes(x = group, y = allmeans_large, fill = factor(type,levels=c("ECcorr", "ECind", "C", "E")))) + 
    geom_bar(stat = "identity")+scale_fill_manual("legend", values = c("E" = "blue", "C" = "red", "ECind" = "orange", "ECcorr" = "purple"))+ggtitle('large subunits') + ylim(-2.25,1)+ scale_x_discrete(limits=c('adj_0_5','adj_1_0','adj_1_5','int_0_5','int_1_0','int_1_5','uni_0_5','uni_1_0','uni_1_5'))
  print(large)
  
  
  }else{
    g_neg_micro<-na.omit(do.call(rbind,(as.matrix(g_neg_micro))))
    g_pos_micro<-do.call(rbind,(as.matrix(g_pos_micro)))
    
    g_all_large_alt<-g_all_large[g_all_large$type!='0',]
    large<-ggplot(data = g_all_large_alt, aes(x = group, y = allmeans_large, fill = factor(type,levels=c("ECcorr", "ECind", "C", "E")))) + 
      geom_bar(stat = "identity")+scale_fill_manual("legend", values = c("E" = "blue", "C" = "red", "ECind" = "orange", "ECcorr" = "purple"))+ggtitle('large subunits') + ylim(-5,5)+ scale_x_discrete(limits=c('adj_0_5','adj_1_0','adj_1_5','int_0_5','int_1_0','int_1_5','uni_0_5','uni_1_0','uni_1_5'))
    print(large)
    
    netgain_micro<-ggplot(data = g_pos_micro, aes(x = group, y = posmean_micro, fill = factor(type,levels=c("EC","E_C", "C", "E","0")))) + 
      geom_bar(stat = "identity")+scale_fill_manual("legend", values = c("0" = "black", "E" = "blue", "C" = "red","E_C" = "orange","EC" = "purple"))+ggtitle('microsites')+ ylim(-5,5)+ scale_x_discrete(limits=c('adj_0_5','adj_1_0','adj_1_5','int_0_5','int_1_0','int_1_5','uni_0_5','uni_1_0','uni_1_5'))
    netloss_micro<-layer(geom='bar',data = g_neg_micro, aes(x = group, y = negmean_micro, fill = factor(type,levels=c("EC", "E_C", "C", "E","0"))), stat = "identity",show.legend=FALSE,position=position_stack(vjust=0))
    #net_small<-netgain_small+netloss_small
    print(netgain_micro+netloss_micro)
    
    g_neg_small<-na.omit(do.call(rbind,(as.matrix(g_neg_small))))
    g_pos_small<-do.call(rbind,(as.matrix(g_pos_small)))
    
    netgain_small<-ggplot(data = g_pos_small, aes(x = group, y = posmean_small, fill = factor(type,levels=c("EC","E_C", "C", "E","0")))) + 
      geom_bar(stat = "identity")+scale_fill_manual("legend", values = c("N" = "black", "E" = "blue", "C" = "red","E_C" = "orange","EC" = "purple"))+ggtitle('small subunits')+ ylim(-5,5)+ scale_x_discrete(limits=c('adj_0_5','adj_1_0','adj_1_5','int_0_5','int_1_0','int_1_5','uni_0_5','uni_1_0','uni_1_5'))
    netloss_small<-layer(geom='bar',data = g_neg_small, aes(x = group, y = negmean_small, fill = factor(type,levels=c("EC","E_C", "C", "E","0"))), stat = "identity",show.legend=FALSE,position=position_stack(vjust=0))
    #net_small<-netgain_small+netloss_small
    print(netgain_small+netloss_small)
    
    g_neg_large<-na.omit(do.call(rbind,(as.matrix(g_neg_large))))
    g_pos_large<-do.call(rbind,(as.matrix(g_pos_large)))
    
    netgain_large<-ggplot(data = g_pos_large, aes(x = group, y = posmean_large, fill = factor(type,levels=c("EC","E_C", "C", "E","0")))) + 
      geom_bar(stat = "identity")+scale_fill_manual("legend", values = c("0" = "black", "E" = "blue", "C" = "red","E_C" = "orange"),"EC" = "purple")+ggtitle('large subunits') + ylim(-5,5)+ scale_x_discrete(limits=c('adj_0_5','adj_1_0','adj_1_5','int_0_5','int_1_0','int_1_5','uni_0_5','uni_1_0','uni_1_5'))
    netloss_large<-layer(geom='bar',data = g_neg_large, aes(x = group, y = negmean_large, fill = factor(type,levels=c("EC","E_C", "C", "E","0"))), stat = "identity",show.legend=FALSE,position=position_stack(vjust=0))
    net_large<-netgain_large+netloss_large
    print(net_large)
  }



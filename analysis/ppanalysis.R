library(spatstat)
library(MuMIn)
library(raster)
library(geoR)
library(pscl)
library(gstat)
library(sp)
library(ggplot2)

#function that calculates the ISAR and Pc at the same freakin time:

AllFunctions <- function(distances, interior, radius) {
  
  focalsps<-unique(rownames(distances))
  allsps<-unique(colnames(distances))
  
  Pc_by_r<-NULL
  ISAR_by_r<-NULL
  for (r in 1:length(radius)) {
    Pc_by_f <- NULL
    ISAR_by_f<-NULL
    for(f in 1:length(focalsps)){
      
      
      allfocalinds<-rownames(distances)==as.character(focalsps[f])
      inds<-which(allfocalinds)
      
      ##DO SOME Pc CALCS##
      dist1<-as.matrix(distances[inds,]) #this should be the same as focalinds in the ISAR code.
      
      if (dim(dist1)[2] == 1) { #if it's a vector... in the form of a matrix
        
        d1subset <- dist1[which(rownames(dist1)==as.character(focalsps[f]))]
        consp <- length(which(d1subset[] <= radius[r] & d1subset[] > 0))
        totaln <- length(which(dist1[] <= radius[r] & dist1[] > 0))
        if (totaln == 0) {
          Pc_by_f[f] <- 0
        }
        else {
          Pc_by_f[f] <- consp / totaln
        }
        
      }
      else {
        pc <- NULL
        for (i in 1:dim(dist1)[1]) {
          d1subset <- dist1[i,which(colnames(dist1)==as.character(focalsps[f]))]
          consp <- length(which(d1subset[] <= radius[r] & d1subset[] > 0))
          totaln <- length(which(dist1[i,] <= radius[r] & dist1[i,] > 0))
          if (totaln == 0) {
            pc[i] <- 0
          }
          else {
            pc[i] <- consp / totaln 
          }
        }
        Pc_by_f[f] <- sum(pc) / dim(dist1)[1]
      }
      
      ##DO SOME ISAR CALCS##
      
      focalinds<-distances[inds,] #this is a matrix of focal individuals as rows
      if(length(inds)>1){
        indexcols<-colnames(focalinds)
        indexrows<-rep(as.character(focalsps[f]),dim(focalinds)[1])}else{
          indexcols<-names(focalinds)
          indexrows<-as.character(focalsps[f])
        }
      if(is.vector(focalinds)==FALSE){
        ns<-NULL
        for(i in 1:dim(focalinds)[1]){
          subinds<-focalinds[i,]<=radius[r]
          tocalc<-indexcols[which(subinds)]
          #including the focal species:
          ns[i]<-length(unique(tocalc))
          #otherwise
          #ns[i]<-length(unique(tocalc))-1
        }
        allns<-length(allsps)
        #if you want to normalize:
        ISAR_by_f[f]<-mean(ns)/allns}
      #otherwise:
      #by.f[f]<-mean(ns)}
      else{ns<-NULL
      subinds<-focalinds<=radius[r]
      tocalc<-indexcols[which(subinds)]
      #including the focal individual:
      ns<-length(unique(tocalc))  #otherwise
      #ns<-length(unique(tocalc))-1
      allns<-length(unique(colnames(distances)))
      #if you want to normalize:
      ISAR_by_f[f]<-mean(ns)/allns}
      #otherwise:
      #by.f[f]<-mean(ns)}
    }
    Pc_by_r[r] <- sum(Pc_by_f) / length(Pc_by_f)
    ISAR_by_r[r]<-sum(ISAR_by_f)/allns
  }
  # par(mar=c(5,4,4,2))
  # pcplot<-data.frame(radius,by_r)
  # plot(by_r~radius,data=pcplot,pch=16,xlim=c(0,1),ylim=c(0,length(unique(dataf$IDf))),
  #      ylab='Proportion Conspecifics (community average)',main='Pc')
  
  output<-list(Pc_by_r,ISAR_by_r)
  return(output)
}

#OLD###########
# PcFunction <- function(distances, interior, radius) {
#   
#   
#   by_r<-NULL
#   for (r in 1:length(radius)) {
#     by_f <- NULL
#     spslist <- unique(interior$ID)
#     for (spsn in 1:length(spslist)) {
#       
#       spstype <- as.character(spslist[spsn])
#       allfocalinds<-rownames(distances)==spstype
#       inds<-which(allfocalinds)
#       dist1<-as.matrix(distances[inds,])
#       
#       if (dim(dist1)[2] == 1) { #if it's a vector... in the form of a matrix
#         
#         d1subset <- dist1[which(rownames(dist1)==spstype)]
#         consp <- length(which(d1subset[] <= radius[r] & d1subset[] > 0))
#         totaln <- length(which(dist1[] <= radius[r] & dist1[] > 0))
#         if (totaln == 0) {
#           by_f[spsn] <- 0
#         }
#         else {
#           by_f[spsn] <- consp / totaln
#         }
#         
#       }
#       else {
#         pc <- NULL
#         for (i in 1:dim(dist1)[1]) {
#           d1subset <- dist1[i,which(colnames(dist1)==spstype)]
#           consp <- length(which(d1subset[] <= radius[r] & d1subset[] > 0))
#           totaln <- length(which(dist1[i,] <= radius[r] & dist1[i,] > 0))
#           if (totaln == 0) {
#             pc[i] <- 0
#           }
#           else {
#             pc[i] <- consp / totaln 
#           }
#         }
#         by_f[spsn] <- sum(pc) / dim(dist1)[1]
#       }
#       
#     }
#     by_r[r] <- sum(by_f) / length(by_f)
#   }
#   # par(mar=c(5,4,4,2))
#   # pcplot<-data.frame(radius,by_r)
#   # plot(by_r~radius,data=pcplot,pch=16,xlim=c(0,1),ylim=c(0,length(unique(dataf$IDf))),
#   #      ylab='Proportion Conspecifics (community average)',main='Pc')
#   return(by_r)
# }
# 
# ISARcomm<-function(data,distances,radius,print=FALSE){
#   focalsps<-unique(rownames(distances))
#   allsps<-unique(colnames(distances))
#   by.r<-NULL
#   for(r in 1:length(radius)){
#     by.f<-NULL
#     for(f in 1:length(focalsps)){
#       allfocalinds<-rownames(distances)==as.character(focalsps[f])
#       inds<-which(allfocalinds)
#       focalinds<-distances[inds,] #this is a matrix of focal individuals as rows
#       if(length(inds)>1){
#         indexcols<-colnames(focalinds)
#         indexrows<-rep(as.character(focalsps[f]),dim(focalinds)[1])}else{
#           indexcols<-names(focalinds)
#           indexrows<-as.character(focalsps[f])
#         }
#       if(is.vector(focalinds)==FALSE){
#         ns<-NULL
#         for(i in 1:dim(focalinds)[1]){
#           subinds<-focalinds[i,]<=radius[r]
#           tocalc<-indexcols[which(subinds)]
#           #including the focal species:
#           ns[i]<-length(unique(tocalc))
#           #otherwise
#           #ns[i]<-length(unique(tocalc))-1
#         }
#         allns<-length(allsps)
#         #if you want to normalize:
#         by.f[f]<-mean(ns)/allns}
#       #otherwise:
#       #by.f[f]<-mean(ns)}
#       else{ns<-NULL
#       subinds<-focalinds<=radius[r]
#       tocalc<-indexcols[which(subinds)]
#       #including the focal individual:
#       ns<-length(unique(tocalc))  #otherwise
#       #ns<-length(unique(tocalc))-1
#       allns<-length(unique(colnames(distances)))
#       #if you want to normalize:
#       by.f[f]<-mean(ns)/allns}
#       #otherwise:
#       #by.f[f]<-mean(ns)}
#     }
#     by.r[r]<-sum(by.f)/allns
#   }
#   return(by.r)
# }
#TOUSE########################################################################
numls<-20
dispersal<-c('adj','int','uni')
ddstrength<-c('_0_5','_1_0','_1_5')
models<-c('neutRdiffGdiffC','neutRdiffGneutC','neutRneutGdiffC','neutRneutGneutC','neutRdiffG#diffC')
ISAR_outputmatrix<-matrix(data=NA,ncol=25,nrow=900,byrow=TRUE)
colnames(ISAR_outputmatrix)<-c('landscape','r1','r2','r3','r4','r5','r6','r7','r8','r9','r10','r11','r12','r13','r14','r15','r16','r17','r18','r19','r20','r21','dispersal','ddstrength','model')
Pc_outputmatrix<-matrix(data=NA,ncol=25,nrow=900,byrow=TRUE)
colnames(Pc_outputmatrix)<-c('landscape','r1','r2','r3','r4','r5','r6','r7','r8','r9','r10','r11','r12','r13','r14','r15','r16','r17','r18','r19','r20','r21','dispersal','ddstrength','model')
run<-1
#names<-c('combomap','efmap','amcmap','neutmap')
#colvec<-c('lightblue','black','blue','slategrey')
#tomap<-list()
#whichsim<-list()

for(m in 1:length(models)){
#m=4
for(d in 1:length(dispersal)){for(s in 1:length(ddstrength)){
#d=2
#s=2
  for(ls in 1:numls){
    
    setwd(paste('~/Documents/submissions/NichePart2020/for Dryad/2020/code_editedsumm2019/neutral/',dispersal[d],'/',ddstrength[s],'/',models[m],sep=''))
    tomap<-read.csv(paste('ppmap_ls',ls,'.csv',sep=''),header = TRUE,sep=',')
    le<-10
    ue<-65-10
    interior<-tomap[c(tomap$X>=le & tomap$Y>=le & tomap$X<=ue & 
                        tomap$Y<=ue),]
    distmat<-round(crossdist.default(X=interior$X,Y=interior$Y,x2=tomap$X,y2=tomap$Y),2)
    colnames(distmat)<-c(as.character(tomap$ID))
    rownames(distmat)<-c(as.character(interior$ID))
    r<-seq(0,10,.5)
    output<-AllFunctions(distmat,interior,radius=r)
    Pc_output<-output[[1]]
    ISAR_output<-output[[2]]
    
    # ISAR_output<-ISARcomm(data=tomap,distances=distmat,radius=r)
    # radius<-r
    # Pc_output<-PcFunction(distmat,interior,radius=r)
    ISAR_outputmatrix[run,]<-c(ls,t(ISAR_output),dispersal[d],ddstrength[s],models[m])
    Pc_outputmatrix[run,]<-c(ls,t(Pc_output),dispersal[d],ddstrength[s],models[m])
    run<-run+1
}}}}


#OLD#############
# setwd(paste('~/Documents/MATLAB/simulation2019Jan/neutral/',dispersal[d],'/',ddstrength[s],'/',models[m],sep=''))
# tomap<-read.csv(paste('ppmap_ls',ls,'.csv',sep=''),header = TRUE,sep=',')
# le<-10
# ue<-65-10
# interior<-tomap[c(tomap$X>=le & tomap$Y>=le & tomap$X<=ue & 
#                     tomap$Y<=ue),]
# distmat<-round(crossdist.default(X=interior$X,Y=interior$Y,x2=tomap$X,y2=tomap$Y),2)
# colnames(distmat)<-c(as.character(tomap$ID))
# rownames(distmat)<-c(as.character(interior$ID))
# 
# r<-seq(0,10,.5)
# #output<-ISARcomm(data=tomap,distances=distmat,radius=r)
# #points(r,output,col='purple',type = 'l',lwd=2)
# 
# 
# 
# radius<-r
# Pcout<-PcFunction(distmat,interior,radius=r)
# points(Pcout~radius,col='black',type='l',lwd=2)
# 

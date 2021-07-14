#N v. I with reduced abiotic heterogeneity:

#pull richness data from the first 5 landscapes
#dispersal = intermediate
#alpha_ij = 1

#from the first initial condition in the regular 
#STD=2.5 landscapes

Nreg<-c(3,4,4,4,4,4,3,4,3,4,4,4,4,3,4,4,4,4,4,4)
Ireg<-c(3,4,4,4,4,4,3,4,3,4,4,4,4,3,4,4,4,4,4,4)

#from the same landscapes and initial conditions with 
#STD=.5
N<-c(3,4,3,4,4,4,3,4,3,4,4,4,4,3,4,4,4,4,4,4)
I<-c(3,2,3,3,3,3,3,3,3,2,3,3,3,3,3,3,2,3,2,3)

#from the same landscapes and initial conditions with 
#STD=.025
Nlow<-c(3,4,3,4,4)
Ilow<-c(2,2,2,2,2)

STD<-c(rep('2.5',40),rep('0.5',40))#,rep('0.025',10))
STD<-factor(STD,levels=c('2.5','0.5'))#,'0.025'))
model<-c(rep('N',20),rep('I',20),rep('N',20),
         rep('I',20))#,rep('N',5),rep('I',5))
model<-factor(model,levels=c('N','I'))
values<-c(Nreg,Ireg,N,I)#,Nlow,Ilow)

dat<-data.frame(model,STD,values)

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


ggplot(dat, aes(x=model, y=values, color=STD,fill=STD)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size=12),legend.background=element_blank())+scale_fill_manual(values=cbbPalette)+scale_color_manual(values=cbbPalette)+
  geom_boxplot()

mod<-lm(values~model*STD,data=dat)

Anova(mod,type=3)

set_theme(legend.item.backcol = 'white',
          axis.textcolor.x='black',axis.textcolor.y='black')

plot_model(mod,type='emm',terms=c('model','STD'),
           col=c(rgb(136/255,86/255,167/255),rgb(158/255,188/255,218/255)))+
  scale_x_discrete("type", expand=c(0,.5))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "black"),axis.text=element_text(size=12))+
  ylim(0,4)
 
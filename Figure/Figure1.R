library(tidyverse)
library(ggpubr)

changeP <- function(p_val,r2) {
  tmp=strsplit(as.character(p_val),"e")[[1]]
  if(length(tmp)>1){
    a=tmp[1];b=tmp[2]
    as.expression(bquote(atop(p~"="~.(a)~"×"~10^~.(b)~","~R^2~"="~.(r2))))
  }else{
    as.expression(bquote(atop(p~"="~.(p_val)~","~R^2~"="~.(r2))))
  }
}

plot_physical_ham=function(data,patient_paper){
  s = summary(lm(`Genetic distance`~`Physical distance`, data=data))
  r2 =signif(s$r.squared,3)
  p_val=signif(s$coefficients[8],3)
  
  p=ggplot(data = data,aes(x=`Physical distance`,
                           y=`Genetic distance`))+
    geom_point(color="lightblue")+
    theme_classic()+
    ggtitle(patient_paper)+xlab("Physical distance")+
    ylab("Genetic distance")+theme_classic()+
    labs(subtitle = changeP(p_val,r2))+
    theme(legend.position = "none",plot.title = element_text(hjust = 0.5))+
    geom_smooth(method='lm',se=T,colour="red",linetype="dashed",
                linewidth=1,fill="pink")+
    theme(text=element_text(size = 6),
          legend.spacing = unit(0,"cm"),
          plot.title=element_text(hjust = 0.5,size=6.5),
          plot.subtitle = element_text(hjust = 0.5,size=6),
          legend.key.width = unit(0.1, 'cm'),
          legend.key.height = unit(0.3, 'cm'),
          line = element_line(size = 0.1),
          legend.box.spacing = unit(0.01,"cm"),
          axis.ticks.length = unit(0.03,"cm"))
  return(p)
}
plot_physical_rna=function(data,patient_paper){
  s = summary(lm(`Transcriptomic distance`~`Physical distance`, data=data))
  r2 =signif(s$r.squared,3)
  p_val=signif(s$coefficients[8],3)
  
  p=ggplot(data = data,aes(x=`Physical distance`,
                           y=`Transcriptomic distance`))+
    geom_point(color="lightblue")+
    theme_classic()+
    ggtitle(patient_paper)+xlab("Physical distance")+
    ylab("Transcriptomic distance")+theme_classic()+
    labs(subtitle = changeP(p_val,r2))+
    theme(legend.position = "none",plot.title = element_text(hjust = 0.5))+
    geom_smooth(method='lm',se=T,colour="red",linetype="dashed",
                linewidth=1,fill="pink")+
    theme(text=element_text(size = 6),
          legend.spacing = unit(0,"cm"),
          plot.title=element_text(hjust = 0.5,size=6.5),
          plot.subtitle = element_text(hjust = 0.5,size=6),
          legend.key.width = unit(0.1, 'cm'),
          legend.key.height = unit(0.3, 'cm'),
          line = element_line(size = 0.1),
          legend.box.spacing = unit(0.01,"cm"),
          axis.ticks.length = unit(0.03,"cm"))
  return(p)
}
plot_ham_rna=function(data,patient_paper){
  
  s = summary(lm(`Transcriptomic distance`~`Genetic distance`, data=data))
  r2 =signif(s$r.squared,3)
  p_val=signif(s$coefficients[8],3)
  
  p=ggplot(data = data,aes(x=`Genetic distance`,
                           y=`Transcriptomic distance`))+
    geom_point(color="lightblue")+
    theme_classic()+
    ggtitle(patient_paper)+xlab("Genetic distance")+
    ylab("Transcriptomic distance")+theme_classic()+
    labs(subtitle = changeP(p_val,r2))+
    theme(legend.position = "none",plot.title = element_text(hjust = 0.5))+
    geom_smooth(method='lm',se=T,colour="red",
                linetype="dashed",linewidth=1,fill="pink")+
    theme(text=element_text(size = 6),
          legend.spacing = unit(0,"cm"),
          plot.title=element_text(hjust = 0.5,size=6.5),
          plot.subtitle = element_text(hjust = 0.5,size=6),
          legend.key.width = unit(0.1, 'cm'),
          legend.key.height = unit(0.3, 'cm'),
          line = element_line(size = 0.1),
          legend.box.spacing = unit(0.01,"cm"),
          axis.ticks.length = unit(0.03,"cm"))
  return(p)
}

distance_table=read_tsv("../Data/Figure1.tsv")
pmulti1=ggplot();pmulti2=ggplot();pmulti3=ggplot()
color=c(RColorBrewer::brewer.pal(12,"Paired"),"darkolivegreen1")
for (patient in unique(distance_table$Patient)){
  data=distance_table%>%filter(Patient==patient)
  pmulti1=pmulti1+geom_smooth(data = data, 
                              aes(x=`Physical distance`,
                                  y=`Genetic distance`,
                                  color=Patient),
                              method='lm',se=FALSE,size=0.6)
  pmulti2=pmulti2+geom_smooth(data = data, 
                              aes(x=`Physical distance`,
                                  y=`Transcriptomic distance`,
                                  color=Patient),
                              method='lm',se=FALSE,size=0.6)
  pmulti3=pmulti3+geom_smooth(data = data, 
                              aes(x=`Genetic distance`,
                                  y=`Transcriptomic distance`,
                                  color=Patient),
                              method='lm',se=FALSE,size=0.6)
}
pmulti1=pmulti1+theme_classic()+
  xlab("Physical distance")+ylab("Genetic distance")+
  scale_color_manual(values=color)+ggtitle("All patients")+
  theme(legend.title = element_blank(),legend.position = "none",
        text=element_text(size = 6),
        legend.spacing = unit(0,"cm"),
        plot.title=element_text(hjust = 0.5,size=6.5),
        plot.subtitle = element_text(hjust = 0.5,size=6),
        legend.key.width = unit(0.1, 'cm'),
        legend.key.height = unit(0.3, 'cm'),
        line = element_line(size = 0.1),
        legend.box.spacing = unit(0.01,"cm"),
        axis.ticks.length = unit(0.03,"cm"))

pmulti2=pmulti2+theme_classic()+
  xlab("Physical distance")+
  ylab("Transcriptomic distance")+
  scale_color_manual(values=color)+ggtitle("All patients")+
  theme(legend.title = element_blank(),legend.position = "none",
        text=element_text(size = 6),
        legend.spacing = unit(0,"cm"),
        plot.title=element_text(hjust = 0.5,size=6.5),
        plot.subtitle = element_text(hjust = 0.5,size=6),
        legend.key.width = unit(0.1, 'cm'),
        legend.key.height = unit(0.3, 'cm'),
        line = element_line(size = 0.1),
        legend.box.spacing = unit(0.01,"cm"),
        axis.ticks.length = unit(0.03,"cm"))

pmulti3=pmulti3+theme_classic()+xlab("Genetic distance")+
  ylab("Transcriptomic distance")+
  scale_color_manual(values=color)+ggtitle("All patients")+
  theme(legend.title = element_blank(),
        text=element_text(size = 6),
        legend.spacing = unit(0,"cm"),
        plot.title=element_text(hjust = 0.5,size=6.5),
        plot.subtitle = element_text(hjust = 0.5,size=6),
        legend.key.width = unit(0.3, 'cm'),
        legend.key.height = unit(0.3, 'cm'),
        line = element_line(size = 0.1),
        legend.box.spacing = unit(0.01,"cm"),
        axis.ticks.length = unit(0.03,"cm"))



patient_paper="SH06"
data=distance_table%>%filter(Patient==patient_paper)
Fig1b1=plot_physical_ham(data,patient_paper)
Fig1b2=plot_physical_rna(data,patient_paper)
Fig1b3=plot_ham_rna(data,patient_paper)
Fig1b=ggarrange(Fig1b1,Fig1b2,Fig1b3,ncol = 3)
Fig1c=pmulti1;Fig1d=pmulti2;Fig1e=pmulti3

Fig1b

Fig1c

Fig1d

Fig1e

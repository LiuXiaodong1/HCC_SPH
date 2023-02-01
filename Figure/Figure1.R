library(tidyverse)

changeP <- function(p_val,r2) {
  tmp=strsplit(as.character(p_val),"e")[[1]]
  if(length(tmp)>1){
    a=tmp[1];b=tmp[2]
    as.expression(bquote(atop(p~"="~.(a)~"??"~10^~.(b)~","~R^2~"="~.(r2))))
  }else{
    as.expression(bquote(atop(p~"="~.(p_val)~","~R^2~"="~.(r2))))
  }
}


plot_physical_fst=function(patient,patient_paper){
  DNA=read.table(paste0(patient,"_fst.tsv"),header = T)
  RNA=read.table(paste0(patient,"_RNA_dis.tsv"),header = T)
  physical=read.table(paste0(patient,"_physical_dis.tsv"),header = T)
  
  data=cbind.data.frame(DNA,RNA=RNA$value,physical=physical$pos_dist)
  s = summary(lm(fst~physical, data=data))
  r2 =signif(s$r.squared,3)
  p_val=signif(s$coefficients[8],3)
  x=min(data$physical) + (max(data$physical)-min(data$physical))*0.25
  y=min(data$fst) + (max(data$fst)-min(data$fst))*0.75
  t=as.data.frame(cbind(r2,x,y,p_val))
  
  p=ggplot(data = data,aes(x=physical,y=fst))+geom_point(color="lightblue")+
    theme_classic()+
    ggtitle(patient_paper)+xlab("Physical distance")+
    ylab("FST")+theme_classic()+
    labs(subtitle = changeP(p_val,r2))+
    theme(legend.position = "none",plot.title = element_text(hjust = 0.5))+
    geom_smooth(method='lm',se=T,colour="red",linetype="dashed",size=1,fill="pink")+
    # annotate("text",parse = T,x=t$x,y=t$y,size=5,color="red",
    #          label=as.expression(bquote(atop(p~"="~.(t$p_val),
    #                                          R^2~"="~.(t$r2)))))+
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
all=read.table("E:/cancer genome/liver/????/distance/all.csv",header = T,sep = ",")
p1=NULL
for (patient in substr(list.files(pattern = "DT.*vcf"),1,4)){
  patient_paper=all[all$Patient==patient,]$Patient_paper[1]
  p1[[patient]]=plot_physical_fst(patient,patient_paper)
}

p1_all=ggpubr::ggarrange(plotlist = p1[c(-6,-9)],nrow = 3,ncol = 4)
ggsave("physical_fst.pdf",p1_all,width = 180,height = 150,units = "mm")

plot_physical_rna=function(patient,patient_paper){
  DNA=read.table(paste0(patient,"_fst.tsv"),header = T)
  RNA=read.table(paste0(patient,"_RNA_dis.tsv"),header = T)
  physical=read.table(paste0(patient,"_physical_dis.tsv"),header = T)
  
  data=cbind.data.frame(DNA,RNA=RNA$value,physical=physical$pos_dist)
  s = summary(lm(RNA~physical, data=data))
  r2 =signif(s$r.squared,3)
  p_val=signif(s$coefficients[8],3)
  x=min(data$physical) + (max(data$physical)-min(data$physical))*0.25
  y=min(data$RNA) + (max(data$RNA)-min(data$RNA))*0.75
  t=as.data.frame(cbind(r2,x,y,p_val))
  
  p=ggplot(data = data,aes(x=physical,y=RNA))+geom_point(color="lightblue")+
    theme_classic()+
    ggtitle(patient_paper)+xlab("Physical distance")+
    ylab("Transcriptomic distance")+theme_classic()+
    labs(subtitle = changeP(p_val,r2))+
    theme(legend.position = "none",plot.title = element_text(hjust = 0.5))+
    geom_smooth(method='lm',se=T,colour="red",linetype="dashed",size=1,fill="pink")+
    # annotate("text",parse = T,x=t$x,y=t$y,size=5,color="red",
    #          label=as.expression(bquote(atop(p~"="~.(t$p_val),
    #                                          R^2~"="~.(t$r2)))))+
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

p2=NULL
for (patient in substr(list.files(pattern = "DT.*vcf"),1,4)){
  patient_paper=all[all$Patient==patient,]$Patient_paper[1]
  p2[[patient]]=plot_physical_rna(patient,patient_paper)
}

p2_all=ggpubr::ggarrange(plotlist = p2[c(-6,-9)],nrow = 3,ncol = 4)
ggsave("physical_rna.pdf",p2_all,width = 180,height = 150,units = "mm")

plot_fst_rna=function(patient,patient_paper){
  DNA=read.table(paste0(patient,"_fst.tsv"),header = T)
  RNA=read.table(paste0(patient,"_RNA_dis.tsv"),header = T)
  physical=read.table(paste0(patient,"_physical_dis.tsv"),header = T)
  
  data=cbind.data.frame(DNA,RNA=RNA$value,physical=physical$pos_dist)
  s = summary(lm(RNA~fst, data=data))
  r2 =signif(s$r.squared,3)
  p_val=signif(s$coefficients[8],3)
  as.character(p_val)
  
  x=min(data$fst) + (max(data$fst)-min(data$fst))*0.25
  y=min(data$RNA) + (max(data$RNA)-min(data$RNA))*0.75
  t=as.data.frame(cbind(r2,x,y,p_val))
  
  p=ggplot(data = data,aes(x=fst,y=RNA))+geom_point(color="lightblue")+
    theme_classic()+
    ggtitle(patient_paper)+xlab("FST")+
    ylab("Transcriptomic distance")+theme_classic()+
    labs(subtitle = changeP(p_val,r2))+
    theme(legend.position = "none",plot.title = element_text(hjust = 0.5))+
    geom_smooth(method='lm',se=T,colour="red",linetype="dashed",size=1,fill="pink")+
    # annotate("text",parse = T,x=t$x,y=t$y,size=5,color="red",
    #          label=as.expression(bquote(atop(p~"="~.(t$p_val),
    #                                          R^2~"="~.(t$r2)))))+
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

p3=NULL
for (patient in sort(substr(list.files(pattern = "DT.*vcf"),1,4))){
  patient_paper=all[all$Patient==patient,]$Patient_paper[1]
  p3[[patient]]=plot_fst_rna(patient,patient_paper)
}

p3_all=ggpubr::ggarrange(plotlist = p3[c(-6,-9)],nrow = 3,ncol = 4)
ggsave("fst_rna.pdf",p3_all,width = 180,height = 150,units = "mm")

pmulti1=ggplot();pmulti2=ggplot();pmulti3=ggplot()

color=c(RColorBrewer::brewer.pal(12,"Paired"),"darkolivegreen1")

for (patient in substr(list.files(pattern = "DT.*vcf"),1,4)){
  DNA=read.table(paste0(patient,"_fst.tsv"),header = T)
  RNA=read.table(paste0(patient,"_RNA_dis.tsv"),header = T)
  physical=read.table(paste0(patient,"_physical_dis.tsv"),header = T)
  data=cbind.data.frame(DNA,RNA=RNA$value,
                        physical=physical$pos_dist)
  data$fst=data$fst
  patient_paper=all[paste0(all$Sample,"_T") == data$sample1[1],]$Patient_paper
  data$patient_paper=patient_paper
  pmulti1=pmulti1+geom_smooth(data = data, 
                              aes(x=physical,y=fst,
                                  color=patient_paper),
                              method='lm',se=FALSE,size=0.6)
  pmulti2=pmulti2+geom_smooth(data = data, 
                              aes(x=physical,y=RNA,
                                  color=patient_paper),
                              method='lm',se=FALSE,size=0.6)
  pmulti3=pmulti3+geom_smooth(data = data, 
                              aes(x=fst,y=RNA,
                                  color=patient_paper),
                              method='lm',se=FALSE,size=0.6)
}

pmulti1=pmulti1+theme_classic()+xlab("Physical distance")+ylab("FST")+
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
ggsave("multiline_physical_FST.pdf",width = 58,height = 45,units = "mm")
pmulti2=pmulti2+theme_classic()+xlab("Physical distance")+ylab("Transcriptomic distance")+
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
ggsave("multiline_physical_RNA.pdf",width = 58,height = 45,units = "mm")
pmulti3=pmulti3+theme_classic()+xlab("FST")+ylab("Transcriptomic distance")+
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
ggsave("multiline_FST_RNA.pdf",width = 58,height = 45,units = "mm")
library(ggpubr)
pp1=ggpubr::ggarrange(p1[[6]],p2[[6]],p3[[6]],
                      ncol = 3)
pp2=ggpubr::ggarrange(p1[[9]],
                      p2[[9]],p3[[9]],
                      ncol = 3)
pp3=ggpubr::ggarrange(pmulti1,
                      pmulti2,
                      pmulti3,ncol = 3,widths = c(0.8,0.8,1))
ppm=ggarrange(pp1,pp2,pp3,nrow = 3,heights = c(1,1,0.8))
ggsave("merge.pdf",width = 180,height = 180,plot = ppm,units = "mm")
ggsave("p1b2.pdf",width = 120,height = 45,units = "mm",plot = pp1)
pp3=ggpubr::ggarrange(pmulti1,
                      pmulti2,
                      pmulti3,ncol = 3,
                      widths = c(0.8,0.8,1),labels = c("c","d","e"),font.label = list(size=8))
ggsave("p1c.pdf",width = 180,height = 45,units = "mm")

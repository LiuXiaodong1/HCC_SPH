library(dplyr)
library(ape)
library(ggplot2)
library(ggpubr)
library(ggnewscale)
setwd("E:/cancer genome/liver/analysis/Fig/2/review/")
colormap=rev(RColorBrewer::brewer.pal(11,'Spectral'))

cal_ITH=function(dist_gene,two_tip){
  dist_gene=as.matrix(dist_gene)
  n_length1=dist_gene[two_tip[1],"N"]
  n_length2=dist_gene[two_tip[2],"N"]
  tlength=dist_gene[two_tip[1],two_tip[2]]
  ITH=tlength*2/(n_length2+n_length1+tlength)
  return(ITH)
}

sampleinfo = read.table("E:/cancer genome/liver/analysis/distance/all.csv",
                        header = T,sep = ",")
set.seed(123)

sample_plot=NULL;sample_plot_nolabel=NULL
sample_plot_nonblock=NULL
for (patient in unique(sampleinfo$Patient)){
  patient_paper=sampleinfo%>%filter(Patient==patient)%>%select(Patient_paper)
  patient_paper=patient_paper[1,1]
  data=sampleinfo%>%filter(Patient_paper==patient_paper)
  dis=read.table(paste0("../../1/distance/",patient,"_physical_dis.tsv"),header = T,sep = " ")
  #tree=read.nexus(paste0("../../1/tree/",patient,".tre"))
  diameter = quantile(dis$pos_dist,probs = seq(0, 1, 0.25))[3]
  select_info=sampleinfo%>%filter(Patient==patient)%>%select(X,Y,Sample_paper,Sample)
  maf_data=read.table(paste0("../../1/tree/",patient,".maf"),header = T)
  tumorNames = colnames(maf_data)[-which(colnames(maf_data) %in% c("mut_id","Variant_Classification"))]
  maf_data$N = rep(0, nrow(maf_data))
  count_mat = maf_data %>% select(tumorNames, "N") %>% mutate_at(c(tumorNames, "N"), as.numeric)
  count_mat = as.data.frame(count_mat)
  rownames(count_mat) = maf_data$mut_id
  binary_mat = count_mat
  binary_mat[binary_mat!=0] = 1
  dist_gene=dist.gene(t(binary_mat))
  
  x1=min(select_info$X)-diameter;x2=max(select_info$X)+diameter
  y1=min(select_info$Y)-diameter;y2=max(select_info$Y)+diameter
  net_x=seq(x1,x2,by=0.02)
  net_y=seq(y1,y2,by=0.02)
  all_point=as.data.frame(x=c(),col.names = c("x","y"))
  for (x in net_x){
    all_point1=as.data.frame(x=c())
    for (y in net_y){
      all_point11=cbind.data.frame(x,y)
      all_point1=rbind(all_point1,all_point11)
    }
    all_point=rbind(all_point,all_point1)
  }
  
  spearman_dist = function(x) {
    dist = 1-cor(x, method="spearman")
    return(dist)
  }
  
  dat=read.table(paste0("../../1/tree/",patient,"_DESeq2.log2_MADtop3000CodingGenes.txt"),header = T)
  tumor <- intersect(colnames(dat),select_info$Sample)
  normal<-setdiff(y=tumor,x=colnames(dat))
  dat$N<-dat[normal[1]]
  dat[normal]<-NULL
  colnames(dat)<-c(colnames(dat)[1:length(colnames(dat))-1],"N")
  rna_dist=spearman_dist(dat)
  f=c();ITH=c();rna_ITH=c()
  for (r in 1:nrow(all_point)){
    distance=apply(select_info[,1:2],1,function(x) sqrt((all_point[r,1]-x[1])**2+(all_point[r,2]-x[2])**2))
    count=length(which(distance<diameter/2))
    if (count>=2){
      f=c(f,r)
      pair=combn(select_info$Sample[which(distance<diameter/2)],2)
      ITH=c(ITH,mean(apply(pair,2,function(x) cal_ITH(dist_gene,x))))
      rna_ITH=c(rna_ITH,mean(apply(pair,2,function(x) cal_ITH(rna_dist,x))))
    }
  }
  
  
  
  all_point_ITH=cbind(all_point[f,],ITH,rna_ITH)
  cluster=readRDS("../review/res_dna_tectonic_best.RDS")
  rna_cluster=readRDS("../review/res_rna_tectonic_best.RDS")
  class=as.character(cluster[[patient]][select_info$Sample])
  rna_class=as.character(rna_cluster[[patient]][select_info$Sample])
  
  if (length(class>0)){
    class=class;rna_class=rna_class}else{class="1";rna_class="1"}
  
  non_block=F
  if (patient %in%c("DT06","DT07","DT12","DT13","DT14")){
    non_block=T;class1=dna_cluster_all[[patient]][select_info$Sample];
    rna_class1=rna_cluster_all[[patient]][select_info$Sample]
  }
  select_info=cbind(select_info,class,rna_class)
  select_info$class=as.character(select_info$class)
  select_info$rna_class=as.character(select_info$rna_class)
  color=RColorBrewer::brewer.pal(12,"Paired")[c(2,4,9)]
  color=c(color,"#E31A1C")
  names(color)=c(1,4,3,2)
  
  r=c()
  for (n in 1:length(ITH)){
    if (sample(c(0,1),1,prob=c(1-ITH[n],ITH[n]))==1){
      r=c(r,n)
    }
    #print(n)
  }
  all_point1_dna=all_point_ITH[r,]
  r=c()
  for (n in 1:length(rna_ITH)){
    if (sample(c(0,1),1,prob=c(1-rna_ITH[n],rna_ITH[n]))==1){
      r=c(r,n)
    }
    #print(n)
  }
  all_point1_rna=all_point_ITH[r,]

  p4=ggplot(data=all_point1_dna,aes(x=x,y=y))+geom_point(alpha=0)+
    stat_density_2d(geom = "raster",aes(fill = ..ndensity..),contour = F)+
    xlab("X (cm)")+ylab("Y (cm)")+
    ggtitle(paste0(patient_paper," spatial genetic heterogeneity"))+
    labs(fill="ITH")+
    scale_fill_gradientn(labels=round(seq(0,max(all_point_ITH$ITH)/1,length.out=5),digits = 2),
                         breaks=c(0,0.25,0.5,0.75,1),colours =  colormap)+
    theme_classic()+ggnewscale::new_scale_fill()+
    geom_point(data = select_info,aes(x=X,y=Y,fill=class),size=2,shape=21,stroke=0.4)+
    geom_text(data = select_info,aes(x=X,y=Y,label=substr(Sample_paper,5,8)),size=1.5,vjust=-1)+
    scale_fill_manual(values =  color[unique(class)])+labs(fill="Spatial block")+
    guides(fill = guide_legend(order=1))+
    theme(text=element_text(size = 6.5),
          legend.spacing = unit(0,"cm"),
          plot.title=element_text(hjust = 0.5,size=7),
          legend.key.width = unit(0.1, 'cm'),
          legend.key.height = unit(0.3, 'cm'),
          line = element_line(size = 0.1),
          legend.box.spacing = unit(0.01,"cm"),
          axis.ticks.length = unit(0.03,"cm"))
  p4_nolab=ggplot(data=all_point1_dna,aes(x=x,y=y))+geom_point(alpha=0)+
    stat_density_2d(geom = "raster",aes(fill = ..ndensity..),contour = F)+
    xlab("X (cm)")+ylab("Y (cm)")+
    ggtitle(paste0(patient_paper," spatial genetic heterogeneity"))+
    labs(fill="ITH")+
    scale_fill_gradientn(labels=round(seq(0,max(all_point_ITH$ITH)/1,length.out=5),digits = 2),
                         breaks=c(0,0.25,0.5,0.75,1),colours =  colormap)+
    theme_classic()+ggnewscale::new_scale_fill()+
    geom_point(data = select_info,aes(x=X,y=Y,fill=class),size=2,shape=21,stroke=0.4)+
    #geom_text(data = select_info,aes(x=X,y=Y,label=substr(Sample_paper,5,8)),size=1.5,vjust=-1)+
    scale_fill_manual(values =  color[unique(class)])+labs(fill="Spatial block")+
    guides(fill = guide_legend(order=1))+
    theme(text=element_text(size = 6.5),
          legend.spacing = unit(0,"cm"),
          plot.title=element_text(hjust = 0.5,size=7),
          legend.key.width = unit(0.1, 'cm'),
          legend.key.height = unit(0.3, 'cm'),
          line = element_line(size = 0.1),
          legend.box.spacing = unit(0.01,"cm"),
          axis.ticks.length = unit(0.03,"cm"))
  p4_rna=ggplot(data=all_point1_rna,aes(x=x,y=y))+
    stat_density_2d(geom = "raster",aes(fill = ..ndensity..),contour = F)+
    xlab("X (cm)")+ylab("Y (cm)")+
    ggtitle(paste0(patient_paper," spatial transcriptomic heterogeneity"))+labs(fill="ITH")+
    scale_fill_gradientn(labels=round(seq(0,max(all_point_ITH$rna_ITH),
                                          length.out=5),digits = 2),
                         breaks=c(0,0.25,0.5,0.75,1),colours =  colormap)+
    theme_classic()+ggnewscale::new_scale_fill()+
    geom_point(data = select_info,aes(x=X,y=Y,fill=rna_class),size=2,shape=21,stroke=0.4)+
    geom_text(data = select_info,aes(x=X,y=Y,label=substr(Sample_paper,5,8)),size=1.5,vjust=-1)+
    scale_fill_manual(values =  color[unique(rna_class)])+labs(fill="Spatial block")+
    guides(fill = guide_legend(order=1))+
    theme(text=element_text(size = 6.5),
          legend.spacing = unit(0,"cm"),
          plot.title=element_text(hjust = 0.5,size=7),
          legend.key.width = unit(0.1, 'cm'),
          legend.key.height = unit(0.3, 'cm'),
          line = element_line(size = 0.1),
          legend.box.spacing = unit(0.01,"cm"),
          axis.ticks.length = unit(0.03,"cm"))
  p4_rna_nolabel=ggplot(data=all_point1_rna,aes(x=x,y=y))+
    stat_density_2d(geom = "raster",aes(fill = ..ndensity..),contour = F)+
    xlab("X (cm)")+ylab("Y (cm)")+
    ggtitle(paste0(patient_paper," spatial transcriptomic heterogeneity"))+labs(fill="ITH")+
    scale_fill_gradientn(labels=round(seq(0,max(all_point_ITH$rna_ITH),
                                          length.out=5),digits = 2),
                         breaks=c(0,0.25,0.5,0.75,1),colours =  colormap)+
    theme_classic()+ggnewscale::new_scale_fill()+
    geom_point(data = select_info,aes(x=X,y=Y,fill=rna_class),size=2,shape=21,stroke=0.4)+
    #geom_text(data = select_info,aes(x=X,y=Y,label=substr(Sample_paper,5,8)),size=1.5,vjust=-1)+
    scale_fill_manual(values =  color[unique(rna_class)])+labs(fill="Spatial block")+
    guides(fill = guide_legend(order=1))+
    theme(text=element_text(size = 6.5),
          legend.spacing = unit(0,"cm"),
          plot.title=element_text(hjust = 0.5,size=7),
          legend.key.width = unit(0.1, 'cm'),
          legend.key.height = unit(0.3, 'cm'),
          line = element_line(size = 0.1),
          legend.box.spacing = unit(0.01,"cm"),
          axis.ticks.length = unit(0.03,"cm"))
  sample_plot[[patient]]=list(p4,p4_rna)
  sample_plot_nolabel[[patient]]=list(p4_nolab,p4_rna_nolabel)
  if (non_block){
    select_info$class1=as.character(class1)
    select_info$rna_class1=as.character(rna_class1)
    p4_nonblock=ggplot(data=all_point1_dna,aes(x=x,y=y))+geom_point(alpha=0)+
      stat_density_2d(geom = "raster",aes(fill = ..ndensity..),contour = F)+
      xlab("X (cm)")+ylab("Y (cm)")+
      ggtitle(paste0(patient_paper," spatial genetic heterogeneity"))+
      labs(fill="ITH")+
      scale_fill_gradientn(labels=round(seq(0,max(all_point_ITH$ITH)/1,length.out=5),digits = 2),
                           breaks=c(0,0.25,0.5,0.75,1),colours =  colormap)+
      theme_classic()+ggnewscale::new_scale_fill()+
      geom_point(data = select_info,aes(x=X,y=Y,fill=class1),size=2,shape=21,stroke=0.4)+
      geom_text(data = select_info,aes(x=X,y=Y,label=substr(Sample_paper,5,8)),size=1.5,vjust=-1)+
      #scale_fill_manual(values =  color[unique(class1)])+
      labs(fill="Cluster")+
      guides(fill = guide_legend(order=1))+
      theme(text=element_text(size = 6.5),
            legend.spacing = unit(0,"cm"),
            plot.title=element_text(hjust = 0.5,size=7),
            legend.key.width = unit(0.1, 'cm'),
            legend.key.height = unit(0.3, 'cm'),
            line = element_line(size = 0.1),
            legend.box.spacing = unit(0.01,"cm"),
            axis.ticks.length = unit(0.03,"cm"))
    p4_rna_nonblock=ggplot(data=all_point1_rna,aes(x=x,y=y))+
      stat_density_2d(geom = "raster",aes(fill = ..ndensity..),contour = F)+
      xlab("X (cm)")+ylab("Y (cm)")+
      ggtitle(paste0(patient_paper," spatial transcriptomic heterogeneity"))+labs(fill="ITH")+
      scale_fill_gradientn(labels=round(seq(0,max(all_point_ITH$rna_ITH),
                                            length.out=5),digits = 2),
                           breaks=c(0,0.25,0.5,0.75,1),colours =  colormap)+
      theme_classic()+ggnewscale::new_scale_fill()+
      geom_point(data = select_info,aes(x=X,y=Y,fill=rna_class1),size=2,shape=21,stroke=0.4)+
      geom_text(data = select_info,aes(x=X,y=Y,label=substr(Sample_paper,5,8)),size=1.5,vjust=-1)+
      #scale_fill_manual(values =  color[unique(rna_class1)])+
      labs(fill="Cluster")+
      guides(fill = guide_legend(order=1))+
      theme(text=element_text(size = 6.5),
            legend.spacing = unit(0,"cm"),
            plot.title=element_text(hjust = 0.5,size=7),
            legend.key.width = unit(0.1, 'cm'),
            legend.key.height = unit(0.3, 'cm'),
            line = element_line(size = 0.1),
            legend.box.spacing = unit(0.01,"cm"),
            axis.ticks.length = unit(0.03,"cm"))
    sample_plot_nonblock[[patient]]=list(p4_nonblock,p4_rna_nonblock)
  }
  ggarrange(p4,p4_rna,ncol = 2,nrow = 1)
  ggsave(paste0("ITH_",patient_paper,".pdf"),width = 5,height = 2)
}
saveRDS(sample_plot_nonblock,file = "sample_plot_nonblock.rds")
sample_plot_nonblock=readRDS("sample_plot_nonblock.rds")
sup1=ggarrange(sample_plot_nonblock[[1]][[1]],sample_plot_nonblock[[2]][[1]],
          sample_plot_nonblock[[3]][[1]],
          sample_plot_nonblock[[4]][[1]],sample_plot_nonblock[[5]][[1]],
          ncol = 5)
sup2=ggarrange(sample_plot_nonblock[[1]][[2]],sample_plot_nonblock[[2]][[2]],
          sample_plot_nonblock[[3]][[2]],
          sample_plot_nonblock[[4]][[2]],sample_plot_nonblock[[5]][[2]],
          ncol = 5)
ggarrange(sup1,sup2,nrow = 2)
ggsave("S_nonblock.pdf",width = 12,height = 4.5)
saveRDS(sample_plot_nolabel,file = "sample_plot_review_nolabel.rds")
saveRDS(sample_plot,file = "sample_plot_review.rds")
ggarrange(sample_plot[[4]][[1]],sample_plot[[5]][[1]],sample_plot[[11]][[1]],
          sample_plot[[10]][[1]],
          ncol = 4)
ggsave("fig2a.pdf",width = 10,height = 2)
ggarrange(sample_plot[[06]][[1]],sample_plot[[13]][[1]])
ggsave("fig2d.pdf",width = 5,height = 2)

sample_plot=readRDS("sample_plot_review.rds")

for (i in 1:length(sample_plot)){
  ggsave(paste0(names(sample_plot)[i],".pdf"),
         plot = sample_plot[[i]][[1]],width = 2.5,height = 2)
}
for (i in 1:length(sample_plot)){
  ggsave(paste0(names(sample_plot)[i],"rna.pdf"),
         plot = sample_plot[[i]][[2]],width = 2.5,height = 2)
}
library(tidyverse)
library(ggrepel)
library(ggpubr)
library(WeightedCluster)
library(ape)
library(ggfortify)
library(ggforce)
setwd("E:/cancer genome/liver/analysis/Fig/2/review/")
nodeHeight = function(edge, Nedge, yy) .C(node_height, 
                                          as.integer(edge[, 1]), 
                                          as.integer(edge[, 2]),
                                          as.integer(Nedge), 
                                          as.double(yy))[[4]]


nodeDepthEdgelength = function(Ntip, Nnode, edge, Nedge,edge.length) .C(node_depth_edgelength,
                                                                        as.integer(edge[,1]), 
                                                                        as.integer(edge[, 2]), 
                                                                        as.integer(Nedge), 
                                                                        as.double(edge.length),
                                                                        double(Ntip + Nnode))[[5]]
changeP <- function(p_val,r2) {
  tmp=strsplit(as.character(p_val),"e")[[1]]
  if(length(tmp)>1){
    a=tmp[1];b=tmp[2]
    as.expression(bquote(atop(p~"="~.(a)~"×"~10^~.(b)~","~R^2~"="~.(r2))))
  }else{
    as.expression(bquote(atop(p~"="~.(p_val)~","~R^2~"="~.(r2))))
  }
}
plot_physical_fst=function(patient,patient_paper){
  dir=getwd()
  setwd("E:/cancer genome/liver/analysis//Fig/1/distance/")
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
  if (patient %in% names(dna_cluster)){
    data$status=""
    for (r in 1:nrow(data)){
      if(dna_cluster[[patient]][substr(data[r,1],1,7)]==dna_cluster[[patient]][substr(data[r,2],1,7)]){
        data[r,6]="within"
      }else{data[r,6]="between"}
    }
    
    
    p=ggplot(data = data,aes(x=physical,y=fst))+
      geom_point(aes(color=status),alpha=0.6)+
      theme_classic()+
      ggtitle(patient_paper)+xlab("Physical distance")+
      ylab("FST")+theme_classic()+
      labs(subtitle =  changeP(p_val,r2))+
      theme(plot.title = element_blank(),legend.title = element_blank())+
      geom_smooth(method='lm',se=T,colour="red",linetype="dashed",size=1,fill="pink")+
      theme(text=element_text(size = 6),
            legend.spacing = unit(0,"cm"),
            plot.title=element_text(hjust = 0.5,size = 7),
            legend.key.width = unit(0.1, 'cm'),
            legend.key.height = unit(0.3, 'cm'),
            line = element_line(size = 0.1),
            legend.box.spacing = unit(0.01,"cm"),
            axis.ticks.length = unit(0.03,"cm"),
            plot.subtitle=element_text(size=5, hjust=0.5, face="italic", color="black"))
  }else{
    p=ggplot(data = data,aes(x=physical,y=fst))+
      geom_point(color="grey")+
      theme_classic()+
      ggtitle(patient_paper)+xlab("Physical distance")+
      ylab("FST")+theme_classic()+
      labs(subtitle =  changeP(p_val,r2))+
      theme(plot.title = element_blank(),legend.position = "none")+
      geom_smooth(method='lm',se=T,colour="red",linetype="dashed",size=1,fill="pink")+
      theme(text=element_text(size = 6),
            legend.spacing = unit(0,"cm"),
            plot.title=element_text(hjust = 0.5,size = 7),
            legend.key.width = unit(0.1, 'cm'),
            legend.key.height = unit(0.3, 'cm'),
            line = element_line(size = 0.1),
            legend.box.spacing = unit(0.01,"cm"),
            axis.ticks.length = unit(0.03,"cm"),
            plot.subtitle=element_text(size=5, hjust=0.5, face="italic", color="black"))
  }
  setwd(dir)
  return(p)
}
all=read.table("E:/cancer genome/liver/analysis//distance/all.csv",header = T,sep = ",")
dna_cluster=readRDS("../review/res_dna_tectonic_best.RDS")
p1_fst=NULL
for (patient in substr(list.files(path = "E:/cancer genome/liver/analysis//Fig/1/distance/",pattern = "DT.*vcf"),1,4)){
  patient_paper=all[all$Patient==patient,]$Patient_paper[1]
  p1_fst[[patient]]=plot_physical_fst(patient,patient_paper)
}


plot_physical_rna=function(patient,patient_paper){
  dir=getwd()
  setwd("E:/cancer genome/liver/analysis//Fig/1/distance/")
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
  
  if (patient %in% names(rna_cluster)){
    data$status=""
    for (r in 1:nrow(data)){
      if(dna_cluster[[patient]][substr(data[r,1],1,7)]==dna_cluster[[patient]][substr(data[r,2],1,7)]){
        data[r,6]="within"
      }else{data[r,6]="between"}}
    p=ggplot(data = data,aes(x=physical,y=RNA))+geom_point(aes(color=status),alpha=0.6)+
      theme_classic()+
      ggtitle(patient_paper)+xlab("Physical distance")+
      ylab("Transcriptomic distance")+theme_classic()+
      labs(subtitle =  changeP(p_val,r2))+
      theme(plot.title = element_text(hjust = 0.5),legend.title = element_blank())+
      geom_smooth(method='lm',se=T,colour="red",linetype="dashed",size=1,fill="pink")+
      theme(text=element_text(size = 6),
            legend.spacing = unit(0,"cm"),
            plot.title=element_text(hjust = 0.5,size=7),
            legend.key.width = unit(0.1, 'cm'),
            legend.key.height = unit(0.3, 'cm'),
            line = element_line(size = 0.1),
            legend.box.spacing = unit(0.01,"cm"),
            axis.ticks.length = unit(0.03,"cm"),
            plot.subtitle=element_text(size=5, hjust=0.5, face="italic", color="black"))
  }else{
    p=ggplot(data = data,aes(x=physical,y=RNA))+geom_point(color="grey")+
      theme_classic()+
      ggtitle(patient_paper)+xlab("Physical distance")+
      ylab("Transcriptomic distance")+theme_classic()+
      labs(subtitle =  changeP(p_val,r2))+
      theme(legend.position = "none",plot.title = element_text(hjust = 0.5))+
      geom_smooth(method='lm',se=T,colour="red",linetype="dashed",size=1,fill="pink")+
      theme(text=element_text(size = 6),
            legend.spacing = unit(0,"cm"),
            plot.title=element_text(hjust = 0.5,size=7),
            legend.key.width = unit(0.1, 'cm'),
            legend.key.height = unit(0.3, 'cm'),
            line = element_line(size = 0.1),
            legend.box.spacing = unit(0.01,"cm"),
            axis.ticks.length = unit(0.03,"cm"),
            plot.subtitle=element_text(size=5, hjust=0.5, face="italic", color="black"))
  }
  setwd(dir)
  return(p)
}

p2_rna=NULL;rna_cluster=readRDS("../review/res_rna_tectonic_best.RDS")
for (patient in substr(list.files(path = "E:/cancer genome/liver/analysis//Fig/1/distance/",pattern = "DT.*vcf"),1,4)){
  patient_paper=all[all$Patient==patient,]$Patient_paper[1]
  p2_rna[[patient]]=plot_physical_rna(patient,patient_paper)
}

all_data=read.csv("../all.csv",header = T)
all_data=all_data[which(!all_data$Patient %in% c("DT03","DT04","DT18")),]
comb=function(x,y){
  res=NULL;a=1
  for (i in x){
    for (ii in y){
      res[[as.character(a)]]=c(i,ii)
      a=a+1
    }
  }
  res=data.frame(res)
  return(res)
}
plot_tree=function(rtrw,cluster,title=""){
  Ntip <- length(rtrw$tip.label)
  Nedge <- dim(rtrw$edge)[1]
  Nnode <- rtrw$Nnode
  yy <- numeric(Ntip + Nnode)
  TIPS <- rtrw$edge[rtrw$edge[, 2] <= Ntip, 2]
  yy[TIPS] <- 1:Ntip
  z <- reorder(rtrw, order = "postorder")
  yy <- nodeHeight(edge = z$edge, Nedge = Nedge, yy = yy)
  xx <- nodeDepthEdgelength(Ntip = Ntip, Nnode = Nnode, edge = z$edge, 
                            Nedge = Nedge, edge.length = z$edge.length)
  tmp <- yy
  yy <- xx
  xx <- tmp - min(tmp) + 1
  edges = rtrw$edge
  dat = tibble(parent_index = edges[, 1],child_index = edges[, 2],
               parent_x = xx[edges[, 1]], parent_y = yy[edges[, 1]], 
               child_x = xx[edges[, 2]], child_y = yy[edges[,2]],
               length = rtrw$edge.length)
  labeling = data.frame(node_name = rtrw$tip.label,node_index = 1:length(rtrw$tip.label))
  labeling1 = data.frame(node_name = "in",
                         node_index = (length(rtrw$tip.label)+1):max(c(dat$parent_index,dat$child_index)))
  labeling1$node_name = paste0(labeling1$node_name,labeling1$node_index)
  labeling = rbind(labeling,labeling1)
  dat = merge(x = dat,y = labeling,by.x = "parent_index",by.y = "node_index",all.x = T)
  dat = dat %>% rename(parent_name = node_name)
  dat = merge(x = dat,y = labeling,by.x = "child_index",by.y = "node_index",all.x = T)
  dat = dat %>% rename(child_name = node_name)
  dat = dat %>% mutate(parent_name = as.character(parent_name),
                       child_name = as.character(child_name))
  
  dat_tips = dat[which(dat$child_name %in% rtrw$tip.label),c("child_x","child_y","child_name")]
  
  dat_tips$child_name = as.character(dat_tips$child_name)
  pon = dat$parent_name[which(dat$child_name=="N")]
  son = dat$child_name[which(dat$parent_name==pon & dat$child_name!="N")]
  son_x = dat$child_x[which(dat$child_name == son)]
  son_y = dat$child_y[which(dat$child_name == son)]
  re_len = dat$length[which(dat$child_name=="N" & dat$parent_name==pon)]
  nnn_x = son_x
  nnn_y = son_y - re_len 
  nnn_dat = tibble(child_index = NA, parent_index = NA, parent_x = son_x, parent_y = son_y,
                   child_x = nnn_x, child_y = nnn_y, length = re_len, parent_name = son, child_name = "N")
  dat1 = rbind(dat,nnn_dat)
  dat1 = dat1[-which(dat1$parent_name==pon & dat1$child_name=="N"),]
  dat1 = dat1[-which(dat1$parent_name==pon & dat1$child_name==son),]
  
  normal_x = dat1$child_x[which(dat1$child_name=="N")] 
  normal_y = dat1$child_y[which(dat1$child_name=="N")] 
  
  dat2 = dat1 
  dat2 = dat2 %>% mutate(parent_x = parent_x - normal_x, parent_y = parent_y - normal_y,
                         child_x = child_x - normal_x, child_y = child_y - normal_y)
  dat1_tips = dat2[which(dat2$child_name %in% rtrw$tip.label),c("child_x","child_y","child_name")]
  dat1_tips$child_name = as.character(dat1_tips$child_name)
  dat1_tips$child_name<-sub("DT[0-9]*_","",dat1_tips$child_name)
  
  
  
  dat1_tips$color<-as.character(c(cluster))
  
  col_tips = dat1_tips$color
  names(col_tips) = dat1_tips$child_name
  color1=RColorBrewer::brewer.pal(12,"Paired")[c(2,4,9)]
  color1=c(color1,"#E31A1C")
  names(color1)=c(1,4,3,2)
  patient_paper=all[all$Patient==patient,]$Patient_paper[1]
  rot_tree = ggplot(data = dat2) +
    geom_segment(aes(x=parent_x, y=parent_y, xend=child_x, yend=child_y),size=0.4) +
    geom_point(data = dat1_tips, 
               aes(x = child_x,y = child_y,fill=color),
               size=5,shape=21,stroke=0.4) +
    scale_fill_manual(values=color1) +
    xlab("")+ggtitle(title)+
    theme_void() + 
    theme(legend.position = "none",plot.margin=unit(c(0,0,0,0),'lines'),
          axis.title.x = element_text(size =6,hjust =0.7),
          plot.title = element_text(hjust = 0.5,size = 7))
  return(rot_tree)
}


dna_cluster=readRDS("../review/res_dna_tectonic_best.RDS")
rna_cluster=readRDS("../review/res_rna_tectonic_best.RDS")
pca_dna=NULL;dna_tree1=NULL
pca_rna=NULL;rna_tree1=NULL;
pca_dna_notectonic=NULL
pca_rna_notectonic=NULL
for (patient in c(unique(all_data$Patient))[c(3,4,8,9,10)]){
  patient_paper=all[all$Patient==patient,]$Patient_paper[1]
  maf= paste0("../../1/tree/",patient,".maf")
  maf_data=read.table(maf,header = T)
  tumorNames = colnames(maf_data)[-which(colnames(maf_data) %in% c("mut_id","Variant_Classification"))]
  count_mat = maf_data %>% select(tumorNames) %>% mutate_at(c(tumorNames), as.numeric)
  count_mat = as.data.frame(count_mat)
  rownames(count_mat) = maf_data$mut_id
  binary_mat = count_mat
  binary_mat[binary_mat!=0] = 1
  tree=ape::read.tree(paste0("../..//1/reveiew/",patient,"par.tre"))
  cluster=dna_cluster[[patient]]
  cluster1=c(cluster,"N"=4)
  cluster1=cluster1[tree$tip.label]
  
  dna_tree1[[patient]]=plot_tree(tree,cluster1,
                                 title = paste0(patient_paper," genetic tree"))
  dna_tree=~plot.phylo(tree,tip.color=cluster1,main=paste0(patient,"DNA"))
  dna_tree=ggplotify::as.grob(dna_tree)
  # color=RColorBrewer::brewer.pal(6,"Set1")
  # names(color)=c(1,2,4,3,5,6)
  color=RColorBrewer::brewer.pal(12,"Paired")[c(2,4,9)]
  color=c(color,"#E31A1C")
  names(color)=c(1,4,3,2)
  
  pca=prcomp(t(binary_mat))
  cluster=cluster[colnames(binary_mat)]
  data=as.data.frame(cluster)
  data$cluster=as.character(data$cluster)
  
  pca_dna[[patient]]=autoplot(pca,data = data,colour="cluster")+
    scale_color_manual(values =  color[as.character(unique(cluster))])+
    ggtitle(paste0(patient_paper," genetic PCA"))+theme_classic()+
    theme(text=element_text(size = 6),
          legend.spacing = unit(0,"cm"),
          legend.title = element_blank(),
          plot.title=element_text(hjust = 0.5,size=7),
          legend.key.width = unit(0.1, 'cm'),
          legend.key.height = unit(0.3, 'cm'),
          line = element_line(size = 0.1),
          legend.box.spacing = unit(0.01,"cm"),
          axis.ticks.length = unit(0.03,"cm"))
  
  
  RNA=read.table(paste0("../../1/distance/",patient,"_DESeq2.log2_top3000CodingGenes.txt"),header = T)  
  
  cluster_rna=rna_cluster[[patient]]
  number_cluster=unique(cluster_rna)
  rna=read.table(paste0("../../1/distance/",patient,"_rna_dis.tsv"),header = T)
  tree=ape::read.tree(paste0("../../1/tree/",patient,"RNA.tre"))
  cluster1_rna=c(cluster_rna,"N"=4)
  cluster1_rna=cluster1_rna[tree$tip.label]
  
  
  rna_tree=~plot.phylo(tree,tip.color=cluster1_rna,main=paste0(patient,"RNA"))
  rna_tree=ggplotify::as.grob(rna_tree)
  rna_tree1[[patient]]=plot_tree(tree,cluster = cluster1_rna,
                                 title = paste0(patient_paper," transcriptomic tree"))

  pca=prcomp(t(RNA))
  cluster_rna=cluster_rna[colnames(RNA)]
  data=as.data.frame(cluster_rna)
  data$cluster=as.character(data$cluster_rna)
  pca_rna[[patient]]=autoplot(pca,data = data,colour="cluster")+
    scale_color_manual(values =  color[as.character(unique(cluster))])+
    ggtitle(paste0(patient_paper," transcriptomic PCA"))+theme_classic()+
    theme(text=element_text(size = 6),
          legend.title = element_blank(),
          legend.spacing = unit(0,"cm"),
          plot.title=element_text(hjust = 0.5,size=7),
          legend.key.width = unit(0.1, 'cm'),
          legend.key.height = unit(0.3, 'cm'),
          line = element_line(size = 0.1),
          legend.box.spacing = unit(0.01,"cm"),
          axis.ticks.length = unit(0.03,"cm"))
}
SH05_fst=p1_fst[["DT09"]]+
  scale_color_manual(values =c("grey","grey"))+
  theme(legend.position = "none")
p1=ggarrange(SH05_fst,dna_tree1[["DT09"]],pca_dna[["DT09"]],ncol = 3,
             labels = c("a","b","c"),font.label = list(size=8))
p_chindex=readRDS("p_chindex.RDS")
sample_plot=readRDS("sample_plot_review_nolabel.rds")

p2=ggarrange(p_chindex,sample_plot$DT09[[1]],sample_plot$DT07[[1]],ncol = 3,
             widths = c(2,1,1),
             labels = c("d","e","f"),font.label = list(size=8))

p3=ggarrange(sample_plot$DT10[[1]],sample_plot$DT16[[1]],sample_plot$DT17[[1]],
             sample_plot$DT19[[1]],ncol=4,
             labels = c("g","",""),font.label = list(size=8))

p4=ggarrange(p2_rna[["DT09"]],rna_tree1[["DT09"]],pca_rna[["DT09"]],ncol=3,
             labels = c("h","i","j"),font.label = list(size=8))
pdf("Figure2_review2.pdf",width = 10,height = 10)
ggarrange(p1,p2,p3,p4,nrow = 4,heights = c(3,2,2,3))
dev.off()

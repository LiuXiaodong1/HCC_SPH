library(tidyverse)
library(ape)
library(ggfortify)


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
  dat = dat %>% rename("node_name"="parent_name"   )
  dat = merge(x = dat,y = labeling,by.x = "child_index",by.y = "node_index",all.x = T)
  dat = dat %>% rename("node_name"="child_name"  )
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
  dat1_tips$child_name<-sub("SH[0-9]*_","",dat1_tips$child_name)
  
  dat1_tips$color<-as.character(c(cluster))
  color1=RColorBrewer::brewer.pal(12,"Paired")[c(2,4,9)]
  color1=c(color1,"#E31A1C")
  names(color1)=c(1,4,3,2)
  rot_tree = ggplot(data = dat2) +
    geom_segment(aes(x=parent_x, y=parent_y,
                     xend=child_x, yend=child_y),size=0.4) +
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

Figure2_data=read_tsv("../Data/Figure2.tsv")
RNA_tree=read.tree("../Data/SH05RNA.tre")
cluster=Figure2_data%>%filter(Patient=="SH05")%>%pull(`RNA block`)
names(cluster)=Figure2_data%>%filter(Patient=="SH05")%>%pull(Sample)
cluster1=c(cluster,"N"=4)
cluster1=cluster1[RNA_tree$tip.label]
Fig2i=plot_tree(RNA_tree,cluster = cluster1,
                title = "SH05 transcriptomic tree")


DNA_tree=read.tree("../Data/SH05DNA.tre")

cluster=Figure2_data%>%filter(Patient=="SH05")%>%pull(`DNA block`)
names(cluster)=Figure2_data%>%filter(Patient=="SH05")%>%pull(Sample)
cluster1=c(cluster,"N"=4)
cluster1=cluster1[DNA_tree$tip.label]
Fig2b=plot_tree(DNA_tree,cluster = cluster1,
                title = "SH05 genetic tree")

color=RColorBrewer::brewer.pal(12,"Paired")[c(2,4,9)]
color=c(color,"#E31A1C")
names(color)=c(1,4,3,2)
binary_mat=read.table("../Data/Figure2c.tsv")
pca=prcomp(t(binary_mat))
cluster=cluster[colnames(binary_mat)]
data=as.data.frame(cluster)
data$cluster=as.character(data$cluster)
Fig2c=autoplot(pca,data = data,colour="cluster")+
  scale_color_manual(values =  color[as.character(unique(cluster))])+
  ggtitle("SH05 genetic PCA")+theme_classic()+
  theme(text=element_text(size = 6),
        legend.spacing = unit(0,"cm"),
        legend.title = element_blank(),
        plot.title=element_text(hjust = 0.5,size=7),
        legend.key.width = unit(0.1, 'cm'),
        legend.key.height = unit(0.3, 'cm'),
        line = element_line(size = 0.1),
        legend.box.spacing = unit(0.01,"cm"),
        axis.ticks.length = unit(0.03,"cm"))

CHindex_info=read_tsv("../Data/Figure2d.tsv")
#Order by CHindex

CHindex_info$Patient=factor(
  CHindex_info$Patient,
  levels=unique(
    CHindex_info$Patient[order(CHindex_info$CHindex,decreasing = T)]
    )
  )
CHindex_info=CHindex_info%>%mutate(Cluster=as.character(Cluster))
Fig2d=ggplot(data=CHindex_info)+
  geom_bar(aes(x=Patient,y=CHindex,fill=Cluster),
           stat = "identity",
           position = position_dodge())+
  scale_y_continuous(expand = c(0,NA))+
  xlab("Patient")+ylab("CHindex")+theme_classic()+
  theme(text = element_text(size=5))

Fig2b

Fig2c

Fig2d

Fig2i

library(tidyverse)

class1=c("DT10_01_T","DT10_02_T","DT10_03_T","DT10_06_T","DT10_07_T",
         "DT10_08_T","DT10_11_T","DT10_12_T")

class2=c("DT10_04_T","DT10_05_T","DT10_09_T","DT10_10_T","DT10_13_T",
         "DT10_14_T","DT10_15_T")

all_seg=read.table("E:/cancer genome/liver/分析/gistic_out/3.0//all_segments.txt",header = T,sep="\t")
cyto=read_tsv("E:/cancer genome/liver/分析/gistic_out/3.0/cytoBand.txt",col_names = F)
cyto=cyto%>%filter(X1%in%c(paste0("chr",c(1:22,"X"))))
genome_len=sum(cyto$X3-cyto$X2)

GII1=c();breakpoint1=c();cnv1=c()
for (s in class1){
  tmp=all_seg%>%filter(sample==s,depth.ratio>0.2|depth.ratio< -0.2)
  len=sum(tmp$end.pos-tmp$start.pos)
  breakpoint1=c(breakpoint1,nrow(tmp))
  cnv1=c(cnv1,mean(abs(tmp$depth.ratio)))
  GII1=c(GII1,len/genome_len)
}
mean(GII1);mean(breakpoint1)

GII2=c();breakpoint2=c();cnv2=c()
for (s in class2){
  tmp=all_seg%>%filter(sample==s,depth.ratio>0.2|depth.ratio< -0.2)
  len=sum(tmp$end.pos-tmp$start.pos)
  breakpoint2=c(breakpoint2,nrow(tmp))
  cnv2=c(cnv2,mean(abs(tmp$depth.ratio)))
  GII2=c(GII2,len/genome_len)
}

library(ggpubr)
mycompare=list(c(1,2))
GII1=cbind.data.frame(GII=GII1,class="GII1")
GII2=cbind.data.frame(GII=GII2,class="GII2")
GII=rbind(GII1,GII2)
GII$class=factor(GII$class,level=c("GII2","GII1"))
p_10_GII=ggplot(data = GII,aes(class,GII))+
  geom_boxplot(aes(fill=class),outlier.alpha = 0,size=0.1)+
  geom_point(size=0.5)+scale_x_discrete(labels=c("Classic","Poliferation"))+
  xlab("Subtype")+ggtitle("SH06")+ggplot2::ylim(0,0.62)+
  scale_fill_manual(values=c(RColorBrewer::brewer.pal(7,"Set2")[c(1,4)]))+
  theme_classic()+stat_compare_means(comparisons = mycompare,method = "t.test",size=2)+
  theme(text=element_text(size = 6.5),
        legend.position = "none",
        plot.title=element_text(hjust = 0.5,size=7),
        line = element_line(size = 0.1),
        axis.ticks.length = unit(0.03,"cm"))

breakpoint1=cbind.data.frame(breakpoint=breakpoint1,class="breakpoint1")
breakpoint2=cbind.data.frame(breakpoint=breakpoint2,class="breakpoint2")
breakpoint=rbind.data.frame(breakpoint1,breakpoint2)


p_10_breakpoint=ggplot(data = breakpoint,aes(class,breakpoint))+
  geom_boxplot(aes(fill=class))+
  geom_point()+scale_x_discrete(labels=c(1,2,3,4))+
  scale_fill_manual(values=c(RColorBrewer::brewer.pal(7,"Set2")[c(1,2,4,5)]))+
  theme_classic()+stat_compare_means(comparisons = mycompare,method = "t.test",size=2)+
  theme(legend.position = "none",
        text = element_text(size = 5),
        plot.title = element_text(hjust = 0.5))
ggarrange(p_10_GII,p_10_breakpoint,ncol = 2)


exp_counts = read.table("E:/cancer genome/liver/分析/DEseq2-normalized_clustering/DT03-DT17/outputs/S01_expCounts_MeanDupCodingGenes.txt",header=TRUE)
exp_counts<-exp_counts[,which(substr(colnames(exp_counts),1,4)=="DT09")]
exp_counts<-exp_counts[,-c(29,30)]
coldata = data.frame(colnames(exp_counts))

library(DESeq2)
library(clusterProfiler)
library(org.Hs.eg.db)
coldata$condition=c(rep("CLASS01", 9), rep("CLASS02", 2),
                    rep("CLASS01", 5), rep("CLASS02", 2),
                    rep("CLASS01", 11),  rep("CLASS02", 1)
)
dds = DESeqDataSetFromMatrix(countData=round(exp_counts), 
                             colData=coldata, design=~condition) # ~1 means no design
dds <- DESeq(dds)
res = results(dds, contrast=c("condition", "CLASS01","CLASS02"))
res = res[order(res$pvalue),]

diff_gene_deseq2_up <-subset(res,padj < 0.05 & (log2FoldChange > 1 ))
dim(diff_gene_deseq2_up)

diff_gene_deseq2_down <-subset(res,padj < 0.05 & (log2FoldChange < -1))
dim(diff_gene_deseq2_down)
x_up=rownames(diff_gene_deseq2_up)
up = bitr(x_up, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
x_down=rownames(diff_gene_deseq2_down)
down = bitr(x_down, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
DT09_up <- enrichGO(up$ENTREZID, "org.Hs.eg.db", keyType = "ENTREZID",ont = 'ALL',
                   pvalueCutoff  = 0.05,pAdjustMethod = "BH",  
                   qvalueCutoff  = 0.1, readable=T)

result1=DT09_up@result[1:5,]
result1$qvalue=-log(result1$qvalue,10)

DT09_down <- enrichGO(down$ENTREZID, "org.Hs.eg.db", keyType = "ENTREZID",ont = 'ALL',
                     pvalueCutoff  = 0.05,pAdjustMethod = "BH",  qvalueCutoff  = 0.1, 
                     readable=T)
result2=DT09_down@result[1:5,]
result2$qvalue=log(result2$qvalue,10)

result1$fill="red"
result2$fill="blue"
result_DT09=rbind(result1,result2)
result_DT09$Description=factor(c(result1$Description,
                            result2$Description),
                          levels = c(rev(result2$Description),
                                     rev(result1$Description)
                          ))
library(ggbreak)
GO_DT09=ggplot(result_DT09)+geom_bar(aes(x=Description,y=qvalue,fill=fill),
                                stat = "identity",width = 0.3)+
  scale_fill_manual(values =  RColorBrewer::brewer.pal(n = 3,"Set2"))+
  scale_y_continuous(expand = c(0,0),limits = c(-20,42),
                     breaks=c(-15,-5,-1,0,1,10,20),labels = c(15,5,1,0,1,10,20))+
  scale_y_break(c(22,38),scales = 0.2,ticklabels=c(40,41))+
  coord_flip()+
  scale_x_discrete()+
  xlab("")+ggtitle("SH05")+
  ylab("Enrichment significance (-log10 (q))")+theme_classic()+
  geom_hline(aes(yintercept=1),color="red",linetype=2)+
  geom_hline(aes(yintercept=0))+
  geom_hline(aes(yintercept=-1),color="red",linetype=2)+
  theme(panel.background = element_blank(),
        legend.position = "none",
        plot.title = element_text(size=6,hjust = 0.8),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x.bottom  = element_blank(),
        axis.text = element_text(size = 6),
        axis.title.x = element_text(size=6,hjust = 0.8),
        axis.line = element_line(size = 0.3),
        axis.ticks = element_line(size=0.3))



exp_counts = read.table("E:/cancer genome/liver/分析/DEseq2-normalized_clustering/DT03-DT17/outputs/S01_expCounts_MeanDupCodingGenes.txt",header=TRUE)
exp_counts<-exp_counts[,which(substr(colnames(exp_counts),1,4)=="DT10")]
exp_counts<-exp_counts[,-c(16,17)]
coldata = data.frame(colnames(exp_counts))


coldata$condition=c(rep("CLASS01", 3), rep("CLASS02", 2),
                    rep("CLASS01", 3), rep("CLASS02", 2),
                    rep("CLASS01", 2), rep("CLASS02", 3)
)
dds = DESeqDataSetFromMatrix(countData=round(exp_counts), 
                             colData=coldata, design=~condition) # ~1 means no design
dds <- DESeq(dds)
#for (i in c("03")){
res = results(dds, contrast=c("condition","CLASS01", "CLASS02"))
res = res[order(res$pvalue),]
head(res)
summary(res)
#write.csv(res,file="All_results.csv")
diff_gene_deseq2_up <-subset(res,padj < 0.05 & (log2FoldChange > 1 ))
dim(diff_gene_deseq2_up)
#write.csv(diff_gene_deseq2,file= "tumor_normal_up.csv")

diff_gene_deseq2_down <-subset(res,padj < 0.05 & (log2FoldChange < -1))
dim(diff_gene_deseq2_down)
# write.csv(diff_gene_deseq2,file= "tumor_normal_down.csv")
#head(diff_gene_deseq2)
x_up=rownames(diff_gene_deseq2_up)
up = bitr(x_up, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
x_down=rownames(diff_gene_deseq2_down)
down = bitr(x_down, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
DT10_up <- enrichGO(up$ENTREZID, "org.Hs.eg.db", keyType = "ENTREZID",ont = 'ALL',
                   pvalueCutoff  = 0.05,pAdjustMethod = "BH",  qvalueCutoff  = 0.1, readable=T)
result1=DT10_up@result[1:5,]
result1$qvalue=-log(result1$qvalue,10)

DT10_down <- enrichGO(down$ENTREZID, "org.Hs.eg.db", keyType = "ENTREZID",ont = 'ALL',
                     pvalueCutoff  = 0.05,pAdjustMethod = "BH",  qvalueCutoff  = 0.1, readable=T)

result2=DT10_down@result[1:5,]
result2$qvalue=log(result2$qvalue,10)

result1$fill="red"
result2$fill="blue"
result_DT10=rbind(result1,result2)
result_DT10$Description=factor(c(result1$Description,
                            result2$Description),
                          levels = c(rev(result2$Description),
                                     rev(result1$Description)
                          ))
library(ggbreak)
GO_DT10=ggplot(result_DT10)+geom_bar(aes(x=Description,y=qvalue,fill=fill),
                                stat = "identity",width = 0.3)+
  scale_fill_manual(values =  RColorBrewer::brewer.pal(n = 3,"Set2"))+
  scale_y_continuous(expand = c(0,0),
                     breaks=c(-5,0,10,20),labels = c(5,0,10,20))+
  scale_y_break(c(5,18),scales = 0.5,ticklabels=c(20,30,40))+
  coord_flip()+
  scale_x_discrete()+
  xlab("")+ggtitle("GO enrichment for SH06")+
  ylab("Enrichment significance (-log10 (q))")+theme_classic()+
  #geom_hline(aes(yintercept=1),color="red",linetype=2)+
  geom_hline(aes(yintercept=0),size=0.1)+
  #geom_hline(aes(yintercept=-1),color="red",linetype=2)+
  theme(panel.background = element_blank(),
        legend.position = "none",
        plot.title =  element_text(size=6,hjust = 0.5),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
       # axis.title.x.bottom  = element_blank(),
        axis.title.x = element_text(size=6,hjust = 0.8),
        axis.text = element_text(size = 6),
        axis.line = element_line(size = 0.3),
        axis.ticks = element_line(size=0.3))

exp_counts = read.table("E:/cancer genome/liver/分析/DEseq2-normalized_clustering/DT03-DT19/outputs/S01_expCounts_MeanDupCodingGenes.txt",header=TRUE)
exp_counts<-exp_counts[,which(substr(colnames(exp_counts),1,4)=="DT16")]
exp_counts<-exp_counts[,-c(28)]
coldata = data.frame(colnames(exp_counts))


coldata$condition=c(rep("CLASS01", 10), rep("CLASS02", 1),
                    rep("CLASS01", 5), rep("CLASS02", 2),
                    rep("CLASS01", 3), rep("CLASS02", 3),
                    rep("CLASS01", 2), rep("CLASS02", 1)
)
dds = DESeqDataSetFromMatrix(countData=round(exp_counts), 
                             colData=coldata, design=~condition) # ~1 means no design
dds <- DESeq(dds)
#for (i in c("03")){
res = results(dds, contrast=c("condition","CLASS01", "CLASS02"))
res = res[order(res$pvalue),]
head(res)
summary(res)
#write.csv(res,file="All_results.csv")
diff_gene_deseq2_up <-subset(res,padj < 0.05 & (log2FoldChange > 1 ))
dim(diff_gene_deseq2_up)
#write.csv(diff_gene_deseq2,file= "tumor_normal_up.csv")

diff_gene_deseq2_down <-subset(res,padj < 0.05 & (log2FoldChange < -1))
dim(diff_gene_deseq2_down)
# write.csv(diff_gene_deseq2,file= "tumor_normal_down.csv")
#head(diff_gene_deseq2)
x_up=rownames(diff_gene_deseq2_up)
up = bitr(x_up, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
x_down=rownames(diff_gene_deseq2_down)
down = bitr(x_down, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
DT16_up <- enrichGO(up$ENTREZID, "org.Hs.eg.db", keyType = "ENTREZID",ont = 'ALL',
                   pvalueCutoff  = 0.05,pAdjustMethod = "BH",  qvalueCutoff  = 0.1, readable=T)
result1=DT16_up@result[1:5,]
result1$qvalue=-log(result1$qvalue,10)

DT16_down <- enrichGO(down$ENTREZID, "org.Hs.eg.db", keyType = "ENTREZID",ont = 'ALL',
                     pvalueCutoff  = 0.05,pAdjustMethod = "BH",  qvalueCutoff  = 0.1, readable=T)

result2=DT16_down@result[1:5,]
result2$qvalue=log(result2$qvalue,10)

result1$fill="red"
result2$fill="blue"
result_DT16=rbind(result1,result2)
result_DT16$Description=factor(c(result1$Description,
                            result2$Description),
                          levels = c(rev(result2$Description),
                                     rev(result1$Description)
                          ))
library(ggbreak)
GO_DT16=ggplot(result_DT16)+geom_bar(aes(x=Description,y=qvalue,fill=fill),
                                stat = "identity",width = 0.3)+
  scale_fill_manual(values =  RColorBrewer::brewer.pal(n = 3,"Set2"))+
  scale_y_continuous(expand = c(0,0),
                     breaks=c(-60,-40,-20,-1,0,1,10,20),
                     labels = c(60,40,20,1,0,1,10,20))+
  scale_y_break(c(5,18),scales = 0.5,ticklabels=c(20,25,30))+
  scale_y_break(c(-5,-25),scales = 0.5)+
  coord_flip()+
  scale_x_discrete()+
  xlab("")+ggtitle("SH10")+
  theme_classic()+
  geom_hline(aes(yintercept=1),color="red",linetype=2)+
  geom_hline(aes(yintercept=0))+
  geom_hline(aes(yintercept=-1),color="red",linetype=2)+
  ylab("Enrichment significance (-log10 (q))")+
  theme(panel.background = element_blank(),
        legend.position = "none",
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(size=6,hjust = 0.8),
        axis.title.x.bottom  = element_blank(),
        axis.title.x = element_text(size=6,hjust = 0.8),
        axis.text = element_text(size = 6),
        axis.line = element_line(size = 0.3),
        axis.ticks = element_line(size=0.3))


library(ape)
library(ggtree)
wb<-read.csv("E:/cancer genome/liver/分析/distance/all.csv")
tree18=treeio::read.nexus(file = "../1/tree/DT18.tre")
tree18$tip.label[-length(tree18$tip.label)]=wb%>%filter(Sample %in% tree18$tip.label)%>%
  dplyr::pull(Sample_paper)%>%substr(.,6,7)


color=c(RColorBrewer::brewer.pal(7,"Set2")[c(5,1,2,4)])
names(color)=c("1","2","3","4")


tree_DT18=ggtree(tree18)+ geom_tiplab(size=1.5,hjust = -0.5)+
  geom_rootedge(rootedge = 10)+
  geom_label2(aes(subset=branch.length==19.625,label="CTNNB1"),size=1)+
  geom_tippoint(aes(color=c(color[2]),subset=(label%in%c("01","02"))),
                size=1.5,show.legend = F)+
  geom_tippoint(aes(color=c(color[1]),subset=(label%in%c("03","04","05"))),
                size=1.5,show.legend = F)+ggtitle("SH12")+
  theme(plot.title = element_text(size = 6,hjust = 0.5))
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


plot_tree=function(rtrw,cluster,title="",label=F){
  
  #cluster=cluster[rtrw$tip.label]
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
  #dat1_tips$child_name<-sub("DT[0-9]*_","",dat1_tips$child_name)
  
  
  
  dat1_tips$color<-as.character(c(cluster))
  
  col_tips = dat1_tips$color
  names(col_tips) = dat1_tips$child_name
  # color1=RColorBrewer::brewer.pal(6,"Set1")
  # names(color1)=c(1,2,4,3,5,6)
  color1=RColorBrewer::brewer.pal(12,"Paired")[c(2,4,9)]
  color1=c(color1,"#E31A1C")
  names(color1)=c(1,4,3,2)
  all=read.csv("E:/cancer genome/liver/分析/distance/all.csv",header = T)
  patient_paper=all[all$Patient==patient,]$Patient_paper[1]
  sample_name=all$Sample_paper
  names(sample_name)=all$Sample
  
  cluster=as.character(cluster)
  rot_tree = ggplot(data = dat2) + 
      geom_segment(aes(x=parent_x, y=parent_y, xend=child_x, yend=child_y),
                   size=0.4,alpha=1) +
      geom_point(data = dat1_tips, 
                 aes(x = child_x,y = child_y,fill=cluster),
                 size=5,shape=21,stroke=0.4) +
      scale_fill_manual(values=c(RColorBrewer::brewer.pal(7,"Set2")[c(5,1)],"white")) +
      xlab("")+ggtitle(title)+
      geom_text(data = dat1_tips,color="black",
                aes(label=child_name,x = child_x,y = child_y),
                hjust=0.5, vjust=0.5, size=2.5) +
    geom_text(aes(x=-0.7,y=160,label="CTNNB1"),size=3)+
      theme_void() +
      theme(legend.position = "none",plot.margin=unit(c(0,0,0,0),'lines'),
            axis.title.x = element_text(size =6,hjust =0.7),
            plot.title = element_text(hjust = 0.5,size = 8))
  return(rot_tree)
}
tree18$edge.length[10]=tree18$edge.length[10]-200
patient="DT16"
tree_DT18_2= plot_tree(tree18,cluster=c(1,1,2,2,2,"3"),title = "SH12")

library(ggforce)

allplot=NULL
class=read.table("../3/class.tsv")
rownames(class)=class$V1
color=c(RColorBrewer::brewer.pal(7,"Set2")[c(5,1,2,4)])
names(color)=c("1","2","3","4")
for (p in unique(wb$Patient) ){
  patient<-wb[which(wb$Patient==p),]
  patient$class=as.character(class[patient$Sample,]$V2)
  mx=max(patient$X)+min(patient$X)
  my=max(patient$Y)+min(patient$Y)
  r=sqrt(mx**2/4+my**2/4)/1.3
  print(r)
  d=as.data.frame(cbind(X=mx/2,Y=my/2,r=r))
  x_alt=2.9-mx/2
  patient$X=patient$X+x_alt
  d$X=d$X+x_alt
  p_paper=patient$Patient_paper[1]
  if (! p %in% c("DT03","DT04")){
    allplot[[p]]=ggplot()+
      geom_circle(data=patient,aes(x0=X, y0=Y,r=0.15,fill=class))+coord_fixed()+
      scale_fill_manual(values=color)+
      #geom_text(data=patient,aes(label=substr(Sample_paper,6,7),x=X,y=Y),size=3)+
      xlab("X(cm)") +#scale_x_continuous(limits=c(0, 5.8))+
      labs(fill="")+
      ylab("Y(cm)") +#scale_y_continuous(limits=c(-1, 5.5))+
      ggtitle(paste0("",p_paper))+
      theme_void()+geom_circle(data=d,aes(x0=X, y0=Y,r=r))+
      theme(plot.title = element_text(hjust = 0.5,size = 6),legend.position = "none")
  }else{
    allplot[[p]]=ggplot()+
      geom_circle(data=patient,aes(x0=X, y0=Y,r=0.25,fill=class))+#coord_fixed()+
      scale_fill_manual(values=color)+
      #geom_text(data=patient,aes(label=substr(Sample_paper,6,7),x=X,y=Y),size=3)+
      xlab("X(cm)") +#scale_x_continuous(limits=c(0, 5.8))+
      labs(fill="")+
      ylab("Y(cm)") +#scale_y_continuous(limits=c(-1, 5.5))+
      ggtitle(paste0("",p_paper))+
      theme_void()+geom_circle(data=d,aes(x0=X, y0=Y,r=r))+
      theme(plot.title = element_text(hjust = 0.5,size = 6),legend.position = "none")
  }
}
pdf("p4_elements1.pdf",width=2,height = 2)
allplot[[5]]
allplot[[6]]
allplot[[11]]
allplot[[12]]
p_10_GII
dev.off()
pdf("p4_elements2.pdf",width=3,height = 2)
GO_DT09
GO_DT10
GO_DT16
tree_DT18
dev.off()



d=rbind(DT16_up@result[1:5,],DT09_up@result[1:5,],DT16_down@result[1:5,],DT09_down@result[1:5,])
d_up=rbind(
DT16_up@result[which(DT16_up@result$Description%in%d$Description),],
DT09_up@result[which(DT09_up@result$Description%in%d$Description),]
)
tmp=d_up%>%group_by(Description)%>%summarise(min(qvalue))%>%arrange(.,`min(qvalue)`)

d_down=rbind(
  DT16_down@result[which(DT16_down@result$Description%in%d$Description),],
  DT09_down@result[which(DT09_down@result$Description%in%d$Description),]
)
tmp1=d_down%>%group_by(Description)%>%summarise(min(qvalue))%>%arrange(.,`min(qvalue)`)

df=rbind(cbind(result_DT09,sample="SH05"),cbind(result_DT16,sample="SH10"))


df$Description=factor(df$Description,levels = rev(c(tmp$Description,tmp1$Description)))

color_1=RColorBrewer::brewer.pal(n = 3,"Set2")
# df$fill1=paste0(df$fill,df$sample,recycle0 = T)

one=names(table(df$Description)[table(df$Description)==1])
df1=df[which(df$Description%in%one),]
df1$qvalue=0
df1$sample=sub("SH05","tmp",df1$sample)
df1$sample=sub("SH10","SH05",df1$sample)
df1$sample=sub("tmp","SH10",df1$sample)
# df1$fill1=sub(05,06,df1$fill1)
df=rbind(df,df1)
#df$fill1=sub(".*SH","SH",df$fill1)

p_go=ggplot(df)+geom_bar(aes(x=Description,y=qvalue,fill=sample),
                             stat = "identity",width = 0.3,position = position_dodge())+
  scale_fill_manual(values =c(RColorBrewer::brewer.pal(12,"Set3")[c(5,6)]))+
  scale_y_continuous(expand = c(0,0),
                     breaks=c(-50,-25,0,25),
                     labels = c(50,25,0,25))+
  #scale_y_break(c(5,18),scales = 0.5,ticklabels=c(20,25,30))+
  #scale_y_break(c(-5,-25),scales = 0.5)+
  coord_flip()+
  scale_x_discrete()+
  xlab("")+ggtitle("GO enrichment for SH05/SH10")+
  theme_classic()+
  #geom_hline(aes(yintercept=1),color="red",linetype=2)+
  geom_hline(aes(yintercept=0),size=0.1)+
  #geom_hline(aes(yintercept=-1),color="red",linetype=2)+
  ylab("Enrichment significance (-log10 (q))")+
  theme(panel.background = element_blank(),
        #legend.position = "none",
        legend.spacing = unit(0,"cm"),
        legend.key.width = unit(0.1, 'cm'),
        legend.key.height = unit(0.3, 'cm'),
        legend.box.spacing = unit(0.01,"cm"),
        legend.title = element_blank(),
        legend.text = element_text(size=5),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(size=6,hjust = 0.8),
        #axis.title.x.bottom  = element_blank(),
        axis.title.x = element_text(size=6,hjust = 0.8),
        axis.text = element_text(size = 6),
        axis.line = element_line(size = 0.3),
        axis.ticks = element_line(size=0.3))
#ggarrange(allplot[[05]],p_go,allplot[[10]],ncol = 3)
pdf("p_go.pdf",width=4,height = 2)
p_go
dev.off()


library(dplyr)
library(ape)
library(ggplot2)
library(ggpubr)
library(ggnewscale)

colormap=rev(RColorBrewer::brewer.pal(11,'Spectral'))

cal_ITH=function(dist_gene,two_tip){
  dist_gene=as.matrix(dist_gene)
  n_length1=dist_gene[two_tip[1],"N"]
  n_length2=dist_gene[two_tip[2],"N"]
  tlength=dist_gene[two_tip[1],two_tip[2]]
  ITH=tlength*2/(n_length2+n_length1+tlength)
  return(ITH)
}

sampleinfo = read.table("E:/cancer genome/liver/分析/distance/all.csv",
                        header = T,sep = ",")
set.seed(123)

sample_plot=NULL
# for (patient in unique(sampleinfo$Patient)){
for (patient in c("DT09","DT10","DT16")){
  patient_paper=sampleinfo%>%filter(Patient==patient)%>%dplyr::select(Patient_paper)
  patient_paper=patient_paper[1,1]
  data=sampleinfo%>%filter(Patient_paper==patient_paper)
  dis=read.table(paste0("../1/distance/",patient,"_physical_dis.tsv"),header = T,sep = " ")
  tree=read.nexus(paste0("../1/tree/",patient,".tre"))
  diameter = quantile(dis$pos_dist,probs = seq(0, 1, 0.25))[3]
  select_info=sampleinfo%>%filter(Patient==patient)%>%dplyr::select(X,Y,Sample_paper,Sample)
  maf_data=read.table(paste0("../1/tree/",patient,".maf"),header = T)
  tumorNames = colnames(maf_data)[-which(colnames(maf_data) %in% c("mut_id","Variant_Classification"))]
  maf_data$N = rep(0, nrow(maf_data))
  count_mat = maf_data %>% dplyr::select(tumorNames, "N") %>% mutate_at(c(tumorNames, "N"), as.numeric)
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
  
  dat=read.table(paste0("../1/tree/",patient,"_DESeq2.log2_MADtop3000CodingGenes.txt"),header = T)
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
      #ITH=c(ITH,mean(apply(pair,2,function(x) cal_ITH(dist_gene,x))))
      rna_ITH=c(rna_ITH,mean(apply(pair,2,function(x) cal_ITH(rna_dist,x))))
    }
  }
  
  
  
  all_point_ITH=cbind(all_point[f,],rna_ITH)
  # cluster=readRDS("dna_cluster.RDS")
  rna_cluster=readRDS("../2/rna_cluster.RDS")
  #class=as.character(cluster[[patient]][select_info$Sample])
  rna_class=as.character(rna_cluster[[patient]][select_info$Sample])
  if (length(rna_class>0)){rna_class=rna_class}else{rna_class="1"}
  #class=read.table("../3/class.tsv")
  rna_subtype=class[substr(class$V1,1,4)==patient,]%>%arrange(.,V1)%>%pull(V2)%>%as.character()
  #if (length(class>0)){class=class}else{class="1"}
  #if (length(rna_class>0)){rna_class=rna_class}else{rna_class="1"}
  select_info=cbind(select_info,rna_class,rna_subtype)
  #select_info$class=as.character(select_info$class)
  select_info$rna_class=as.character(select_info$rna_class)
  select_info$rna_subtype=as.character(select_info$rna_subtype)
  select_info$rna_subtype=sub("1","Wnt",select_info$rna_subtype)
  select_info$rna_subtype=sub("2","Classic",select_info$rna_subtype)
  select_info$rna_subtype=sub("3","Hybrid",select_info$rna_subtype)
  select_info$rna_subtype=sub("4","Poliferation",select_info$rna_subtype)
    
  color=c(RColorBrewer::brewer.pal(7,"Set2")[c(5,1,2,4)])
  names(color)=c("Wnt","Classic","Hybrid","Poliferation")

  r=c()
  for (n in 1:length(rna_ITH)){
    if (sample(c(0,1),1,prob=c(1-rna_ITH[n],rna_ITH[n]))==1){
      r=c(r,n)
    }
    print(n)
  }
  all_point1_rna=all_point_ITH[r,]
  

  p4_rna=ggplot(data=all_point1_rna,aes(x=x,y=y))+
    stat_density_2d(geom = "raster",aes(fill = ..ndensity..),contour = F)+
    #stat_density_2d(geom = "polygon",aes(fill=..level..),bins=80)+
    xlab("X (cm)")+ylab("Y (cm)")+
    ggtitle(paste0(patient_paper," spatial transcriptomic heterogeneity"))+labs(fill="ITH")+
    scale_fill_gradientn(labels=round(seq(0,max(all_point_ITH$rna_ITH),
                                          length.out=5),digits = 2),
                         breaks=c(0,0.25,0.5,0.75,1),colours =  colormap)+
    theme_classic()+ggnewscale::new_scale_fill()+
    geom_point(data = select_info,aes(x=X,y=Y,fill=rna_subtype),size=2,shape=21,stroke=0.4)+
    scale_fill_manual(values =  color[unique(select_info$rna_subtype)])+
    labs(fill="Subtype",shape="")+
    #scale_shape_manual(values = c(21,22,23))+
    # geom_text(data = select_info,color="white",aes(x=X,y=Y,
    #                                  label=substr(Sample_paper,6,7)))+
    guides(fill = guide_legend(order=1))+
    theme(text=element_text(size = 6.5),
          legend.spacing = unit(0,"cm"),
          plot.title=element_text(hjust = 0.5,size=7),
          legend.key.width = unit(0.1, 'cm'),
          legend.key.height = unit(0.3, 'cm'),
          line = element_line(size = 0.1),
          legend.box.spacing = unit(0.01,"cm"),
          axis.ticks.length = unit(0.03,"cm"))
  
  a=MASS::kde2d(all_point1_rna$x,all_point1_rna$y,n=100)
  data_density=matrix(ncol=3)
  for (n in 1:length(a$x)){
    data_density=rbind(data_density,cbind(a$z[n,],a$x[n],a$y))
  }
  data_density=as.data.frame(data_density)
  data_density=na.omit(data_density)
  data_density1=data_density[which(data_density$V1>quantile(data_density$V1)),]
  sample_plot[[patient]]=p4_rna
  #ggarrange(p2,p4,p1,p3,p2_rna,p4_rna,p1_rna,p3_rna,ncol = 4,nrow = 2)
  #ggsave(paste0("ITH_4",patient_paper,".pdf"),width = 8,height = 4)
}


set_absolute_pos <- function(segs, cyto_length){
  chrom <- c(seq(1,22),"X")
  # Add cumulative length and get absolute position
  cum_length <- 0
  for (chr in chrom){
    segs[which(segs$chromosome == chr), 'Start'] <- segs[which(segs$chromosome == chr), 'Start'] + cum_length
    segs[which(segs$chromosome == chr), 'End'] <- segs[which(segs$chromosome == chr), 'End'] + cum_length
    cum_length <- cum_length + cytoband_length[which(cytoband_length$chromosome == chr), 'length', drop=T]
  }
  return(segs)
}


cytoband_length <- cyto%>% group_by(X1) %>% summarise(length = max(X3))
cytoband_length$chromosome=sub("chr","",cytoband_length$X1)
all_segs=all_seg%>%dplyr::rename(End="end.pos")%>%dplyr::rename(Start="start.pos")
all_segs=set_absolute_pos(all_segs,genome_len)
cytoband_pos <- cytoband_length
cytoband_pos <- cytoband_pos %>% arrange(factor(cytoband_pos$chromosome, levels = seq(1,22)))
# Get absolute position of chrom end
cytoband_pos <- cytoband_pos %>% mutate(pos=cumsum(length))
cytoband_arm <- cyto %>% dplyr::filter(X5=="acen")
colnames(cytoband_arm)=c("chromosome","Start","End","X4","X5")
cytoband_arm$chromosome=sub("chr","",cytoband_arm$chromosome)
cytoband_arm <- set_absolute_pos(cytoband_arm, cytoband_length)

# Set label for p and q
cytoband_arm <- cytoband_arm %>% group_by(chromosome) %>% summarise(Start=min(Start), End=max(End))
cytoband_arm <- cytoband_arm %>% mutate(start_label = paste0(chromosome, "p"),
                                        end_label = paste0(chromosome, "q"))

# Set label for mid point of p and q
cytoband_arm <- cytoband_arm %>% mutate(midpoint=(Start + End)/2)
cytoband_arm <- cytoband_arm %>% mutate(size=End-Start)

cytoband_pos$chromosome
all_segs$class=class[substr(all_segs$sample,1,7),1]
library(scales)
plot_genome=ggplot(data=all_segs) +
  geom_segment(aes(x=sample, xend=sample,
                   y=Start, yend=End, 
                   color=depth.ratio), size=4) +
  geom_hline(yintercept=cytoband_pos$pos,size=0.01,col="dimgray") +
  #geom_hline(yintercept = cytoband_arm$midpoint, alpha=0.5, linetype=2) +
  scale_y_continuous(trans = "reverse",breaks=as.numeric(cytoband_arm$midpoint),
                     labels = cytoband_arm$chromosome,expand = c(0, 0)) +
  scale_color_gradient2(low="#4393c3", high="#d6604d",mid="#f7f7f7",
                        name="Adjusted CN",limits=c(-1,1),oob=squish)+
  theme_bw() + 
  labs(x = "", y="") + 
  theme(axis.text.y = element_text(size = 8,vjust =0.5,hjust=1),panel.grid = element_blank(),axis.ticks = element_blank(),
        axis.text.x = element_blank(),axis.ticks.x = element_blank(),legend.direction = 'vertical',
        legend.position = 'right')+
  facet_grid(~class,scales='free')

all_segs_DT10=all_segs[substr(all_segs$sample,1,4) =="DT10",]
all_segs_DT10$class=sub(2,"Classic",all_segs_DT10$class)
all_segs_DT10$class=sub(4,"Poliferation",all_segs_DT10$class)
library(RColorBrewer)
plot_DT10=ggplot(data=all_segs_DT10) +
  geom_segment(aes(x=sample, xend=sample,
                   y=Start, yend=End, 
                   color=depth.ratio), size=5) +
  geom_hline(yintercept=cytoband_pos$pos,size=0.01,col="dimgray") +
  #geom_hline(yintercept = cytoband_arm$midpoint, alpha=0.5, linetype=2) +
  scale_y_continuous(trans = "reverse",breaks=as.numeric(cytoband_arm$midpoint),
                     labels = cytoband_arm$chromosome,expand = c(0, 0)) +
  scale_color_gradient2(low="#4393c3", high="#d6604d",mid="#f7f7f7",
                        name="Adjusted CN",limits=c(-0.25,0.25),oob=squish) +
  theme_bw() + 
  labs(x = "", y="") + ggtitle("SH06")+
  theme(axis.text.y = element_text(size = 3.5,vjust =0.5,hjust=1),
        panel.grid = element_blank(),axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.border = element_rect(size=0.2),
        legend.direction = 'vertical',
        legend.spacing = unit(0,"cm"),
        legend.key.width = unit(0.1, 'cm'),
        legend.key.height = unit(0.3, 'cm'),
        legend.box.spacing = unit(0.01,"cm"),
        legend.title = element_text(size=6),
        legend.text = element_text(size=5),
        legend.position = 'right',
        plot.title = element_text(size=6,hjust=0.5),
        strip.text = element_text(size=5),
        strip.background = element_rect(size=0.2),
        plot.margin = unit(c(0,0,0,0),"mm"),
        panel.spacing = unit(0.2,"mm"))+
  facet_grid(~class,scales='free')



all_segs_DT18=all_segs[substr(all_segs$sample,1,4) =="DT18",]
plot_DT18=ggplot(data=all_segs_DT18) +
  geom_segment(aes(x=sample, xend=sample,
                   y=Start, yend=End, 
                   color=depth.ratio), size=5) +
  geom_hline(yintercept=cytoband_pos$pos,size=0.01,col="dimgray") +
  #geom_hline(yintercept = cytoband_arm$midpoint, alpha=0.5, linetype=2) +
  scale_y_continuous(trans = "reverse",breaks=as.numeric(cytoband_arm$midpoint),
                     labels = cytoband_arm$chromosome,expand = c(0, 0)) +
  scale_color_gradient2(low = "#2600D1", high = "#D60C00",mid="white",midpoint = 0,limits=c(-1.5,1.5)) +
  theme_bw() + 
  labs(x = "", y="") + ggtitle("SH06")+
  theme(axis.text.y = element_text(size = 3.5,vjust =0.5,hjust=1),
        panel.grid = element_blank(),axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.border = element_rect(size=0.2),
        legend.direction = 'vertical',
        legend.spacing = unit(0,"cm"),
        legend.key.width = unit(0.1, 'cm'),
        legend.key.height = unit(0.3, 'cm'),
        legend.box.spacing = unit(0.01,"cm"),
        legend.title = element_blank(),
        legend.text = element_text(size=5),
        legend.position = 'right',
        plot.title = element_text(size=6,hjust=0.5),
        strip.text = element_text(size=5),
        strip.background = element_rect(size=0.2),
        plot.margin = unit(c(0,0,0,0),"mm"),
        panel.spacing = unit(0.2,"mm"))+
  facet_grid(~class,scales='free')



exp_counts = read.table("E:/cancer genome/liver/分析/DEseq2-normalized_clustering/DT03-DT19//outputs/S01_expCounts_MeanDupCodingGenes.txt",header=TRUE)
exp_counts<-exp_counts[,which(substr(colnames(exp_counts),1,4)=="DT18")]
exp_counts<-exp_counts[,-c(6)]
coldata = data.frame(colnames(exp_counts))


coldata$condition=c(rep("CLASS01", 2), rep("CLASS02", 3))
dds = DESeqDataSetFromMatrix(countData=round(exp_counts), 
                             colData=coldata, design=~condition) # ~1 means no design
dds <- DESeq(dds)
#for (i in c("03")){
res = results(dds, contrast=c("condition","CLASS01", "CLASS02"))
res = res[order(res$pvalue),]
head(res)
summary(res)
#write.csv(res,file="All_results.csv")
diff_gene_deseq2_up <-subset(res,padj < 0.05 & (log2FoldChange > 1 ))
dim(diff_gene_deseq2_up)
#write.csv(diff_gene_deseq2,file= "tumor_normal_up.csv")

diff_gene_deseq2_down <-subset(res,padj < 0.05 & (log2FoldChange < -1))
dim(diff_gene_deseq2_down)
# write.csv(diff_gene_deseq2,file= "tumor_normal_down.csv")
#head(diff_gene_deseq2)
x_up=rownames(diff_gene_deseq2_up)
up = bitr(x_up, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
x_down=rownames(diff_gene_deseq2_down)
down = bitr(x_down, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# DT18_up <- enrichGO(up$ENTREZID, "org.Hs.eg.db", keyType = "ENTREZID",ont = 'ALL',
#                     pvalueCutoff  = 0.05,pAdjustMethod = "BH",  qvalueCutoff  = 0.1, readable=T)
DT18_up<- enrichKEGG(up$ENTREZID, "hsa",
                        pvalueCutoff  = 0.05,pAdjustMethod = "BH",  qvalueCutoff  = 0.1)

result1=DT18_up@result[1:7,]
result1$qvalue=-log(result1$qvalue,10)

# DT18_down <- enrichGO(down$ENTREZID, "org.Hs.eg.db", keyType = "ENTREZID",ont = 'ALL',
#                       pvalueCutoff  = 0.05,pAdjustMethod = "BH",  qvalueCutoff  = 0.1, readable=T)
DT18_down <- enrichKEGG(down$ENTREZID, "hsa",
                      pvalueCutoff  = 0.05,pAdjustMethod = "BH",  qvalueCutoff  = 0.1)

result2=DT18_down@result[1:5,]
result2$qvalue=log(result2$qvalue,10)

result1$fill="red"
result2$fill="blue"
result_DT18=rbind(result1,result2)
result_DT18$Description=factor(c(result1$Description,
                                 result2$Description),
                               levels = c(rev(result2$Description),
                                          rev(result1$Description)
                               ))
library(ggbreak)
GO_DT18=ggplot(result_DT18)+geom_bar(aes(x=Description,y=qvalue,fill=fill),
                                     stat = "identity",width = 0.3)+
  scale_fill_manual(values =  RColorBrewer::brewer.pal(n = 3,"Set2"))+
  scale_y_continuous(expand = c(0,0),
                     breaks=c(-5,-1,0,1,10,20),labels = c(5,1,0,1,10,20))+
  #scale_y_break(c(5,18),scales = 0.5,ticklabels=c(20,30,40))+
  coord_flip()+
  scale_x_discrete()+
  xlab("")+ggtitle("KEGG enrichment for SH12")+
  ylab("Enrichment significance (-log10 (q))")+theme_classic()+
  #geom_hline(aes(yintercept=1),color="red",linetype=2)+
  geom_hline(aes(yintercept=0),size=0.1)+
  #geom_hline(aes(yintercept=-1),color="red",linetype=2)+
  theme(panel.background = element_blank(),
        legend.position = "none",
        plot.title =  element_text(size=6,hjust = 0.5),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        #axis.title.x.bottom  = element_blank(),
        axis.title.x = element_text(size=6,hjust = 0.8),
        axis.text = element_text(size = 6),
        axis.line = element_line(size = 0.3),
        axis.ticks = element_line(size=0.3))
pdf("p4_v12.pdf",width = 8,height = 6)
p1=ggarrange(sample_plot[["DT09"]],p_go,sample_plot[["DT16"]],
             widths = c(1,1.5,1), ncol = 3,labels=c("a","b","c"),
             font.label = list(size=8))
p2=ggarrange(sample_plot[["DT10"]],GO_DT10,plot_DT10,
             widths = c(1,1.3,1.2), ncol = 3,labels=c("d","e","f"),
             font.label = list(size=8))
p3=ggarrange(p_10_GII,tree_DT18_2,GO_DT18,widths = c(1,1,1.5), ncol = 3,
             labels=c("g","h","i"),font.label = list(size=8))

ggarrange(p1,p2,p3,nrow = 3)
dev.off()
# save.image(file="p4.Rdata")
# pdf("subtype_position.pdf",width = 2.36,height = 1.97)
# print(sample_plot)
# dev.off()

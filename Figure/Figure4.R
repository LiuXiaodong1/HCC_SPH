library(tidyverse)
library(ggpubr)
library(DESeq2)
library(clusterProfiler)
library(org.Hs.eg.db)

all_seg=read_tsv("../Data/Figure4g.tsv")
class1=c("SH06T01", "SH06T02", "SH06T03", "SH06T06",
         "SH06T07", "SH06T08", "SH06T11", "SH06T12")
class2=c("SH06T04", "SH06T05", "SH06T09", "SH06T10",
         "SH06T13", "SH06T14", "SH06T15")

# cyto=read_tsv("../data/cytoBand.txt",col_names = F)
# cyto=cyto%>%filter(X1%in%c(paste0("chr",c(1:22,"X"))))
# genome_len=sum(cyto$X3-cyto$X2)
genome_len=3036303846

GII1=c()
for (s in class1){
  tmp=all_seg%>%filter(sample==s,depth.ratio>0.2|depth.ratio< -0.2)
  len=sum(tmp$end.pos-tmp$start.pos)
  GII1=c(GII1,len/genome_len)
}

GII2=c()
for (s in class2){
  tmp=all_seg%>%filter(sample==s,depth.ratio>0.2|depth.ratio< -0.2)
  len=sum(tmp$end.pos-tmp$start.pos)
  GII2=c(GII2,len/genome_len)
}

mycompare=list(c(1,2))
GII1=cbind.data.frame(GII=GII1,class="GII1")
GII2=cbind.data.frame(GII=GII2,class="GII2")
GII=rbind(GII1,GII2)
GII$class=factor(GII$class,level=c("GII2","GII1"))

class_number=cbind.data.frame(class=c("GII1","GII2"),n=c(8,7),y=0)
Fig4g=ggplot(data = GII,aes(class,GII))+
  geom_boxplot(aes(fill=class),outlier.alpha = 0,linewidth=0.1)+
  geom_text(data=class_number,
            aes(x=class,y=0,label=paste0("n=",n)),
            size=1.5)+
  geom_point(size=0.5)+scale_x_discrete(labels=c("Classic","Poliferation"))+
  xlab("Subtype")+ggtitle("SH06")+ggplot2::ylim(0,0.62)+
  scale_fill_manual(values=c(RColorBrewer::brewer.pal(7,"Set2")[c(1,4)]))+
  theme_classic()+
  stat_compare_means(
    comparisons = mycompare,
    method = "wilcox.test",size=2)+
  theme(text=element_text(size = 6.5),
        legend.position = "none",
        plot.title=element_text(hjust = 0.5,size=7),
        line = element_line(size = 0.1),
        axis.ticks.length = unit(0.03,"cm"))

exp_counts = read.table("E:/cancer genome/liver/analysis/DEseq2-normalized_clustering/DT03-DT17/outputs/S01_expCounts_MeanDupCodingGenes.txt",header=TRUE)

exp_counts<-exp_counts[,which(substr(colnames(exp_counts),1,4)=="DT09")]
exp_counts<-exp_counts[,-c(29,30)]
coldata = data.frame(colnames(exp_counts))

coldata$condition=c(rep("CLASS01", 9), rep("CLASS02", 2),
                    rep("CLASS01", 5), rep("CLASS02", 2),
                    rep("CLASS01", 11),  rep("CLASS02", 1)
)
dds = DESeqDataSetFromMatrix(countData=round(exp_counts), 
                             colData=coldata, design=~condition) 
dds <- DESeq(dds)
res = results(dds, contrast=c("condition", "CLASS01","CLASS02"))
res = res[order(res$pvalue),]

diff_gene_deseq2_up <-subset(res,padj < 0.05 & (log2FoldChange > 1 ))

diff_gene_deseq2_down <-subset(res,padj < 0.05 & (log2FoldChange < -1))
x_up=rownames(diff_gene_deseq2_up)
up = bitr(x_up, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
x_down=rownames(diff_gene_deseq2_down)
down = bitr(x_down, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
DT09_up <- enrichGO(up$ENTREZID, "org.Hs.eg.db", keyType = "ENTREZID",ont = 'ALL',
                    pvalueCutoff  = 0.05,pAdjustMethod = "BH",  
                    qvalueCutoff  = 0.1, readable=T)
DT09_down <- enrichGO(down$ENTREZID, "org.Hs.eg.db", keyType = "ENTREZID",ont = 'ALL',
                      pvalueCutoff  = 0.05,pAdjustMethod = "BH",  qvalueCutoff  = 0.1, 
                      readable=T)


exp_counts = read.table("E:/cancer genome/liver/analysis/DEseq2-normalized_clustering/DT03-DT17/outputs/S01_expCounts_MeanDupCodingGenes.txt",header=TRUE)
exp_counts<-exp_counts[,which(substr(colnames(exp_counts),1,4)=="DT10")]
exp_counts<-exp_counts[,-c(16,17)]
coldata = data.frame(colnames(exp_counts))


coldata$condition=c(rep("CLASS01", 3), rep("CLASS02", 2),
                    rep("CLASS01", 3), rep("CLASS02", 2),
                    rep("CLASS01", 2), rep("CLASS02", 3)
)
dds = DESeqDataSetFromMatrix(countData=round(exp_counts), 
                             colData=coldata, design=~condition) 
dds <- DESeq(dds)
res = results(dds, contrast=c("condition","CLASS01", "CLASS02"))
res = res[order(res$pvalue),]
diff_gene_deseq2_up <-subset(res,padj < 0.05 & (log2FoldChange > 1 ))
diff_gene_deseq2_down <-subset(res,padj < 0.05 & (log2FoldChange < -1))
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

Fig4e=ggplot(result_DT10)+geom_bar(aes(x=Description,y=qvalue,fill=fill),
                                     stat = "identity",width = 0.3)+
  scale_fill_manual(values =  RColorBrewer::brewer.pal(n = 3,"Set2"))+
  scale_y_continuous(expand = c(0,0),
                     breaks=c(-5,0,10,20),labels = c(5,0,10,20))+
  scale_y_break(c(5,18),scales = 0.5,ticklabels=c(20,30,40))+
  coord_flip()+
  scale_x_discrete()+
  xlab("")+ggtitle("GO enrichment for SH06")+
  ylab("Enrichment significance (-log10 (q))")+theme_classic()+
  geom_hline(aes(yintercept=0),size=0.1)+
  theme(panel.background = element_blank(),
        legend.position = "none",
        plot.title =  element_text(size=6,hjust = 0.5),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size=6,hjust = 0.8),
        axis.text = element_text(size = 6),
        axis.line = element_line(size = 0.3),
        axis.ticks = element_line(size=0.3))

exp_counts = read.table("E:/cancer genome/liver/analysis/DEseq2-normalized_clustering/DT03-DT19/outputs/S01_expCounts_MeanDupCodingGenes.txt",header=TRUE)
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
res = results(dds, contrast=c("condition","CLASS01", "CLASS02"))
res = res[order(res$pvalue),]
diff_gene_deseq2_up <-subset(res,padj < 0.05 & (log2FoldChange > 1 ))
diff_gene_deseq2_down <-subset(res,padj < 0.05 & (log2FoldChange < -1))
x_up=rownames(diff_gene_deseq2_up)
up = bitr(x_up, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
x_down=rownames(diff_gene_deseq2_down)
down = bitr(x_down, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
DT16_up <- enrichGO(up$ENTREZID, "org.Hs.eg.db", keyType = "ENTREZID",ont = 'ALL',
                    pvalueCutoff  = 0.05,pAdjustMethod = "BH",  qvalueCutoff  = 0.1, readable=T)
DT16_down <- enrichGO(down$ENTREZID, "org.Hs.eg.db", keyType = "ENTREZID",ont = 'ALL',
                      pvalueCutoff  = 0.05,pAdjustMethod = "BH",  qvalueCutoff  = 0.1, readable=T)

SH05_up=DT09_up@result%>%arrange(qvalue)
path_up1=SH05_up%>%head(5)%>%pull(Description)
SH16_up=DT16_up@result%>%arrange(qvalue)
path_up2=SH16_up%>%head(5)%>%pull(Description)
path_up=c(path_up1,path_up2)
tmp=SH05_up%>%filter(Description%in%path_up)%>%
  dplyr::select(Description,qvalue)%>%
  mutate(sample="SH05",fill="blue",qvalue= -log10(qvalue))
tmp1=SH16_up%>%filter(Description%in%path_up)%>%
  dplyr::select(Description,qvalue)%>%
  mutate(sample="SH10",fill="red",qvalue= -log10(qvalue))
df=rbind(tmp,tmp1)
order1=df%>%group_by(Description)%>%summarise(qvalue=max(qvalue))%>%
  arrange(desc(qvalue))%>%head(5)%>%pull(Description)
SH05_down=DT09_down@result%>%arrange(qvalue)
path_down1=SH05_down%>%head(5)%>%pull(Description)
SH16_down=DT16_down@result%>%arrange(qvalue)
path_down2=SH16_down%>%head(5)%>%pull(Description)
path_down=c(path_down1,path_down2)
tmp=SH05_down%>%filter(Description%in%path_down)%>%
  dplyr::select(Description,qvalue)%>%
  mutate(sample="SH05",fill="blue",qvalue= log10(qvalue))
tmp1=SH16_down%>%filter(Description%in%path_down)%>%
  dplyr::select(Description,qvalue)%>%
  mutate(sample="SH10",fill="red",qvalue= log10(qvalue))
df1=rbind(tmp,tmp1)
order2=df1%>%group_by(Description)%>%summarise(qvalue=min(qvalue))%>%
  arrange(qvalue)%>%head(5)%>%pull(Description)
order=c(rev(order2),rev(order1))
df=rbind(df,df1)%>%filter(Description%in%order)
df$Description=factor(df$Description,levels =order)
Fig4b=ggplot(df)+geom_bar(aes(x=Description,y=qvalue,fill=sample),
                         stat = "identity",width = 0.3,position = position_dodge())+
  scale_fill_manual(values =c(RColorBrewer::brewer.pal(12,"Set3")[c(5,6)]))+
  scale_y_continuous(expand = c(0,0),
                     breaks=c(-50,-25,0,25),
                     labels = c(50,25,0,25))+
  coord_flip()+
  scale_x_discrete()+
  xlab("")+ggtitle("GO enrichment for SH05/SH10")+
  theme_classic()+
  geom_hline(aes(yintercept=0),size=0.1)+
  ylab("Enrichment significance (-log10 (q))")+
  theme(panel.background = element_blank(),
        legend.spacing = unit(0,"cm"),
        legend.key.width = unit(0.1, 'cm'),
        legend.key.height = unit(0.3, 'cm'),
        legend.box.spacing = unit(0.01,"cm"),
        legend.title = element_blank(),
        legend.text = element_text(size=5),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(size=6,hjust = 0.8),
        axis.title.x = element_text(size=6,hjust = 0.8),
        axis.text = element_text(size = 6),
        axis.line = element_line(size = 0.3),
        axis.ticks = element_line(size=0.3))

exp_counts = read.table("E:/cancer genome/liver/analysis/DEseq2-normalized_clustering/DT03-DT19//outputs/S01_expCounts_MeanDupCodingGenes.txt",header=TRUE)
exp_counts<-exp_counts[,which(substr(colnames(exp_counts),1,4)=="DT18")]
exp_counts<-exp_counts[,-c(6)]
coldata = data.frame(colnames(exp_counts))
coldata$condition=c(rep("CLASS01", 2), rep("CLASS02", 3))
dds = DESeqDataSetFromMatrix(countData=round(exp_counts), 
                             colData=coldata, design=~condition) 
dds <- DESeq(dds)
res = results(dds, contrast=c("condition","CLASS01", "CLASS02"))
res = res[order(res$pvalue),]
diff_gene_deseq2_up <-subset(res,padj < 0.05 & (log2FoldChange > 1 ))
diff_gene_deseq2_down <-subset(res,padj < 0.05 & (log2FoldChange < -1))
x_up=rownames(diff_gene_deseq2_up)
up = bitr(x_up, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
x_down=rownames(diff_gene_deseq2_down)
down = bitr(x_down, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

DT18_up<- enrichKEGG(up$ENTREZID, "hsa",
                     pvalueCutoff  = 0.05,pAdjustMethod = "BH",  qvalueCutoff  = 0.1)

result1=DT18_up@result[1:7,]
result1$qvalue=-log(result1$qvalue,10)
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
Fig4i=ggplot(result_DT18)+geom_bar(aes(x=Description,y=qvalue,fill=fill),
                                     stat = "identity",width = 0.3)+
  scale_fill_manual(values =  RColorBrewer::brewer.pal(n = 3,"Set2"))+
  scale_y_continuous(expand = c(0,0),
                     breaks=c(-5,-1,0,1,10,20),labels = c(5,1,0,1,10,20))+
  coord_flip()+
  scale_x_discrete()+
  xlab("")+ggtitle("KEGG enrichment for SH12")+
  ylab("Enrichment significance (-log10 (q))")+theme_classic()+
  geom_hline(aes(yintercept=0),size=0.1)+
  theme(panel.background = element_blank(),
        legend.position = "none",
        plot.title =  element_text(size=6,hjust = 0.5),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size=6,hjust = 0.8),
        axis.text = element_text(size = 6),
        axis.line = element_line(size = 0.3),
        axis.ticks = element_line(size=0.3))

Fig4i

Fig4b

Fig4g

Fig4e

library(tidyverse)
library(ComplexHeatmap)
library(RColorBrewer)
library(survival)
library(survminer)
library(ConsensusClusterPlus)
library(circlize)
library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)
library(ggpubr)
library(ggsci)
library(dendextend)
library(colorRamps)
set.seed(12345)


changeP <- function(p_val,r2) {
  tmp=strsplit(as.character(p_val),"e")[[1]]
  if(length(tmp)>1){
    a=tmp[1];b=tmp[2]
    as.expression(bquote(atop(p~"="~.(a)~"??"~10^~.(b)~","~R^2~"="~.(r2))))
  }else{
    as.expression(bquote(atop(p~"="~.(p_val)~","~R^2~"="~.(r2))))
  }
}
marker=read.table("E:/cancer genome/liver/analysis//decoder/lihc/HCC_primary/Final_marker_genes.txt",sep = "\t",header = T)
# colnames(marker)=gsub("_.*","",colnames(marker))%>%sub("X5.5","X5.5(Wnt)",.)%>%sub("X5.1","X5.1(metabolism)",.)%>%
#   sub("X5.3","X5.3(ECM)",.)%>%sub("X7.5","X7.5(Cell Cycle)",.)
colnames(marker)=gsub("_.*","",colnames(marker))%>%sub("X5.5","Wnt",.)%>%
  sub("X5.1","metabolism",.)%>%
  sub("X5.3","ECM",.)%>%
  sub("X7.5","Cell Cycle",.)
marker=marker[1:100,]
survival=read_tsv("E:/cancer genome/liver/analysis//RNA_cluster/NMF/lihc_tcga_clinical_data.tsv")
survival$`Disease Free Status`=as.numeric(substr(survival$`Disease Free Status`,1,1))
survival$`Overall Survival Status`=as.numeric(substr(survival$`Overall Survival Status`,1,1))
gene_count=read_tsv("E:/cancer genome/liver/analysis//decoder/lihc/HCC_primary/deseq_log2_primary.tsv")
#colnames(gene_count)=gsub("\\.","-",colnames(gene_count))


other_studies=read.table("../../4/literatre_pathway_tumor.tsv")
other_studies=other_studies[colnames(gene_count)[-1],]

maffile=list.files(path = "..",pattern = "*.maf$",all.files = T,
                   recursive = T,full.names = T)
maf1=read_tsv(maffile[1],comment = "#")
maf2=read_tsv(maffile[2],comment = "#")
maf3=read_tsv(maffile[3],comment = "#")
maf4=read_tsv(maffile[4],comment = '#')
maf=dplyr::bind_rows(maf1[,c(1,16)],maf2[,c(1,16)],maf3[,c(1,16)],maf4[,c(1,16)])
maf=unique(maf)
maf=maf%>%filter(Hugo_Symbol=="CTNNB1")
mut_sample=gsub(pattern = "-",replacement = ".",
                substr(maf$Tumor_Sample_Barcode,1,16))
CTNNB1_mut = colnames(gene_count)[-1]%in%mut_sample

ha_fun_s=function(other_studies){
  anno = HeatmapAnnotation(Lee_2004=other_studies$Lee,
                           Boyault_2007=other_studies$Boyault,
                           Chiang_2008=other_studies$Chiang,
                           Hoshida_2009=other_studies$Hoshida,
                           "CTNNB1 mutation"=CTNNB1_mut,
                           show_annotation_name = T,show_legend = F,
                           annotation_name_side = "right",
                           annotation_name_gp= gpar(fontsize =5),
                           annotation_legend_param = list(title_gp = gpar(fontsize =5),
                                                          labels_gp = gpar(fontsize =5),
                                                          grid_width = unit(1, "mm"),
                                                          grid_height=unit(2,"mm"),
                                                          legend_height = unit(1, "cm")),
                           simple_anno_size  = unit(2,"mm"),border=F,
                           col = list(Lee_2004=c("SURVIVAL_DN"="hotpink4","SURVIVAL_UP"="khaki3"),
                                      Hoshida_2009=c("S1"=brewer.pal(8,"Set1")[3],"S2"=brewer.pal(8,"Set1")[4],"S3"=brewer.pal(8,"Set1")[5]),
                                      Chiang_2008=c("CTNNB1"=brewer.pal(5,"Accent")[1],"INTERFERON"=brewer.pal(5,"Accent")[2],
                                                    "POLYSOMY7"=brewer.pal(5,"Accent")[3],"PROLIFERATION"=brewer.pal(5,"Accent")[4],
                                                    "UNANNOTATED"=brewer.pal(5,"Accent")[5]),
                                      "CTNNB1 mutation"=c("TRUE"="brown1","FALSE"="darkgreen"),
                                      Boyault_2007=c("G12"="pink","G3"="beige","G56"="dodgerblue")
                           ))
  return(anno)
}
purity=read.table("E:/cancer genome/liver/analysis/RNA_cluster/gene_symbol/TCGA_data.tsv",header = T)
purity=purity%>%dplyr::select(tumor_sample_barcode,purity)%>%unique()
purity$tumor_sample_barcode=gsub("-",".", purity$tumor_sample_barcode)
rownames(purity)=purity$tumor_sample_barcode
sample_select=intersect(rownames(purity),colnames(gene_count))
weight=read.table("E:/cancer genome/liver/analysis/decoder/lihc/HCC_primary/Final_sample_weights.txt",header = T)
rownames(weight)=weight$sampleID
colnames(weight)=gsub("_.*","",colnames(weight))
weight_s=weight[sample_select,]
purity=purity[sample_select,]
data=cbind.data.frame(weight_s,purity)
data$tumor_sample_barcode=NULL
pic=NULL;factor_name=c("Wnt signaling","Metabolism","ECM","Cell cycle")
for (i in 2:(ncol(data)-1)){
  s = summary(lm(data[,ncol(data)]~data[,i]))
  r2 =signif(s$r.squared,3)
  p_val=signif(s$coefficients[8],3)
  t=as.data.frame(cbind(r2,p_val))
  t$x=(max(data$purity)+min(data$purity))/2;t$y=max(data[,i]+0.5)
  pic[[i-1]]=ggplot(data = data,aes_string(x="purity",y=colnames(weight)[i]))+
    geom_point(color="grey")+
    geom_smooth(method='lm',se=FALSE,colour="red",linetype="dashed")+
    labs(subtitle = changeP(p_val,r2))+
    xlab("Tumor cellularity")+ylab(paste0(factor_name[i-1]," factor weight"))+
    theme_classic()+ylim(c(min(data[,i]-0.25),max(data[,i]+0.5)))+
    theme(text=element_text(size = 6.5),
          legend.spacing = unit(0,"cm"),
          plot.title=element_text(hjust = 0.5),
          legend.key.width = unit(0.1, 'cm'),
          legend.key.height = unit(0.3, 'cm'),
          line = element_line(size = 0.1),
          legend.box.spacing = unit(0.01,"cm"),
          axis.ticks.length = unit(0.03,"cm"),
          plot.subtitle=element_text(size=5, hjust=0.5, face="italic", color="black"))
}
col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

plot_fun=function(marker,data,dis="euclidean",ncluster=2,n_col=2){
  label=stringr::str_to_sentence(colnames(marker))
  label=ifelse(label=="Ecm","ECM",label)
  select_data= data%>%filter(gene_id %in% c(unlist(marker)))
  id=select_data$gene_id
  select_data$gene_id=NULL
  rownames(select_data)=id
  
  gene=cbind.data.frame(gene=c(unlist(marker)),
                        class=rep(colnames(marker),each= nrow(marker)))
  gene=gene[!duplicated(gene$gene),]
  rownames(gene)=gene$gene
  gene=gene[rownames(select_data),]
  #marker_col=c(rainbow(ncol(marker)))
  marker_col=ggsci::pal_d3("category20")(ncol(marker))
  names(marker_col)=colnames(marker)
  gene_anno = rowAnnotation(factor = gene$class,
                            show_annotation_name=F,
                            col=list(factor=marker_col))
  
  gene_weight_scale=t(scale(t(select_data)))
  
  result_Sp_wd=ConsensusClusterPlus(gene_weight_scale,maxK=3,reps=500,pItem=0.8,
                                    seed = 123,
                                    pFeature=0.8,plot = "pdf",
                                    innerLinkage="ward.D2",
                                    distance=dis)
  ha=ha_fun_s(other_studies )
  
  cluster=cutree(result_Sp_wd[[ncluster]]$consensusTree,k = 2)
  if (mean(gene_weight_scale[,which(cluster==1)])>
      mean(gene_weight_scale[,which(cluster==2)])){
    cluster[which(cluster==1)]="High";cluster[which(cluster==2)]="Low"
  }else{
    cluster[which(cluster==1)]="Low";cluster[which(cluster==2)]="High"
  }
  cluster_anno=HeatmapAnnotation(`Metabolism cluster`=cluster,
                                 simple_anno_size  = unit(2,"mm"),
                                 annotation_name_gp= gpar(fontsize =6),
                                 show_legend = F,show_annotation_name = F,
                                 annotation_legend_param = list(title_gp = gpar(fontsize =5),
                                                                labels_gp = gpar(fontsize =5),title="",
                                                                grid_width = unit(1, "mm"),
                                                                grid_height=unit(2,"mm"),
                                                                legend_height = unit(1, "cm")
                                 ),col = list("Metabolism cluster"=c("High"="orange",
                                                                     "Low"="lightblue")))
  ord=order.dendrogram(as.dendrogram(result_Sp_wd[[ncluster]]$consensusTree))
  if(cluster[ord][1]=="Low"){
    r=2
  }else{
    r=1
  }
  column_tree=color_branches(as.dendrogram(result_Sp_wd[[ncluster]]$consensusTree),
                             k=2,
                             col = c("1"=ifelse(r==1,"orange","lightblue"),
                                     "2"=ifelse(r==1,"lightblue","orange")))
  
  
  p=Heatmap(gene_weight_scale,show_row_names = F,show_row_dend = F,
            row_title = label,row_title_side = "left",
            column_dend_side = "bottom",
            col=col_fun,name = "Z-score",
            cluster_columns = column_tree,
            clustering_distance_rows = "spearman", 
            gap = unit(0.5,"mm"),column_gap = unit(0.5,"mm"),
            row_title_gp = gpar(fontsize=6),
            column_names_gp = gpar(fontsize=6),
            column_title_gp = gpar(fontsize=6),
            column_dend_height = unit(5,"mm"),
            bottom_annotation = c(ha,cluster_anno),
            column_title  = switch(r,c("High","Low"),c("Low","High")),
            column_title_side = "bottom",
            heatmap_legend_param = list(title_gp = gpar(fontsize =5),
                                        labels_gp = gpar(fontsize =5),
                                        grid_width = unit(1, "mm"),
                                        legend_height=unit(1,"cm")),
            clustering_method_rows = "ward.D2",
            column_split =ncluster,show_column_names = F)
  
  p1=draw(p)
  
  
  #survival=read_tsv("E:/cancer genome/liver/????/RNA_cluster/NMF/lihc_tcga_clinical_data.tsv")
  for (n in 1:ncluster){
    assign(paste0("class",n),colnames(gene_weight_scale)[column_order(p1)[[n]]])
    assign(paste0("stage_class",n),
           survival[which(survival$`Sample ID` %in% 
                            gsub("A$","",gsub("\\.","-",get(paste0("class",n))))),]$
             `Neoplasm Disease Stage American Joint Committee on Cancer Code`)
  }
  
  for (n in 1:ncluster){
    assign(paste0("class",n),
           gsub("-01A","",gsub("\\.","-",get(paste0("class",n)))))
  }
  
  sample_class=data.frame()
  for (n in 1:ncluster){
    sample_class=rbind.data.frame(sample_class,
                                  cbind("Patient ID"=get(paste0("class",n)),
                                        class=n,
                                        stage=get(paste0("stage_class",n))))
  }
  
  survival = left_join(survival,sample_class,by="Patient ID")
  
  survival$Cluster=as.numeric(survival$class)
  
  fit <- survfit(Surv(time=survival$`Overall Survival (Months)`/12,
                      event=survival$`Overall Survival Status`) ~Cluster,
                 data = survival)
  
  p_os=ggsurvplot(fit, data = survival,pval = T,pval.method = T,
                  palette = switch(r,c("#FFA500","#ADD8E6"),
                                   c("#ADD8E6","#FFA500")),
                  size=0.4,censor.size=1,
                  legend.title=label,
                  legend.labs=switch(r,c("High","Low"),
                                     c("Low","High")) ,
                  ggtheme = theme_classic(base_size = 5,
                                          base_line_size = 0.15),
                  pval.size=1.5,title=element_blank())+
    ylab("Overall suvival probability")+xlab("Time(year)")
  
  survival$`Disease Free Status`=
    as.numeric(substr(survival$`Disease Free Status`,1,1))
  
  fit <- survfit(Surv(survival$`Disease Free (Months)`/12,
                      survival$`Disease Free Status`) ~ Cluster,  
                 data = survival)
  p_dfs=ggsurvplot(fit, data = survival,pval = T,title="",
                   ggtheme = theme_survminer() + 
                     theme(plot.title = element_text(hjust = 0.5,
                                                     face = "bold")))+
    ylab("Disease Free Suvival Probability")+xlab("Time(year)")
  
  sample_class$stage=gsub("[A,B,C]$","",sample_class$stage)
  sample_class=sample_class[which(!is.na(sample_class$stage)),]
  
  sample_class$stage=factor(
    sample_class$stage,levels=c("Stage I",
                                "Stage II",
                                "Stage III",
                                "Stage IV")
  )
  
  tmp=sample_class%>%group_by(class)%>%summarise(table(stage))
  tmp1=tmp$`table(stage)`
  
  f_p=fisher.test(matrix(tmp1,nrow = 4,byrow = F))
  f_p= format(f_p$p.value,digit=3,scientific = T)
  f_p_d=as.data.frame(matrix(c(f_p,"1","2",1.05),ncol = 4))
  colnames(f_p_d)=c("p","group1","group2","y.position")
  color=c("#FDDF71",
          "#F28F2F",
          "#DE0918",
          "#8B0000")
  names(color)=c("Stage I",
                 "Stage II",
                 "Stage III",
                 "Stage IV")
  f_p_d$y.position=as.numeric(f_p_d$y.position)
  sample_number=cbind.data.frame(y=0,table(sample_class$class))
  
  p_stage=ggplot(data=sample_class)+
    geom_bar(aes(x=class,fill=stage),position = "fill",width = 0.6) +
    theme_classic()+theme()+scale_fill_manual(values = color,na.translate=F)+
    geom_text(data=sample_number, 
              aes(x=Var1,y=y,label=paste0("n=",Freq)),vjust=2,size=1.5)+
    scale_y_continuous(limits=c(0,1.1),breaks =c(0.00,0.25,0.50,0.75,1.00),
                       expand = c(0.1,0.01) )+
    xlab(paste0(label," cluster"))+ylab("Proportion")+labs(fill="Stage")+
    scale_x_discrete(labels=switch(r,c("High","Low"),c("Low","High")))+
    stat_pvalue_manual(f_p_d,step.increase = 0.04,size=1.5,bracket.size = 0.15)+
    theme(text=element_text(size = 6),
          legend.position = "right",
          legend.spacing = unit(0,"cm"),
          legend.key.width = unit(0.1, 'cm'),
          legend.key.height = unit(0.3, 'cm'),
          legend.box.spacing = unit(0.01,"cm"),
          plot.title=element_text(hjust = 0.5),
          line = element_line(size = 0.1),
          axis.ticks.length = unit(0.03,"cm"))
  
  p_purity=pic[[n_col]]
  
  marker_all=read_tsv("E:/cancer genome/liver/analysis/decoder/lihc/HCC_primary/Final_marker_genes.txt")
  v1_marker=c(na.exclude(marker_all[,n_col]))[[1]]
  gene.df <- bitr(v1_marker, fromType = "SYMBOL",
                  toType = c("ENTREZID"),
                  OrgDb = org.Hs.eg.db)
  
  marker_go_1=enrichGO(gene.df$ENTREZID, 
                       "org.Hs.eg.db", ont = 'ALL',
                       pvalueCutoff  = 0.05,pAdjustMethod = "BH",  
                       qvalueCutoff  = 0.01, readable=T,pool = 1)
  
  result=marker_go_1@result
  result1=result[1:10,]
  result1$Description=sub(" process","",result1$Description)
  result1$Description=factor(result1$Description,levels =rev(result1$Description))
  library(scales)
  reverselog_trans <- function(base = exp(1)) {
    trans <- function(x) -log(x, base)
    inv <- function(x) base^(-x)
    trans_new(paste0("reverselog-", format(base)), trans, inv, 
              log_breaks(base = base), 
              domain = c(1e-100, Inf))
  }
  
  
  fmt_dcimals <- function(decimals=0){
    function(x) format(x,digits = decimals,scientific = T)
  }
  if (length(which(marker_go_1@result$p.adjust<0.01))>=1){
    p_go=ggplot(result1)+geom_bar(aes(x=Description,y=p.adjust),
                                  stat = "identity",fill="#FC8D62")+
      scale_y_continuous(trans = reverselog_trans(base=10),expand = c(0,0))+
      coord_flip()+
      scale_x_discrete()+
      xlab(paste0(label," factor genes' GO term"))+ylab("Adjusted pvalue")+theme_classic()+
      theme(text=element_text(size = 6),
            legend.spacing = unit(0,"cm"),
            plot.title=element_text(hjust = 0.5),
            legend.key.width = unit(0.1, 'cm'),
            legend.key.height = unit(0.3, 'cm'),
            line = element_line(size = 0.1),
            legend.box.spacing = unit(0.01,"cm"),
            axis.ticks.length = unit(0.03,"cm"))
  }else{
    p_go=NULL
  }
  kk <- enrichKEGG(gene  = gene.df$ENTREZID,
                   organism   = 'hsa')
  result=kk@result
  result1=result[1:10,]
  result1$Description=sub(" process","",result1$Description)
  result1$Description=factor(result1$Description,levels =rev(result1$Description))
  if (length(which(kk@result$p.adjust<0.01))>=1){
    p_kegg=ggplot(result1)+geom_bar(aes(x=Description,y=p.adjust),
                                    stat = "identity",fill="#FC8D62")+
      scale_y_continuous(trans = reverselog_trans(base=10),expand = c(0,0))+
      coord_flip()+
      scale_x_discrete()+
      xlab(paste0(label," factor genes' KEGG term"))+ylab("Adjusted pvalue")+theme_classic()+
      theme(text=element_text(size = 6),
            legend.spacing = unit(0,"cm"),
            plot.title=element_text(hjust = 0.5),
            legend.key.width = unit(0.1, 'cm'),
            legend.key.height = unit(0.3, 'cm'),
            line = element_line(size = 0.1),
            legend.box.spacing = unit(0.01,"cm"),
            axis.ticks.length = unit(0.03,"cm"))
  }else{
    p_kegg=NULL
  }
  p1_g=grid.grabExpr(draw(p1,padding = unit(c(2, 2, 2, 9), "mm"),))
  pl1=ggarrange(p_purity,p_go,ncol = 2,labels = c("a","b"),widths = c(0.8,1.2))
  pl2=ggarrange(p1_g,labels = "c")
  pl3=ggarrange(p_os$plot,p_stage,ncol = 2,labels = c("d","e"))
  final_p1=ggarrange(pl1,pl2,pl3,ncol = 3,widths=c(2,2,1.5))
  final_p=final_p1
  p_list=list(p_purity,p_go,p1_g,p_os$plot,p_stage,p_kegg)
  return(p_list)
}


p2=plot_fun(marker= marker[c(2)],n_col=2, data = gene_count,ncluster=2)


marker=read.table("../../4/Final_marker_genes.txt",sep = "\t",header = T)
colnames(marker)=sub("_.*","",colnames(marker))
data=read.table("E:/cancer genome/liver/analysis/decoder/lihc/HCC_primary/deseq_log2_primary.tsv",header = T,row.names = 1)

marker_top=marker[1:100,c(1,2,3,4)]



other_studies=read.table("../../4/literatre_pathway.tsv")
other_studies=other_studies[colnames(data),]
ha=HeatmapAnnotation(Lee_2004=other_studies$Lee.2004,
                     Boyault_2007=other_studies$Boyault.2007,
                     Chiang_2008=other_studies$Chiang.2008,
                     Hoshida_2009=other_studies$Hoshida.2009,
                     "CTNNB1 mutation"=CTNNB1_mut,
                     annotation_name_gp= gpar(fontsize =5),
                     annotation_legend_param = list(title_gp = gpar(fontsize =5),
                                                    labels_gp = gpar(fontsize =5),
                                                    grid_width = unit(1, "mm"),
                                                    grid_height=unit(2,"mm"),
                                                    legend_height = unit(1, "cm")),
                     simple_anno_size  = unit(2,"mm"),
                     col = list(Lee_2004=c("SURVIVAL_DN"="hotpink4","SURVIVAL_UP"="khaki3"),
                                Hoshida_2009=c("S1"=brewer.pal(8,"Set1")[3],"S2"=brewer.pal(8,"Set1")[4],"S3"=brewer.pal(8,"Set1")[5]),
                                Chiang_2008=c("CTNNB1"=brewer.pal(5,"Accent")[1],"INTERFERON"=brewer.pal(5,"Accent")[2],
                                              "POLYSOMY7"=brewer.pal(5,"Accent")[3],"PROLIFERATION"=brewer.pal(5,"Accent")[4],
                                              "UNANNOTATED"=brewer.pal(5,"Accent")[5]),
                                "CTNNB1 mutation"=c("TRUE"="brown1","FALSE"="darkgreen"),
                                Boyault_2007=c("G12"="pink","G3"="bisque","G56"="dodgerblue")))



col_fun = colorRamp2(c(-2,0,2),c("blue","white","red"))

marker_top=apply(marker_top,2,function(x) intersect(x,rownames(data)))
factor=c()
for (i in 1:ncol(marker_top)){
  factor=c(factor,rep(i,length(marker_top[,i])))
}


factor=unlist(factor)
marker_col=rainbow(length(unique(factor)))
factor=factor%>%sub("X5.5","Wnt",.)%>%sub("X5.1","metabolism",.)%>%
  sub("X5.3","ECM",.)%>%sub("X7.5","Cell Cycle",.)
names(marker_col)=unique(factor)
gene_anno = rowAnnotation(factor = factor,
                          show_annotation_name=F,
                          simple_anno_size  = unit(2,"mm"),
                          annotation_name_gp= gpar(fontsize =5),
                          show_legend=F,
                          annotation_legend_param = list( title_gp = gpar(fontsize =5),
                                                          labels_gp = gpar(fontsize =5),
                                                          grid_width = unit(1, "mm"),
                                                          grid_height=unit(2,"mm")),
                          col=list(factor=marker_col))

gene=unlist(c(marker_top))
select_data=data[gene,]
select_data=t(scale(t(select_data)))

# result_Sp_wd1=ConsensusClusterPlus(t(select_data),maxK=8,reps=500,pItem=0.8,
#                                    pFeature=0.8,plot = "pdf",seed=123,
#                                    innerLinkage="ward.D2",
#                                    distance="spearman",title = "geneConsensus")
# saveRDS(result_Sp_wd1,file = "result_Sp_wd1.RDS")
result_Sp_wd1=readRDS("../result_Sp_wd1.RDS")

# result_Sp_wd=ConsensusClusterPlus(select_data,maxK=8,reps=500,pItem=0.8,
#                                   pFeature=0.8,plot = "pdf",seed = 123,
#                                   innerLinkage="ward.D2",
#                                   distance="spearman",title = "consensusplot")
# saveRDS(result_Sp_wd,file = "result_Sp_wd.RDS")
result_Sp_wd=readRDS("../result_Sp_wd.RDS")

select_data1=select_data
colnames(select_data1)=sub("^TCGA.*$","",colnames(select_data))
column_tree=color_branches(as.dendrogram(result_Sp_wd[[4]]$consensusTree),
                           k=4,
                           col = c(RColorBrewer::brewer.pal(7,"Set2")[c(5,1,2,4)]))
split=factor(cutree(result_Sp_wd[[4]]$consensusTree,4),levels=c(2,1,4,3))

colordend=c(RColorBrewer::brewer.pal(7,"Set2")[c(5,1,2,4)])[split]
p_heatmap=Heatmap(select_data1,cluster_columns =  T,
                  clustering_distance_columns = "spearman",
                  column_split=split,col=col_fun,
                  cluster_rows =  result_Sp_wd1[[4]]$consensusTree,
                  row_split =  4,show_row_names = F,show_column_names =F,
                  row_title_gp = gpar(fontsize=6),name = "Z-score",
                  gap = unit(0.5,"mm"),
                  show_row_dend = F,
                  column_title_gp = gpar(fontsize=6),
                  row_title = c("ECM","Wnt","Metabolism","Cell Cycle"),
                  row_title_side = "right",
                  column_dend_height = unit(5,"mm"),
                  column_title = c("Wnt","Classic","Hybrid","Proliferation"),
                  column_title_side = "bottom",
                  heatmap_legend_param = list(title_gp = gpar(fontsize =6),
                                              labels_gp = gpar(fontsize =5),
                                              grid_width = unit(1, "mm"),
                                              grid_height=unit(2,"mm"),
                                              legend_height=unit(1.5,"cm")),
                  bottom_annotation =
                    c(ha,HeatmapAnnotation(
                      subtype=anno_block(
                        gp=gpar(fill=c(
                          RColorBrewer::brewer.pal(7,"Set2")[c(5,1,2,4)])),
                        height = unit(2.5,"mm")),border = F)),
                  right_annotation = gene_anno)
p_heatmap=draw(p_heatmap,column_title="", padding = unit(c(2, 2, 2, 2), "mm"), 
               annotation_legend_side = "bottom")

# pdf("test.pdf",paper = "a4",height = 4)
# p_heatmap
# dev.off()
V1=unlist(lapply(column_order(p_heatmap),function(x){colnames(select_data)[x]}))
V2=c()
for (i in 1:4){V2=c(V2,rep(i,length(column_order(p_heatmap)[[i]])))}
class=cbind.data.frame(V2,V1)

class_tcga=class[which(substr(class$V1,1,4)=="TCGA"),]
#write.table(class_tcga,file = "class_tcga.tsv")
class_tcga$V1=gsub("\\.","-",class_tcga$V1)
rownames(class_tcga)=class_tcga$V1

survival=read_tsv("E:/cancer genome/liver/analysis//RNA_cluster/NMF/lihc_tcga_clinical_data.tsv")
survival$`Disease Free Status`=as.numeric(substr(survival$`Disease Free Status`,1,1))
survival$`Overall Survival Status`=as.numeric(substr(survival$`Overall Survival Status`,1,1))
purity=read.table("E:/cancer genome/liver/analysis/RNA_cluster/gene_symbol/TCGA_data.tsv",header = T)
purity=purity%>%dplyr::select(tumor_sample_barcode,purity)%>%unique()
#purity$tumor_sample_barcode=gsub("-",".", purity$tumor_sample_barcode)
rownames(purity)=purity$tumor_sample_barcode
sample_select=intersect(rownames(purity),class_tcga$V1)

survival$`Sample ID`=paste0(survival$`Sample ID`,"A")
survival_select=survival[survival$`Sample ID`%in%class_tcga$V1,]
survival_select$Cluster=class_tcga[survival_select$`Sample ID`,1]
library(survival)
library(survminer)
fit <- survfit(Surv(time=survival_select$`Overall Survival (Months)`/12,
                    event=survival_select$`Overall Survival Status`) ~Cluster,
               data = survival_select)

p_os=ggsurvplot(fit, data = survival_select,pval = T,pval.method = T,
                #log.rank.weights="n",
                palette = c(RColorBrewer::brewer.pal(7,"Set2")[c(5,1,2,4)]),
                title=element_blank(),pval.size=1.5,
                legend.title="Subtype",legend.labs=c("Wnt","Classic","Hybrid","Proliferation"),
                size=0.4,censor.size=0.4,
                ggtheme = theme_classic(base_size = 5,base_line_size = 0.15) + 
                  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
                        plot.margin = unit(c(2,2,2,2),"mm")))+
  ylab("Overall suvival probability")+xlab("Time(year)")

survival_select$`Disease Free Status`=as.numeric(substr(survival_select$`Disease Free Status`,1,1))
fit <- survfit(Surv(survival_select$`Disease Free (Months)`/12,
                    survival_select$`Disease Free Status`) ~ Cluster,  
               data = survival_select)
p_dfs=ggsurvplot(fit, data = survival_select,pval = T,
                 title="",pval.method = T,log.rank.weights="n",
                 ggtheme = theme_survminer() + 
                   theme(plot.title = element_text(hjust = 0.5, face = "bold"),
                         plot.margin = unit(c(2,20,2,2),"mm")))+
  ylab("Disease free suvival probability")+xlab("Time(year)")

sample_class=survival_select%>%
  dplyr::select(`Neoplasm Disease Stage American Joint Committee on Cancer Code`,`Sample ID`,Cluster)%>%
  dplyr::rename("stage"=`Neoplasm Disease Stage American Joint Committee on Cancer Code`,"class"=Cluster)
sample_class$stage=gsub("[A,B,C]$","",sample_class$stage)
sample_class=sample_class[which(!is.na(sample_class$stage)),]

sample_class$stage=factor(
  sample_class$stage,levels=c("Stage I",
                              "Stage II",
                              "Stage III",
                              "Stage IV")
)

s=chisq.test(sample_class$stage,sample_class$class)

tmp=sample_class%>%group_by(class)%>%summarise(table(stage))
tmp1=tmp$`table(stage)`

f_p=fisher.test(matrix(tmp1,nrow = 4,byrow = F)[,c(1,2)])
f_p= format(f_p$p.value,digit=2,scientific = T)
f_p_d=as.data.frame(matrix(c(f_p,"1","2",1.05),ncol = 4))
colnames(f_p_d)=c("p","group1","group2","y.position")

f_p=fisher.test(matrix(tmp1,nrow = 4,byrow = F)[,c(2,3)])
f_p= format(f_p$p.value,digit=2,scientific = T)
f_p_d=rbind(f_p_d,c(f_p,"2","3",1.05))

f_p=fisher.test(matrix(tmp1,nrow = 4,byrow = F)[,c(3,4)])
f_p= format(f_p$p.value,digit=2,scientific = T)
f_p_d=rbind(f_p_d,c(f_p,"3","4",1.05))

f_p=fisher.test(matrix(tmp1,nrow = 4,byrow = F)[,c(1,3)])
f_p= format(f_p$p.value,digit=2,scientific = T)
f_p_d=rbind(f_p_d,c(f_p,"1","3",1.15))

f_p=fisher.test(matrix(tmp1,nrow = 4,byrow = F)[,c(2,4)])
f_p= format(f_p$p.value,digit=2,scientific = T)
f_p_d=rbind(f_p_d,c(f_p,"2","4",1.15))

f_p=fisher.test(matrix(tmp1,nrow = 4,byrow = F)[,c(1,4)])
f_p= format(f_p$p.value,digit=2,scientific = T)
f_p_d=rbind(f_p_d,c(f_p,"1","4",1.25))
f_p_d$y.position=as.numeric(f_p_d$y.position)

color=c("#FDDF71",
        "#F28F2F",
        "#DE0918",
        "#8B0000")
names(color)=c("Stage I",
               "Stage II",
               "Stage III",
               "Stage IV")
f_p_d$y.position=as.numeric(f_p_d$y.position)
sample_number=cbind.data.frame(y=0,table(sample_class$class))
sample_class$class=as.character(sample_class$class)
f_p_d=f_p_d[c(1:4,6),]
p_stage=ggplot(data=sample_class)+
  geom_bar(aes(x=class,fill=stage),position = "fill",width = 0.6) +
  theme_classic()+theme()+
  scale_fill_manual(values = color,na.translate=F)+
  scale_x_discrete(labels=c("Wnt","Classic","Hybrid","Proliferation"))+
  geom_text(data=sample_number, 
            aes(x=Var1,y=y,label=paste0("n=",Freq)),vjust=1.8,size=1.5)+
  scale_y_continuous(limits=c(0,1.15),expand = c(0.1,0.01),
                     breaks =c(0.00,0.25,0.50,0.75,1.00))+
  xlab("Subtype")+ylab("Proportion")+labs(fill="Stage")+
  geom_text(data=data.frame(),
            aes(x=2.5,y=1.1,
                label=paste0("p=",format(s$p.value,digit=2,scientific = T))),
            size=1.5,color="#4D4D80")+
  #stat_pvalue_manual(f_p_d,step.increase = 0.04,size = 1.5,bracket.size = 0.15)+
  theme(text=element_text(size = 5),
        legend.spacing = unit(0,"cm"),
        plot.title=element_text(hjust = 0.5),
        legend.key.width = unit(0.1, 'cm'),
        legend.key.height = unit(0.3, 'cm'),
        line = element_line(size = 0.1),
        legend.box.spacing = unit(0.01,"cm"),
        axis.ticks.length = unit(0.03,"cm"),
        axis.text.x = element_text(angle = 45,vjust=1,hjust=1))

f1=ggarrange(p2[[2]],p2[[3]],ncol=2,widths = c(1,1.5),
             labels = c("a","b"),
             font.label = list(size=8))
f2=ggarrange(p2[[5]],p2[[4]],p_stage,p_os$plot,
             ncol=4,widths = c(1,1,1,1.2),labels = c("c","d","f","g"),
             font.label = list(size=8))
f=ggarrange(f1,f2,
            grid.grabExpr(draw(p_heatmap)),
            nrow = 3,heights = c(1.5,1.2,2.5),labels = c("","","e"),
            font.label = list(size=8))


purity=read.table("E:/cancer genome/liver/analysis/RNA_cluster/gene_symbol/TCGA_data.tsv",header = T)
samples=intersect(class_tcga$V1,purity$tumor_sample_barcode)
purity=purity[purity$tumor_sample_barcode%in%samples,]
purity$class=class_tcga[purity$tumor_sample_barcode,1]
data=purity
data$class=as.character(data$class)

mycompare=list(c(1,2),c(2,3),c(3,4),c(1,4))
p_purity=ggplot(data = data,aes(class,purity))+
  geom_violin(aes(fill=class),trim=FALSE,size=0.15)+
  scale_x_discrete(labels=c("Wnt","Classic","Hybrid","Proliferation"))+
  scale_fill_manual(values=c(RColorBrewer::brewer.pal(7,"Set2")[c(5,1,2,4)]))+
  geom_boxplot(width=0.2,outlier.size = 0.5,size=0.15)+
  stat_compare_means(comparisons = mycompare,size=1.5,bracket.size = 0.15)+
  theme_classic()+xlab("Subtype")+ylab("Purity")+
  theme(legend.position = "none",text = element_text(size=5),
        plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size = 5),
        strip.text = element_text(size=5),
        strip.background = element_rect(size=0.1),
        line = element_line(size = 0.1),
        legend.box.spacing = unit(0.01,"cm"),
        axis.ticks.length = unit(0.03,"cm"),
        axis.text.x = element_text(angle = 45,vjust = 0.5))

immune_score=read.table("../liver_hepatocellular_carcinoma_RNAseqV2.txt",header = T)
immune_score$ID=paste0(immune_score$ID,"A")
immune_samples=intersect(class_tcga$V1,immune_score$ID)
immune=immune_score[match(immune_samples,immune_score$ID),"Immune_score"]
stromal=immune_score[match(immune_samples,immune_score$ID),"Stromal_score"]
purity_e=immune_score[match(immune_samples,immune_score$ID),"ESTIMATE_score"]
class=class_tcga[immune_samples,1]
data=cbind.data.frame(class,immune,stromal,purity_e)
data$class=as.character(data$class)

mycompare=list(c(1,2),c(2,3),c(1,3),c(1,4))
p_immune=ggplot(data = data,aes(class,immune))+
  geom_violin(aes(fill=class),trim=FALSE,size=0.15)+
  scale_x_discrete(labels=c("Wnt","Classic","Hybrid","Proliferation"))+
  scale_fill_manual(values=c(RColorBrewer::brewer.pal(7,"Set2")[c(5,1,2,4)]))+
  geom_boxplot(width=0.2,outlier.size = 0.5,size=0.15)+
  stat_compare_means(comparisons = mycompare,size=1.5,bracket.size = 0.15,method = "t.test")+
  theme_classic()+xlab("Subtype")+ylab("Immune score")+
  theme(legend.position = "none",text = element_text(size=5),
        plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size = 5),
        strip.text = element_text(size=5),
        strip.background = element_rect(size=0.1),
        line = element_line(size = 0.1),
        legend.box.spacing = unit(0.01,"cm"),
        axis.ticks.length = unit(0.03,"cm"),
        axis.text.x = element_text(angle = 45,vjust = 1,hjust=1))
ggsave("immune_score.pdf",width = 39.511,height = 46.704,units = "mm")


file.dir="E:/cancer genome/liver/TCGA/CNV/"
files=list.files(path = "E:/cancer genome/liver/TCGA/CNV/",
                 pattern="seg.v2.txt",recursive = T)

info=read.table("E:/cancer genome/liver/TCGA/CNV/gdc_sample_sheet.2021-07-02.tsv",header = T,sep="\t")
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

cytoband_ifile <- "E:/cancer genome/liver/TCGA/CNV/cytoBand.txt"

cytoband <- read_tsv(cytoband_ifile, col_names = F) %>%
  dplyr::filter(X1 %in% paste0("chr",c(1:22,"X")))

colnames(cytoband) = c("chromosome", "Start", "End", "arm", "Misc")
cytoband <- cytoband %>% mutate(chromosome = gsub("chr(.*)", "\\1", chromosome))

cytoband_length <- cytoband %>% group_by(chromosome) %>% summarise(length = max(End))
cytoband_pos <- cytoband_length
cytoband_pos <- cytoband_pos %>% arrange(factor(cytoband_pos$chromosome, levels = seq(1,22)))
# Get absolute position of chrom end
cytoband_pos <- cytoband_pos %>% mutate(pos=cumsum(length))

# Get arm position
cytoband_arm <- cytoband %>% dplyr::filter(Misc=="acen")
cytoband_arm <- set_absolute_pos(cytoband_arm, cytoband_length)

# Set label for p and q
cytoband_arm <- cytoband_arm %>% group_by(chromosome) %>% summarise(Start=min(Start), End=max(End))
cytoband_arm <- cytoband_arm %>% mutate(start_label = paste0(chromosome, "p"),
                                        end_label = paste0(chromosome, "q"))

# Set label for mid point of p and q
cytoband_arm <- cytoband_arm %>% mutate(midpoint=(Start + End)/2)
cytoband_arm <- cytoband_arm %>% mutate(size=End-Start)

select=info[which(info$Data.Type=="Masked Copy Number Segment"),]
select=select%>%filter(Sample.Type=="Primary Tumor")


files=paste(select$File.ID,select$File.Name,sep = "/")


all_seg=data.frame()
for (i in 1:length(files)){
  seg=read.table(paste0(file.dir,files[i]),sep = "\t",header = T)
  seg$GDC_Aliquot=select[i,7]
  all_seg=rbind(all_seg,seg)
}
all_seg=all_seg%>%dplyr::rename("chromosome"="Chromosome")%>%
  mutate(chromosome=sub("chr","",chromosome))%>%
  filter(chromosome%in%c(1:22,"X"))

all_segs=all_seg%>%dplyr::rename("Tumor_Sample_Barcode"="GDC_Aliquot")
all_segs$adjustedCN=all_segs$Segment_Mean

all_segs=all_segs[which(all_segs$Tumor_Sample_Barcode%in%class_tcga$V1),]
all_segs$Tumor_Sample_Barcode=factor(all_segs$Tumor_Sample_Barcode,levels = class_tcga$V1)

all_seg=all_segs
all_seg$adjustedCN=NULL
all_segs1=all_seg[which(all_seg$Tumor_Sample_Barcode %in% class_tcga$V1[which(class_tcga$V2==1)]),]
all_segs2=all_seg[which(all_seg$Tumor_Sample_Barcode %in% class_tcga$V1[which(class_tcga$V2==2)]),]
all_segs3=all_seg[which(all_seg$Tumor_Sample_Barcode %in% class_tcga$V1[which(class_tcga$V2==3)]),]
all_segs4=all_seg[which(all_seg$Tumor_Sample_Barcode %in% class_tcga$V1[which(class_tcga$V2==4)]),]
all_segs <- set_absolute_pos(all_segs, cytoband_length)

rownames(class_tcga)=class_tcga$V1
all_segs$class=class_tcga[all_segs$Tumor_Sample_Barcode,2]




# all_segs1=read.table("tcga_seg_class1.tsv",header = T)
# all_segs2=read.table("tcga_seg_class2.tsv",header = T)
# all_segs3=read.table("tcga_seg_class3.tsv",header = T)
# all_segs4=read.table("tcga_seg_class4.tsv",header = T)
cyto=read_tsv("E:/cancer genome/liver/TCGA/CNV/cytoBand.txt",col_names = F)
cyto=cyto%>%filter(X1%in%c(paste0("chr",c(1:22,"X"))))
genome_len=sum(cyto$X3-cyto$X2)
segs=all_segs1
samples=unique(segs$Tumor_Sample_Barcode)
GII1=c();breakpoint1=c();cnv1=c()
for (s in samples){
  tmp=segs%>%filter(Tumor_Sample_Barcode==s,Segment_Mean>0.2|Segment_Mean< -0.2)
  len=sum(tmp$End-tmp$Start)
  breakpoint1=c(breakpoint1,nrow(tmp))
  cnv1=c(cnv1,mean(abs(tmp$Segment_Mean)))
  GII1=c(GII1,len/genome_len)
}
mean(GII1)

segs=all_segs2
samples=unique(segs$Tumor_Sample_Barcode)
GII2=c();breakpoint2=c();cnv2=c()
for (s in samples){
  tmp=segs%>%filter(Tumor_Sample_Barcode==s,Segment_Mean>0.2|Segment_Mean< -0.2)
  len=sum(tmp$End-tmp$Start)
  breakpoint2=c(breakpoint2,nrow(tmp))
  cnv2=c(cnv2,mean(abs(tmp$Segment_Mean)))
  GII2=c(GII2,len/genome_len)
}
mean(GII2)

segs=all_segs3
samples=unique(segs$Tumor_Sample_Barcode)
GII3=c();breakpoint3=c();cnv3=c()
for (s in samples){
  tmp=segs%>%filter(Tumor_Sample_Barcode==s,Segment_Mean>0.2|Segment_Mean< -0.2)
  len=sum(tmp$End-tmp$Start)
  breakpoint3=c(breakpoint3,nrow(tmp))
  cnv3=c(cnv3,mean(abs(tmp$Segment_Mean)))
  if(len>genome_len){
    GII3=c(GII3,len/(2*genome_len))
  }else{
    GII3=c(GII3,len/genome_len)
  }
  
}
mean(GII3)

segs=all_segs4
samples=unique(segs$Tumor_Sample_Barcode)
GII4=c();breakpoint4=c();cnv4=c()
for (s in samples){
  tmp=segs%>%filter(Tumor_Sample_Barcode==s,Segment_Mean>0.2|Segment_Mean< -0.2)
  len=sum(tmp$End-tmp$Start)
  breakpoint4=c(breakpoint4,nrow(tmp))
  if(len>genome_len){
    GII4=c(GII4,len/(2*genome_len))
  }else{
    GII4=c(GII4,len/genome_len)
  }
  cnv4=c(cnv4,mean(abs(tmp$Segment_Mean)))
}
mean(GII4)
GII1=cbind.data.frame(GII=GII1,class="GII1")
GII2=cbind.data.frame(GII=GII2,class="GII2")
GII3=cbind.data.frame(GII=GII3,class="GII3")
GII4=cbind.data.frame(GII=GII4,class="GII4")
GII=rbind(GII1,GII2,GII3,GII4)
library(ggpubr)

my_comparisons=list(c(1,2),c(2,3),c(2,4))
p_GII=ggplot(data = GII,aes(class,GII))+
  geom_violin(aes(fill=class),trim=FALSE,size=0.15)+
  geom_boxplot(width=0.2,outlier.size = 0.5,size=0.15)+xlab("Subtype")+
  scale_x_discrete(labels=c("Wnt","Classic","Hybrid","Proliferation"))+
  scale_fill_manual(values=c(RColorBrewer::brewer.pal(7,"Set2")[c(5,1,2,4)]))+
  theme_classic()+stat_compare_means(comparisons = my_comparisons,size=1.5,bracket.size = 0.15)+
  theme(legend.position = "none",text = element_text(size=5),
        plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size = 5),
        strip.text = element_text(size=5),
        strip.background = element_rect(size=0.1),
        line = element_line(size = 0.1),
        legend.box.spacing = unit(0.01,"cm"),
        axis.ticks.length = unit(0.03,"cm"),
        axis.text.x = element_text(angle = 45,vjust = 1,hjust=1))

breakpoint1=cbind.data.frame(breakpoint=breakpoint1,class="breakpoint1")
breakpoint2=cbind.data.frame(breakpoint=breakpoint2,class="breakpoint2")
breakpoint3=cbind.data.frame(breakpoint=breakpoint3,class="breakpoint3")
breakpoint4=cbind.data.frame(breakpoint=breakpoint4,class="breakpoint4")
breakpoint=rbind(breakpoint1,breakpoint2,breakpoint3,breakpoint4)
p_breakpoint=ggplot(data = breakpoint,aes(class,breakpoint))+
  geom_violin(aes(fill=class),trim=FALSE,size=0.15)+
  geom_boxplot(width=0.2,outlier.size = 0.5,size=0.15)+
  scale_x_discrete(labels=c(1,2,3,4))+
  scale_fill_manual(values=c(RColorBrewer::brewer.pal(7,"Set2")[c(5,1,2,4)]))+
  theme_classic()+
  stat_compare_means(comparisons = my_comparisons,size=1.5,bracket.size = 0.15)+
  theme(legend.position = "none",text = element_text(size=5),
        plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size = 5),
        line = element_line(size = 0.1),
        legend.box.spacing = unit(0.01,"cm"),
        axis.ticks.length = unit(0.03,"cm"))

ggarrange(p_purity,p_GII,p_breakpoint,ncol = 3,nrow = 1)


immune=read.csv("E://cancer genome/liver/analysis/cibersort/CIBERSORTx_Job4_Results.csv",header = T)
immune=immune[,-c(24:26)]
immune_matrix=t(as.matrix(immune[,-1]))
colnames(immune_matrix)=immune$Mixture

data=tidyr::gather(immune,key = celltype,value = value,-Mixture)


class_tcga$V1=gsub("-",".",class_tcga$V1)
#data$Mixture=factor(data$Mixture,levels = class$V1)

class1=class_tcga$V1[which(class_tcga$V2==1)]
class2=class_tcga$V1[which(class_tcga$V2==2)]
class3=class_tcga$V1[which(class_tcga$V2==3)]
class4=class_tcga$V1[which(class_tcga$V2==4)]
immune_matrix1=immune_matrix[,class1]
immune_matrix2=immune_matrix[,class2]
immune_matrix3=immune_matrix[,class3]
immune_matrix4=immune_matrix[,class4]
tmp=cbind(immune_matrix1,immune_matrix3)
tmp=t(scale(t(tmp)))
Heatmap(tmp,cluster_columns = F)
p=c()
for( r in 1:22 ){
  p=c(p,t.test(immune_matrix2[r,],immune_matrix4[r,])$p.value)
}


diff=rownames(immune_matrix)[which(p.adjust(p)<0.05)]

select_data=data[(data$celltype%in%diff)&(data$Mixture%in%c(class2,class4)),]
select_data$class=as.character(ifelse(select_data$Mixture %in% class2,2,4))

library(ggpubr)

select_data=data[data$celltype %in% diff,]
class_name=class_tcga$V2
names(class_name)=class_tcga$V1
select_data$class=as.character(class_name[select_data$Mixture])
select_data= select_data[select_data$class %in% c("2","4"),]
select_data$celltype[which(select_data$celltype=="T.cells.regulatory..Tregs.")]="Treg"

select_data$celltype = sub("^.*\\.","",select_data$celltype)
select_data$celltype[which(select_data$celltype=="resting")]="CD4+ resting T"
select_data$celltype=factor(select_data$celltype,levels=c("CD4+ resting T","Monocytes","M1","M0","Treg"))
immune_plot=ggplot(data = select_data,aes(class,value))+
  scale_x_discrete(labels=c("Classic","Proliferation"))+
  geom_violin(aes(fill=class),trim = FALSE,size=0.15)+
  facet_wrap(~ celltype,scales = "free",ncol = 5)+xlab("Subtype")+ylab("Proportion")+
  geom_boxplot(width=0.2,outlier.size = 0.5,size=0.15)+theme_classic()+
  scale_fill_manual(values=c(RColorBrewer::brewer.pal(7,"Set2")[c(1,4)]))+
  stat_compare_means(comparisons = list(c("2","4")),
                     size=1.5,bracket.size = 0.15)+
  theme(legend.position = "none",text = element_text(size=5),
        plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size = 5),
        strip.text = element_text(size=5),
        strip.background = element_blank(),
        line = element_line(size = 0.1),
        legend.box.spacing = unit(0.01,"cm"),
        axis.ticks.length = unit(0.03,"cm"),
        axis.text.x = element_text(angle = 45,vjust=1,hjust=1))
#ggsave("immune_3subtype.pdf",width = 4,height = 2)
p_g=ggarrange(p_GII,p_immune,immune_plot,ncol = 3,widths = c(1,1,2.5),labels = c("h","i","j"),
              font.label = list(size=8))
ff=ggarrange(f,p_g,nrow = 2,heights = c(5,1.2))
pdf("p3_v14.pdf",width = 6.8,height = 9.5)
ff
dev.off()





pacman::p_load(data.table,tidyverse,fgsea)
setDTthreads(10)

## this script returns multiple gene signatures including Sia et al, inflammaotry, IFNG, exhausted t-cell by calculating GSVA scores using gene lists
expression <-read.table("E:/cancer genome/liver/????/RNA_cluster/combine_tcga_DT/deseq_log2_tcga_DT_tumor.tsv")

Sia_immune_class <- fread("geneSets/Sia_et_al_immune_class_signature.txt") %>% 
  as.list()
names(Sia_immune_class) <- "ImmuneClass"
immune_signatures <- gmtPathways("geneSets/Immune_signatures.txt")

MVI_6_gene <- list(MVI_6_gene=c("ROS1", "UGT2B7", "FAS", "ANGPTL7", "GMNN", "MKI67"))

GS_full <- c(Sia_immune_class,immune_signatures,MVI_6_gene)
# PD1_18gene=list(PD1_18gene=c("CXCR6","TIGIT","CD27","PDCD1LG2","CD274",
#                              "CD8A","LAG3","NKG7","CCL5","CMKLR1","PSMB10","CXCL9","IDO1","HLA-DQA1",
#                              "CD276","STAT1","HLA-DRB1","HLA-E"))
# GS_full=c(GS_full,PD1_18gene)

#es.max<- GSVA::gsva(as.matrix(expression), GS_full, mx.diff=TRUE, verbose=FALSE, parallel.sz=1)

# write_tsv(es.max %>% t() %>% as.data.frame() %>% 
#             tibble::rownames_to_column("sample") ,"outputs/DT_GeneSignatures.tsv")



library(msigdbr)
library(GSVA)
#pathways = msigdbr(species = "Homo sapiens",category = c("H"))
pathways = msigdbr(species = "Homo sapiens")
#pathways_select= pathways[grep("METASTASIS",pathways$gs_name),]
#tmp=pathways[grep("METASTASIS",pathways$gs_name),]
#tmp=pathways[grep("AIZARANI_LIVER",pathways$gs_name),]


pathway=c("HO_LIVER_CANCER_VASCULAR_INVASION",
          "MINGUEZ_LIVER_CANCER_VASCULAR_INVASION_DN",
          "MINGUEZ_LIVER_CANCER_VASCULAR_INVASION_UP",
          'ROESSLER_LIVER_CANCER_METASTASIS_DN',
          "ROESSLER_LIVER_CANCER_METASTASIS_UP",
          "HALLMARK_E2F_TARGETS",
          "HALLMARK_G2M_CHECKPOINT",
          "BUDHU_LIVER_CANCER_METASTASIS_DN",
          "BUDHU_LIVER_CANCER_METASTASIS_UP",
          "ROESSLER_LIVER_CANCER_METASTASIS_DN",
          "ROESSLER_LIVER_CANCER_METASTASIS_UP",
          "HALLMARK_ANGIOGENESIS")
tmp=pathways[pathways$gs_name %in% pathway,]

pathways_select = tmp %>% split(x = .$gene_symbol, f = .$gs_name)
pathways_select=c(pathways_select,GS_full)

library(ComplexHeatmap)
gsva_res <- gsva(as.matrix(expression), pathways_select, parallel.sz=1,
                 min.sz=1, max.sz=500, verbose=TRUE)
#rownames(gsva_res)=sub(" ","",gsub("_"," ",sub("[^_]*","",rownames(gsva_res))))

ht=Heatmap(gsva_res,column_split = 3,row_split = 2)

# pdf("select_GSEA.pdf",width = 40,height = 10)
# h2<-draw(ht, padding = unit(c(2, 2, 2, 20), "mm"))
# dev.off()

library(tidyverse)
library(ggsci)
library(ggpubr)
info=read.table("../../../distance/all.csv",header = T,sep = ",")
subtype=read.table("../../3/class_name_tcga_DT_tumor.tsv")
subtype_vec=subtype$V2
subtype_vec=subtype_vec%>%gsub("1","Wnt",.)%>%gsub("2","Classic",.)%>%
  gsub("3","Hybrid",.)%>%gsub("4","Proliferation",.)
names(subtype_vec)=subtype$V1


GEP=read.table("outputs/combine_GEP_score.tsv",header = T)
impres=read.table("outputs/combine_IMPRES_scores.txt",header = T)
tide=read.csv("../TIDE_score_combine.csv")
tide=tide%>%arrange(Patient)%>%filter(Patient%in%colnames(gsva_res))
gsva_res=rbind(gsva_res,"GEP score"=GEP$GEP,IMPRES=impres$IMPRES,TIDE=tide$TIDE)
melt=gsva_res%>%as.data.frame(.)%>%gather()
melt=cbind(melt,pathway=rownames(gsva_res))
melt$subtype=factor(subtype_vec[melt$key],
                    levels = c("Wnt","Classic","Hybrid","Proliferation"))
###############compare subtype
two.way <- aov(value ~ subtype + pathway, data = melt)

library(ggpubr)
p1=ggplot(melt)+geom_violin(aes(y=value,x=subtype,fill=subtype))+
  geom_boxplot(aes(x=subtype,y=value),width=0.2,outlier.size = 0.5,size=0.15)+
  stat_compare_means(aes(y=value,x=subtype))+
  scale_fill_manual(values=c(RColorBrewer::brewer.pal(7,"Set2")[c(5,1,2,4)]))+
  facet_wrap(~pathway,scales = "free",nrow = 6,ncol = 4)+theme_bw()
ggsave("compare_subtype1.pdf",p1,width = 20,height = 20)
fill=rownames(gsva_res)
names(fill)=c("Angiogenesis","E2F","G2M","HO Vascular invasion",
              "Vascular invasion down","Vascular invasion",
              "Metastasis down","Metastasis","ImmuneClass","Exhaused T cell",
              "Interferon gamma","Antigen presenting cells","Inflammatory",
              "Gajewski inflammatory","Cytolytic score","MVI 6 gene","GEP score","IMPRES",
              "TIDE")
pathway_s=c("MINGUEZ_LIVER_CANCER_VASCULAR_INVASION_UP",
            "ROESSLER_LIVER_CANCER_METASTASIS_UP",
            "ImmuneClass","Cytolytic_score",
            "Inflammatory","GEP score","IMPRES","TIDE")
pathway_s%in%melt$pathway
# pathway_s=c("MINGUEZ_LIVER_CANCER_VASCULAR_INVASION_UP",
#             "ROESSLER_LIVER_CANCER_METASTASIS_UP",
#             "HALLMARK_ANGIOGENESIS",
#             "PD1_18gene","Inflammatory")
melt_DT_target_pathway=melt[which(melt$pathway %in% pathway_s & substr(melt$key,1,2)=="DT"),]

melt_DT_target_pathway$pathway=plyr::revalue(melt_DT_target_pathway$pathway,
                                             c(ImmuneClass="Immune class",
                                               MINGUEZ_LIVER_CANCER_VASCULAR_INVASION_UP="Vascular invasion",
                                               ROESSLER_LIVER_CANCER_METASTASIS_UP="Metastasis",
                                               Cytolytic_score="Cytolytic score"))

pcolor=c(RColorBrewer::brewer.pal(12,"Paired"),"darkolivegreen1")
patient_trans = unique(info$Patient_paper)
names(patient_trans)=unique(info$Patient)
pdf("select_pathway1.pdf",width = 3,height = 6)
melt_DT_target_pathway$patient_paper=patient_trans[substr(melt_DT_target_pathway$key,1,4)]
select_pathway_p=ggplot(melt_DT_target_pathway)+
  geom_violin(aes(x=patient_paper,y=value,fill=patient_paper),
              size=0.2,scale="width",adjust=1.2)+
  scale_y_continuous(limits = function(x) c(x[1],x[2]+(x[2]-x[1])/5))+
  geom_boxplot(aes(x=patient_paper,y=value),width=0.2,outlier.size = 0.2,size=0.15)+
  stat_compare_means(aes(y=value,x=patient_paper),method = "anova",size=1.5)+
  facet_wrap(pathway~.,ncol = 1,scales = "free_y",strip.position = "left")+
  scale_fill_manual(values=pcolor)+
  labs(GEP="gep",x="Patient",y="Score")+
  theme_bw()+theme(legend.position = "none",text = element_text(size=5),
                   plot.title = element_text(hjust = 0.5,size = 6),
                   axis.title = element_text(size = 5),
                   strip.text = element_text(size=5),
                   strip.background = element_rect(size=0.1),
                   line = element_line(size = 0.1),
                   legend.box.spacing = unit(0.01,"cm"),
                   axis.ticks.length = unit(0.03,"cm"),
                   axis.text.x = element_text(angle = 45,vjust=0.5))
print(select_pathway_p)
dev.off()

melt_DT=melt[which(substr(melt$key,1,2)=="DT"),]
melt_DT$patient_paper=patient_trans[substr(melt_DT$key,1,4)]
p3=ggplot(melt_DT)+geom_violin(aes(y=value,x=patient_paper))+
  geom_boxplot(aes(x=patient_paper,y=value),width=0.2,outlier.size = 0.5,size=0.15)+
  stat_compare_means(aes(y=value,x=patient_paper),method = "anova")+
  facet_wrap(pathway~.,ncol = 1)+theme_bw()
ggsave("patient_compare1.pdf",width = 25,height = 30)
compare_means(value~patient_paper,data = melt_DT,method = "anova")
res=aov(value~patient_paper,data = melt_DT)
summary(res)
###############compare interference


info=read.table("../../../distance/all.csv",header = T,sep = ",")

interference=c()
for (p in sort(unique(info$Patient_paper))){
  data=info[which(info$Patient_paper==p),]
  o_x=mean(data$X)
  o_y=mean(data$Y)
  dis=sqrt((data$X-o_x)**2+(data$Y-o_y)**2)
  data$dis=dis
  interference=c(interference,data[dis>max(dis)/2,]$Sample)
}
#melt_DT$out="in" 
#melt_DT$out[melt_DT$key %in% interference]="out"

melt_DT_target_pathway$out="Inside"
melt_DT_target_pathway$out[melt_DT_target_pathway$key %in% interference]="Outside"

out_count= melt_DT_target_pathway%>%
  group_by(pathway)%>%summarise(count=table(out))%>%
  mutate(out=names(count))

ggplot(melt_DT_target_pathway)+geom_violin(aes(y=value,x=out,fill=out))+
  xlab("Position")+scale_y_continuous(expand =expansion(mult = c(0.1,0.1)))+
  ylab("Value")+
  geom_boxplot(aes(x=out,y=value),
               width=0.2,outlier.size = 0.2,size=0.15)+
  stat_compare_means(aes(y=value,x=out),method = "anova",
                     size=2,vjust = -0.25,label.x.npc =  0.2)+
  facet_wrap(~pathway,scales = "free",nrow = 2,ncol = 4)+theme_classic()+
  theme(text = element_text(size = 7),legend.position = "none")
ggsave("patient_out2.pdf",width = 8,height = 4)
ggsave("Supplementary Fig9.pdf",width = 8,height = 4)
plot=NULL;plot_sup=NULL
for (p in sort(unique(info$Patient_paper))){
  for (f in rownames(gsva_res)){
    data=info[which(info$Patient_paper==p),]
    gsva_sel=gsva_res[,which(colnames(gsva_res) %in% data$Sample)]
    data=cbind(data,fill=gsva_sel[f,])
    plot[[f]]=c(plot[[f]],
                list(p=ggplot(data = data)+geom_point(aes(x=X,y=Y,fill=fill),
                                                      size=2,shape=21,stroke=0.4)+
                       scale_x_continuous(limits=function(x) c(x[1]-(x[2]-x[1])/10,x[2]+(x[2]-x[1])/10))+
                       scale_y_continuous(limits = function(x) c(x[1]-0.5,x[2]+0.5))+
                       scale_fill_gsea()+
                       ggtitle(paste0(p," ",names(fill)[fill==f]))+
                       theme_classic()+
                       labs(fill="Score",x="X (cm)",y="Y (cm)")+
                       theme(text=element_text(size = 7),
                             legend.spacing = unit(0,"cm"),
                             plot.title=element_text(hjust = 0.5),
                             legend.key.width = unit(0.1, 'cm'),
                             legend.key.height = unit(0.3, 'cm'),
                             line = element_line(size = 0.1),
                             legend.box.spacing = unit(0.01,"cm"),
                             axis.ticks.length = unit(0.03,"cm"),
                             strip.background = element_rect(size=0.1),
                             panel.background = element_rect(fill = "#5E4FA2"),
                             axis.text.x = element_text())
                ))
    plot_sup[[f]]=c(plot_sup[[f]],
                list(p=ggplot(data = data)+geom_point(aes(x=X,y=Y,fill=fill),
                                                      size=2,shape=21,stroke=0.4)+
                       scale_x_continuous(limits=function(x) c(x[1]-(x[2]-x[1])/10,x[2]+(x[2]-x[1])/10))+
                       scale_y_continuous(limits = function(x) c(x[1]-0.5,x[2]+0.5))+
                       scale_fill_gsea()+
                       ggtitle(paste0(p," ",names(fill)[fill==f]))+
                       theme_classic()+
                       labs(fill="Score",x="X (cm)",y="Y (cm)")+
                       theme(text=element_text(size = 7),
                             legend.spacing = unit(0,"cm"),
                             plot.title=element_text(hjust = 0.5),
                             legend.key.width = unit(0.1, 'cm'),
                             legend.key.height = unit(0.3, 'cm'),
                             line = element_line(size = 0.1),
                             legend.box.spacing = unit(0.01,"cm"),
                             axis.ticks.length = unit(0.03,"cm"),
                             strip.background = element_rect(size=0.1),
                             panel.background = element_rect(fill = "#81B5E1"),
                             axis.text.x = element_text())
                ))
    
  }
}

patient_pathway=NULL;patient_pathway1=NULL
for (i in c(1,2,3,4,7,8,9,12)){
  for (p in pathway_s){
    patient_pathway=c(patient_pathway,plot_sup[[p]][i])
  }
}
for (i in c(5,6,10,11,13)){
  for (p in pathway_s){
    patient_pathway1=c(patient_pathway1,plot_sup[[p]][i])
  }
}
pdf("sup8_1.pdf",width = c(16),height = c(10))
ggarrange(plotlist = patient_pathway1,
          ncol = 8,nrow = 5,align = "h")
dev.off()
pdf("sup8_2.pdf",width = c(16),height = c(16))
ggarrange(plotlist = patient_pathway,
          ncol = 8,nrow = 8,align = "h")
dev.off()
fig_a=NULL
fig_c=NULL
fig_d=NULL
pdf("pic1.pdf",width = 15,height = 8)

for (i in unique(melt_DT$pathway)){
  melt_DT_sel=melt_DT[melt_DT$pathway==i,]
  pcolor=c(RColorBrewer::brewer.pal(12,"Paired"),"darkolivegreen1")
  names(pcolor)=unique(melt_DT_sel$patient_paper)
  fig_a[[i]]=ggplot(melt_DT_sel)+
    geom_violin(aes(y=value,x=patient_paper,fill=patient_paper),size=0.2,
                scale="width",adjust=1.2)+
    scale_fill_manual(values=pcolor)+xlab("Patient")+ylab("Score")+
    geom_boxplot(aes(x=patient_paper,y=value),width=0.15,
                 outlier.size = 0.5,size=0.15)+
    stat_compare_means(aes(y=value,x=patient_paper),
                       method = "anova",size=2,hjust=-1.7)+
    theme_classic()+ggtitle(names(fill)[fill==i])+
    theme(legend.position = "none",text = element_text(size = 7),
          plot.title = element_text(hjust = 0.5),
          axis.title = element_text(),
          strip.text = element_text(),
          strip.background = element_rect(size=0.1),
          line = element_line(size = 0.1),
          legend.box.spacing = unit(0.01,"cm"),
          axis.ticks.length = unit(0.03,"cm"),
          axis.text.x = element_text(angle = 45,vjust=0.5))
  fig_c[[i]]=ggplot(melt[which(melt$pathway==i),])+
    geom_violin(aes(y=value,x=subtype,fill=subtype))+
    geom_boxplot(aes(x=subtype,y=value),width=0.2,outlier.size = 0.5,size=0.15)+
    stat_compare_means(aes(y=value,x=subtype),size=1.5,hjust=-4)+
    scale_fill_manual(values=c(RColorBrewer::brewer.pal(7,"Set2")[c(5,1,2,4)]))+
    theme_classic()+ggtitle(names(fill)[fill==i])+
    theme(legend.position = "none",text = element_text(size=7),
          plot.title = element_text(hjust = 0.5),
          axis.title = element_text(),
          strip.text = element_text(),
          strip.background = element_rect(size=0.1),
          line = element_line(size = 0.1),
          legend.box.spacing = unit(0.01,"cm"),
          axis.ticks.length = unit(0.03,"cm"),
          axis.text.x = element_text())
  
  fig_d[[i]]=ggarrange(plotlist = plot[[i]][c(5,6,10)],
                       ncol = 3)
  
  
  
  f_t1=ggarrange(fig_a[[i]],fig_c[[i]],ncol = 2,widths = c(2,1))
  f_t2=ggarrange(f_t1,fig_d[[i]],nrow = 2)
  print(f_t2)
  
}
dev.off()


sum_squres=function(x){sum((x-mean(x))**2)}
SS=sum_squres(melt_DT$value)
sse=melt_DT%>%group_by(pathway)%>%summarise(sum_squres(value))
ss_s=(melt_DT$value-mean(melt_DT$value))**2

melt_DT$sss=ss_s
ssp=melt_DT%>%group_by(pathway)%>%summarise(sum(sss))
f_val=(ssp$`sum(sss)`-sse$`sum_squres(value)`)/sse$`sum_squres(value)`

pdf("pathway_score1.pdf",width = 10,height = 8)
ComplexHeatmap::Heatmap(gsva_res,cluster_rows = F,cluster_columns = F)

mean_value=gsva_res%>%as.data.frame()%>%mutate(pathway=rownames(.))%>%
  gather(key = "sample",value = "value",-pathway)%>%
  mutate(patient=substr(sample,1,4))%>%
  group_by(patient,pathway)%>%summarise(mean=mean(value))

mean_sd=gsva_res%>%as.data.frame()%>%mutate(pathway=rownames(.))%>%
  gather(key = "sample",value = "value",-pathway)%>%
  mutate(patient=substr(sample,1,4))%>%
  group_by(patient,pathway)%>%summarise(var=var(value))

mean_value_s=mean_value%>%spread(.,patient,mean)%>%column_to_rownames("pathway")

library(ComplexHeatmap)
Heatmap(t(mean_value_s),cluster_rows = F,cluster_columns = F, name = "mean",
        column_title = "Patients' mean value",
        bottom_annotation = 
          columnAnnotation(var=f_val,
                           col=list(var=circlize::colorRamp2(c(min(na.rm = T,f_val),
                                                               max(na.rm = T ,f_val)/2,
                                                               max(na.rm = T ,f_val)),
                                                             c("blue", "white", "red")))))

mean_sd_s=mean_sd%>%spread(.,patient,var)%>%column_to_rownames("pathway")
ComplexHeatmap::Heatmap(scale(mean_sd_s),cluster_rows = F,cluster_columns = F,
                        name = "variance",column_title = "Patients' variance")

dev.off()

tmp=melt_DT_sel[which(melt_DT$pathway=="MINGUEZ_LIVER_CANCER_VASCULAR_INVASION_UP"),]
t=summary(lm(tmp$value~tmp$patient_paper))
t$r.squared
t$r.squared

r=melt_DT_target_pathway%>%group_by(pathway)%>%
  summarise(r=1-summary(lm(value~patient_paper))$r.squared)
spread_DT_target_pathway=melt_DT_target_pathway%>%select(pathway,key,value)%>%
  spread(key = pathway,value = value)
mean_sd_select=melt_DT_target_pathway%>%group_by(patient_paper,pathway)%>%
  summarise(var=var(value))%>%spread(.,patient_paper,var)%>%column_to_rownames("pathway")%>%
  t(.)%>%scale(.)


column_order=c("Vascular invasion","Metastasis","GEP score",
               "Inflammatory","Immune class","Cytolytic score","IMPRES","TIDE")

mean_sd_select=mean_sd_select[,column_order]
subtype_DT_t=read.table("../../3/class.tsv")
m=matrix(data = "",nrow = 13,ncol = 4)

subtype_DT_prop=subtype_DT_t%>%mutate(patient=substr(V1,1,4))%>%
  filter(substr(V1,1,2)=="DT")%>%mutate(V2=as.character(V2))%>%
  select(patient,V2)
rownames(m)=sort(unique(subtype_DT_prop$patient))
colnames(m)=c(1,2,3,4)
for (p in rownames(m)){
  tmp=subtype_DT_prop[which(subtype_DT_prop$patient==p),]
  a=table(tmp)/nrow(tmp)
  m[rownames(a),colnames(a)]=a
}
m=matrix(as.numeric(m),ncol = 4)
m[is.na(m)]=0

subtype_anno = rowAnnotation("Subtype" = anno_barplot(m,gp = gpar(fill = c(RColorBrewer::brewer.pal(7,"Set2")[c(5,1,2,4)])), 
                                                      axis = F,
                                                      bar_width = 1,width = unit(5, "mm")),
                             annotation_name_gp= gpar(fontsize =6),annotation_name_rot = 45)


IBD_DNA_R2=c(0,	0.06,	0.19,	0.12,	0.25,	0.435,	0.2,	
             0.169,	0,	0.24,	0.268,	0.539,	0.339)
IBD_RNA_R2=c(0.16,	0.11,	0.14,	0.14,	0.2,	0.364,	0.384,
             0.258,	0.0345,	0.0755,	0.281,	0.335,	0.441)
DNA_RNA_R2=c(0.21,	0.64,	0.43,	0.81,	0.94,	0.72,	0.566,	
             0.352,	0.204,	0.359,	0.736,	0.634,	0.578)


foo=HeatmapAnnotation(foo = anno_block(gp = gpar(fill = 2:4)),
                      annotation_height  = unit(2,"mm"))
var=ComplexHeatmap::Heatmap(mean_sd_select,cluster_rows = F,cluster_columns = F,
                            name = "z-score",column_split = c(1,1,2,2,2,2,2,2),
                            column_gap = unit(0, "mm"), border = TRUE,
                            #right_annotation = right_annotation,
                            left_annotation = subtype_anno,
                            top_annotation = c(foo),
                            column_title = c("Diagnose","Treatment"),
                            heatmap_legend_param = list(title_gp = gpar(fontsize =5),
                                                        labels_gp = gpar(fontsize =5),
                                                        grid_width = unit(1, "mm"),
                                                        legend_height=unit(1,"mm")),
                            row_title_gp = gpar(fontsize=6),
                            column_names_gp = gpar(fontsize=6),
                            column_names_rot = 45,
                            row_names_gp = gpar(fontsize=6),
                            column_title_gp = gpar(fontsize=6))
lgd = Legend(at = c("Wnt","Classic","Hybrid","Proliferation"), 
             title = "Transcritomic\nSubtype", 
             grid_height = unit(1, "mm"),
             grid_width = unit(1, "mm"),labels_gp = gpar(fontsize = 5),
             title_gp = gpar(fontsize = 5),
             legend_gp = gpar(fill = c(RColorBrewer::brewer.pal(7,"Set2")[c(5,1,2,4)]))
             )
plot_var=grid::grid.grabExpr(draw(var,padding = unit(c(2, 7, 2, 2), "mm"),
                                  merge_legends=T,gap=unit(0,"mm"),
                                  heatmap_legend_list = list(lgd),
                                  column_title="Variance of biomarkers",
                                  column_title_gp = gpar(fontsize=8)))

gsva_res_select=gsva_res[pathway_s,]
rownames(gsva_res_select)=names(fill)[match(rownames(gsva_res_select),fill)]
split=factor(subtype_vec[colnames(gsva_res_select)],levels = c("Wnt","Classic","Hybrid","Proliferation"))
heatmap_sel= ComplexHeatmap::Heatmap(t(scale(t(gsva_res_select))),cluster_rows = T,
                                     cluster_columns = T,cluster_column_slices = F,
                                     name = "z-score",column_split = split,show_row_dend = F,
                                     column_gap = unit(0, "mm"), border = TRUE,show_column_names = F,
                                     show_column_dend = F,
                                     heatmap_legend_param = list(title_gp = gpar(fontsize =5),
                                                                 labels_gp = gpar(fontsize =5),
                                                                 grid_width = unit(1, "mm"),
                                                                 legend_height=unit(1,"mm")),
                                     row_title_gp = gpar(fontsize=6),
                                     column_names_gp = gpar(fontsize=6),
                                     row_names_gp = gpar(fontsize=6),
                                     column_title_gp = gpar(fontsize=6))
heatmap_sel= grid.grabExpr(draw(heatmap_sel))

melt_sel=gsva_res_select%>%as.data.frame(.)%>%gather()
melt_sel=cbind(melt_sel,pathway=rownames(gsva_res_select))
melt_sel$subtype=factor(subtype_vec[melt_sel$key],
                        levels = c("Wnt","Classic","Hybrid","Proliferation"))

two.way <- aov(value ~ subtype + pathway, data = melt_sel)


p1_sel=ggplot(melt_sel)+geom_violin(aes(y=value,x=subtype,fill=subtype))+
  geom_boxplot(aes(x=subtype,y=value),width=0.2,outlier.size = 0.5,size=0.15)+
  stat_compare_means(aes(y=value,x=subtype))+
  scale_fill_manual(values=c(RColorBrewer::brewer.pal(7,"Set2")[c(5,1,2,4)]))+
  facet_wrap(~pathway,scales = "free",nrow = 6,ncol = 4)+theme_bw()

i="MINGUEZ_LIVER_CANCER_VASCULAR_INVASION_UP"
fig_d=ggarrange(plotlist = list(heatmap_sel,plot[["GEP score"]][5]$p,
                                plot[["Cytolytic_score"]][5]$p),
                ncol = 3,labels = c("c","d","e"),font.label = list(size=8),
                widths = c(3.77,2,2))

tmp=ggarrange(fig_a[[i]],plot_var,ncol  = 2,widths = c(1,1.1),
              labels = c("a","b"),font.label = list(size=8))

pdf("fig6_v6_a4.pdf",width = 8.27,height =5)
ggarrange(tmp,fig_d,nrow =2,heights = c(1.5,1))+
  theme(plot.margin = unit(c(0.1,0.25,0.1,0.25),"inch"))
dev.off()


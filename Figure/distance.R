library(tidyverse)

for (patient in substr(list.files(pattern = "DT.*3000.*"),1,4)){
  RNA=read.table(paste0(patient,"_DESeq2.log2_top3000CodingGenes.txt"),header = T)
  dist = 1-cor(RNA, method="spearman")
  #dist=as.matrix(dist(t(RNA),method = "euclidean",diag = TRUE,upper = TRUE))
  RNA_distance=dist%>%as.data.frame()%>%mutate(sample1=rownames(dist))%>%
    gather(sample2,value,-sample1)%>%arrange(sample1)
  sample_list=c();tmp=c()
  for (r in 1:nrow(RNA_distance)){
    if (paste0(RNA_distance[r,1:2],collapse = ":") %in% sample_list|RNA_distance[r,3]==0){
      tmp=c(tmp,r)
    }else{
      sample_list=c(sample_list,paste0(RNA_distance[r,2:1],collapse = ":"))
    }
  }
  RNA_distance=RNA_distance[-tmp,]
  write.table(RNA_distance,file = paste0(patient,"_RNA_dis.tsv"),row.names = F,quote = F)
}

physical=read.csv("E:/cancer genome/liver/·ÖÎö/distance/all.csv",header = T)
for (p in unique(physical$Patient)){
  
  select_data=physical%>%filter(Patient==p)%>%select(Sample,X,Y)
  
  pos_dist=c();sample1=c();sample2=c()
  for (i in 1:nrow(select_data)){
    if (i != nrow(select_data)){
      for (ii in (i+1):nrow(select_data)){
        pos_dist=c(pos_dist,sqrt((select_data$X[i]-select_data$X[ii])**2+
                                   (select_data$Y[i]-select_data$Y[ii])**2))
        sample1=c(sample1,select_data$Sample[i])
        sample2=c(sample2,select_data$Sample[ii])
      }
    }
  }
  res=cbind.data.frame(sample1,sample2,pos_dist)
  write.table(res,file=paste0(p,"_physical_dis.tsv"),row.names = F,quote = F)
}


library(lsa)

clone_files=list.files("E:/cancer genome/liver/analysis/mapscape/citup1/citup_all/",
                       pattern = "*cluster_prev.txt",full.names = T)
for (n in 1:length(clone_files)){
  patient=sub("_.*$","", sub("^.*/","",clone_files[n]))
  samples=read.table(paste0("E:/cancer genome/liver/analysis/mapscape/citup1/citup_all/",
                            patient,"_sample_order.txt"))[,1]
  clone_prev=read.table(clone_files[n])
  colnames(clone_prev)=samples
  
  comb=combn(samples,2)
  
  clone_dist=NULL
  for (col in 1:ncol(comb)){
    sample_pair=comb[,col]
    d=dist(t(clone_prev[sample_pair]))
    #d=1-cosine(clone_prev[,sample_pair[1]],clone_prev[,sample_pair[2]])
    clone_dist=rbind(clone_dist,c(sample_pair,d))
  }
  
  colnames(clone_dist)=c("sample1","sample2","clone.dist")
  write.table(clone_dist,file = paste0(patient,"_clone_dist.txt"),quote = F,
              col.names = T,row.names = F)
  
}

maf.files=list.files("../tree/",pattern = "*maf",full.names = T)
for (n in 1:length(maf.files)){
  maf_data=read.table(maf.files[n],header = T)
  patient=sub(".maf","",sub("^.*/","",maf.files[n]))
  tumorNames = colnames(maf_data)[-which(colnames(maf_data) %in% c("mut_id","Variant_Classification"))]
  count_mat = maf_data %>% select(tumorNames) %>% mutate_at(c(tumorNames), as.numeric)
  count_mat = as.data.frame(count_mat)
  rownames(count_mat) = maf_data$mut_id
  binary_mat = count_mat
  binary_mat[binary_mat!=0] = 1
  samples=colnames(binary_mat)
  
  comb=combn(samples,2)
  
  ham_dist=NULL
  for (col in 1:ncol(comb)){
    sample_pair=comb[,col]
    d=length(which(binary_mat[,sample_pair[1]]!=binary_mat[,sample_pair[2]]))
    ham_dist=rbind(ham_dist,c(sample_pair,d))
  }
  colnames(ham_dist)=c("sample1","sample2","ham.dist")
  write.table(ham_dist,file = paste0(patient,"_ham_dist.txt"),quote = F,
              col.names = T,row.names = F)
  
}

tree_files=list.files("../tree/",pattern = "*newick.tre$",full.names = T)
library(ape)

for (n in 1:length(tree_files)){
  patient=sub("newick.*$","", sub("^.*/","",tree_files[n]))
  tree=read.tree(tree_files[n])
  par_tree=read.tree(paste0(patient,"par.tre"))
  nj_dist=cophenetic(tree)
  par_dist=cophenetic(par_tree)
  samples=tree$tip.label
  comb=combn(samples,2)
  
  nj_tree_dist=NULL;par_tree_dist=NULL
  for (col in 1:ncol(comb)){
    sample_pair=comb[,col]
    d=nj_dist[sample_pair[1],sample_pair[2]]
    
    d1=par_dist[sample_pair[1],sample_pair[2]]
    nj_tree_dist=rbind(nj_tree_dist,c(sample_pair,d))
    par_tree_dist=rbind(par_tree_dist,c(sample_pair,d1))
  }
  
  colnames(nj_tree_dist)=c("sample1","sample2","nj.dist")
  colnames(par_tree_dist)=c("sample1","sample2","par.dist")
  
  write.table(nj_tree_dist,file = paste0(patient,"_nj_dist.txt"),quote = F,
              col.names = T,row.names = F)
  write.table(par_tree_dist,file = paste0(patient,"_par_dist.txt"),quote = F,
              col.names = T,row.names = F)
}
plot(nj_tree_dist[,3],par_tree_dist[,3])

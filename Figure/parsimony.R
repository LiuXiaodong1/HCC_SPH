library(ape)
library(tidyverse)
library(phangorn)
set.seed(123)
plot_DNA_tree = function(maf){
  maf_data=read.table(maf,header = T)
  tumorNames = colnames(maf_data)[-which(colnames(maf_data) %in% c("mut_id","Variant_Classification"))]
  maf_data$N = rep(0, nrow(maf_data))
  count_mat = maf_data %>% select(tumorNames, "N") %>% mutate_at(c(tumorNames, "N"), as.numeric)
  count_mat = as.data.frame(count_mat)
  rownames(count_mat) = maf_data$mut_id
  binary_mat = count_mat
  sequence_mat=count_mat
  alt=unlist(lapply(str_split(rownames(count_mat),"_"),function(x)x[[6]]))
  ref=unlist(lapply(str_split(rownames(count_mat),"_"),function(x)x[[5]]))
  binary_mat[binary_mat!=0] = 1

  sequence_mat=apply(sequence_mat,2,function(x){
    x[which(x!=0)]=alt[which(x!=0)];
    x[which(x==0)]=ref[which(x==0)];return(x)})
  sequence_list=apply(sequence_mat,2,function(x){paste0(x,collapse = "")})
  seq_file=sub("^.*/","",sub(".maf","par.fasta",maf))
  if (file.exists(seq_file)){
    file.remove(seq_file)
  }
  for (l in names(sequence_list)){
    seqname=paste0(">",l)
    line=c(seqname,sequence_list[[l]])
    write_lines(line,file = seq_file,append = T)
  }
  
  
  mut_count = binary_mat %>% transmute(count=rowSums(.))
  trunk_length = mut_count %>% filter(count == length(tumorNames)) %>% nrow()
  private_mut = mut_count %>% filter(count == 1) %>% nrow()
  binary_mat=t(binary_mat)
  
  
  phyDat=phyDat(t(sequence_mat))
  trw=NJ(dist.hamming(phyDat,ratio = F))
  rtrw=root(trw,outgroup = "N", resolve.root = T)

  par_tree=pratchet(phyDat,start = rtrw, method = "fitch", maxit = 10000,
                    minit = 10, k = 100, trace = 1, all = FALSE,
                    rearrangements = "SPR",
                    perturbation = "ratchet")
  par_tree_ratchet=acctran(par_tree,phyDat)
  rpar = root(par_tree_ratchet,outgroup = "N",resolve.root = T)

  treeSPR <- optim.parsimony(rtrw, phyDat)
  treeSPR <- acctran(treeSPR,phyDat)
  #rpar=root(treeSPR,outgroup="N",resolve.root=T)
  parsimony(c(trw,rtrw,rpar,par_tree_ratchet,treeSPR),phyDat)
  
  write.tree(rpar,file = sub("^.*/","",sub(".maf","par.tre",maf)))

}

for (maf in paste0("../tree/",c("DT03","DT04","DT06","DT07","DT09","DT10",
              "DT12","DT13","DT14","DT16","DT17","DT18","DT19"),".maf")){
  plot_DNA_tree(maf)
}


pdf(paste0("compare_trees.pdf"),width = 16,height = 8)
for (patient in c("DT03","DT04","DT06","DT07","DT09","DT10",
                  "DT12","DT13","DT14","DT16","DT17","DT18","DT19")){

par_tree=read.tree(paste0(patient,"par.tre"))
NJ_tree=read.tree(paste0("../tree/",patient,"newick.tre"))

compare1=phytools::cophylo(NJ_tree,par_tree,rotate = T)


plot(compare1)
mtext("NJ tree",side = 2,line = -1)
mtext("Parsimony tree",side = 4,line = -1)
}
dev.off()


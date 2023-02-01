suppressPackageStartupMessages(library(mapscape))
suppressPackageStartupMessages(library(htmlwidgets))
suppressPackageStartupMessages(library(tidyverse))

args = commandArgs(trailingOnly = TRUE)

patient <- args[1]
qip <- args[2]
square_dir <- args[3]

options(viewer = NULL)

tree <- list.files("./", "*tree_[0-9]+.adj")
tree_num <- gsub(".*tree_([0-9]+).adj", "\\1", tree)

clonal_prev <- read_tsv(paste0(patient, "_tree_", tree_num, "_clone_freq.txt"),
                        col_types = cols())
clonal_prev <- clonal_prev %>% rename("sample_id" = "sample")

sample_order <- read_tsv(paste0(patient, "_sample_order.txt"), col_names = F,
                         col_types = cols())
sample_order$X2 <- as.integer(rownames(sample_order)) - 1
colnames(sample_order) <- c("sample_name", "sample_id")

# Set sample names
clonal_prev$sample_id <- sapply(clonal_prev$sample_id, 
                               function(x) sample_order[sample_order$sample_id == x, 'sample_name', drop=T])

tree <- read_tsv(paste0(patient, "_tree_", tree_num, ".adj"), col_names = F,
                 col_types = cols())
colnames(tree) <- c("source", "target")

# Backmapping of clusters
if (qip=="qip"){
	var_assignment <- read.delim(paste0(patient, "_tree_", tree_num, "_cluster_assignment.txt"),
				     header = F)
} else {
	var_assignment <- read.delim(paste0(patient, "_tree_", tree_num, "_variant_assignment.txt"),
				     header = F)
}
var_assignment$V2 <- rownames(var_assignment)
colnames(var_assignment) <- c("subclone_assigned", "original_cluster_id")

# Get mutations
mutations <- read_tsv(paste0(patient, "_loci_filtered.tsv"), col_types = cols())
 
# Assign mutations to new clone ID
mutations$cluster_id <- sapply(mutations$cluster_id, 
                               function(x) var_assignment[var_assignment$original_cluster_id == x, 'subclone_assigned', drop=T])

mutations <- separate(mutations, col=mutation_id, sep=":|_", into=c("Hugo_Symbol", "chrom", "coord"))
mutations <- rename(mutations, "VAF" = "variant_allele_frequency", "clone_id" = "cluster_id")

#original
#x = sample(50:150, nrow(sample_order))
#y = sample(50:150, nrow(sample_order))

#in vertical line
#x = rep(100, nrow(sample_order))
#y = seq(30,170,length=nrow(sample_order))

#in horizontal line
x = seq(30,170,length=nrow(sample_order))
y = rep(100, nrow(sample_order))

#binxin's version
#interval = 100/nrow(samples)
#x=list()
#y=list()
#for (i in 1:nrow(samples))
#{
## The samples are in a horizontal line
#    x[[i]]=50+interval*(i-1)
#    y[[i]]=100
#}

 
sample_info<-read.csv(paste0(square_dir,"/all.csv"))
sample_locations<-sample_info[which(sample_info$Patient==patient),][,c(1,11,12)]
sample_locations$X<-sample_locations$X * 100
sample_locations$Y<- sample_locations$Y * 100
sample_locations$sample_id<-paste0(sample_locations$Sample,"_T")
sample_locations<-sample_locations%>%mutate(Sample=NULL)%>%
  mutate(location_id=sample_id)%>%rename(x=X)%>%rename(y=Y)
#sample_locations <- data_frame(sample_id=unique(clonal_prev$sample_id),
#                               location_id=sample_order$sample_name,
#                               x=x, y=y)

#img_ref <- "/mnt/projects/lailhh/workspace/Liver_TCR/PyClone/PyClone_package/pyclone_nf_pipelines/squares.png"
#img_ref <- "../scripts/squares.png"
img_ref <- paste0(square_dir,"/squares.png")

clonal_prev$sample_id<-substr(clonal_prev$sample_id,3,7)
sample_locations$sample_id<-substr(sample_locations$sample_id,3,7)
mutations$sample_id<-substr(mutations$sample_id,3,7)
savetree <- mapscape(clonal_prev = clonal_prev, tree_edges = tree,
                     sample_locations = sample_locations, img_ref = img_ref)

saveWidget(savetree, paste0(patient, "_mapscape.html"),selfcontained = F)


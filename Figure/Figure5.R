library(ComplexHeatmap)

data=matrix(data=1,nrow = 3,ncol = 13)
colnames(data)=c(paste0("SH0",1:9),paste0("SH",10:13))
data[1,]=c(1,0,1,1,1,0,1,0,1,0,1,1,1)
data[2,]=c(1,1,0,1,1,1,1,1,1,1,1,1,1)
data[3,]=c(0,1,0,0,1,1,1,1,1,0,1,1,1)

data=as.data.frame(data)

data=sapply(data, function(x) as.character(x))
rownames(data)=c("Imbalance DNA tree",
                 "Imbalance RNA tree",
                 "MOBSTER test")

pdf("Figure5d.pdf",width = 2.95,height = 1.96)
p=Heatmap(data,cluster_rows = F,cluster_columns = F,
          column_split = c(2,2,2,2,1,1,2,2,2,1,2,1,2),
          rect_gp = gpar(col = "black", lwd = 0.2),
          row_names_gp = gpar(fontsize=5),
          column_names_gp = gpar(fontsize=5),
          #width = unit(6,"cm"),height = unit(5,"cm"),
          column_title = c("SCD","Non-SCD"),
          column_title_gp = gpar(fontsize=6),
          col = c("#377EB8","#E41A1C"),
          heatmap_legend_param = list(
            at = c(1, 0),
            title="",grid_width = unit(2, "mm"),
            labels_gp = gpar(fontsize = 5),
            nrow=1,
            #legend_direction="horizontal",
            labels = c("Significant","Not significant"),
            legend_height = unit(4, "cm")
          ))
draw(p,heatmap_legend_side = "bottom",
     annotation_legend_side = "bottom",merge_legend=T)
dev.off()



#######################################
library(Seurat)
library(ggplot2)

source("/ceph/project/holab/jwoo/10X_PW/Xenium/scripts/niche_functions.R")


output_folder <- "/merged_Niche_Analysis_expression_K35/"

seurat<-readRDS("/post_merge_obj.RDS")

seurat<-JoinLayers(seurat,assay = "XENIUM")

cnt<-seurat@assays$XENIUM$counts

##### Convert Cnt into assay object 3 
options(Seurat.object.assay.version = "v3")
seurat_n <- CreateAssayObject(counts = cnt)

seurat[["XENIUM3"]]<-seurat_n


seurat@active.ident<-as.factor(seurat$PerSlide_prediction)




options(future.globals.maxSize = 8000 * 1024^2)

setwd(output_folder)

seurat$Cluster <- seurat@active.ident
###niche analysis - expression based
niche.exprs <- BuildNicheExpressionAssay(seurat, neighbors.k = 35, assay.counts = "XENIUM3")
seurat <- FindNiches(seurat, t(niche.exprs), niche.name = "Expr.Niches_k35", niche.assay.name = "niche.exprs_k35")


df <- CellTypeNicheEnrichment(seurat$Expr.Niches_k35, seurat$Cluster)
cairo_pdf("niche_cell_type_enrichment_barplots_per_niche_expression_method.pdf", onefile = T)

for( niche in levels( unique(seurat$Expr.Niches_k35)) ){
  
  print(PlotCellTypeNicheEnrichmentBarPlot(df, niche=niche))
  
}
dev.off()

cairo_pdf("niche_cell_type_enrichment_heatmap_expression_method.pdf")
print(PlotCellTypeNicheEnrichmentHeatmap(df))
dev.off()


write.csv(df,"Enrichment_test_k35.csv")

cairo_pdf("niche_overlay_spatial_slide_expression_method_k35.pdf", onefile = T)
for( image in Images(seurat)){
  
  print(ImageDimPlot(seurat,fov=image, border.size = NA, group.by = "Expr.Niches_k35") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))
  
}
dev.off()



#############################################################
Idents(seurat) <- "Expr.Niches_k35"
seurat@active.assay <- "SCT"
mk <- FindAllMarkers(seurat, max.cells.per.ident = 500)
genes <- lapply(unique(seurat$Expr.Niches_k35), FUN=function(x){ dplyr::filter(mk, cluster == x & avg_log2FC > 0 &
                                                                                 !startsWith(gene, "mt-"))[1:6, "gene"]})
genes <- unique(unlist(genes))
genes <- genes[!is.na(genes)]
plot <- Seurat::DotPlot(seurat, features=genes, group.by = "Expr.Niches_k35")
mat <- reshape2::acast(plot$data[, 3:5], id~features.plot, value.var="avg.exp.scaled")
row_order <- hclust(dist(mat))$order
col_order <- hclust(dist(t(mat)))$order
plot$data$features.plot <- factor(plot$data$features.plot, levels=levels(plot$data$features.plot)[col_order])
plot$data$id <- factor(plot$data$id, levels=levels(plot$data$id)[row_order])

cairo_pdf("bubbleplot_niche_marker_genes_expression_method.pdf", height=5, width=14)
print(ggplot(plot$data, aes(id, features.plot, fill=avg.exp.scaled, size=pct.exp)) + geom_point(shape=21
) + scale_fill_distiller(palette="BrBG") + theme_light(base_size = 10
) + labs(x="", y="", fill="Average Exp.",
         size="Percent\nPositive")  + theme(axis.text.x = element_text(angle=90)) +coord_flip())

dev.off()

write.table(mk, file="niche_marker_genes_expression_method.xls", sep="\t", quote=FALSE, row.names = FALSE)

setwd(output_folder)
saveRDS(seurat, file="/ceph/project/holab/jwoo/IPF_HT/Xenium/post_niche_analysis_seurat.RDS")
#######################################################


############## Run this bit later for color adjust for spatial niche color ####
library(scales)

seurat <- readRDS("/ceph/project/holab/jwoo/IPF_HT/Xenium/post_niche_analysis_seurat.RDS")


show_col(hue_pal()(8))
hue_pal()(8)

my_cols <- c('0'="#F8766D",'1'="#CD9600" ,'2'="#7CAE00",'3'="#00BE67",'4'="#00BFC4" ,'5'="#00A9FF",'6'="#C77CFF",'7'="#FF61CC" ,'8'="#E76BF3",'9'="#FF62BC")

cairo_pdf("niche_overlay_spatial_slide_expression_method_k35.pdf", onefile = T)
for( image in Images(seurat)){
  
  print(ImageDimPlot(seurat,fov=image, border.size = NA, group.by = "Expr.Niches_k35",cols=my_cols) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))
  
}
dev.off()


#############. Niche by Niche Spatial Plot ############
setwd("/merged_Niche_Analysis_expression_k35/Niche_Highlight")

Niches<-levels(unique(seurat$Expr.Niches_k35))
for (i in 1:length(Niches)) {
  seurat$Niche_ID<-seurat$Expr.Niches_k35
  seurat$Niche_ID<-ifelse(seurat$Niche_ID ==Niches[i],Niches[i],NA)
  
  cairo_pdf(paste0(Niches[i],"_overlay_spatial_slide_expression_method_k35.pdf"), onefile = T)
  for( image in Images(seurat)){
    
    print(ImageDimPlot(seurat, fov=image, group.by = "Niche_ID",border.size = NA,cols=my_cols[i]) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))
    
  }
  dev.off()
  
}


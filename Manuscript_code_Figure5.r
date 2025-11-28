
library(Seurat)
library(ggplot2)
library(scales)


source("/niche_functions.R")

################# Fig 5B. Niche Overlay on spatial ######################
# Load seurat object after niche identification
seurat<-readRDS("IPF_xenium.RDS")


show_col(hue_pal()(8))
hue_pal()(8)

my_cols <- c('0'="#F8766D",'1'="#CD9600" ,'2'="#7CAE00",'3'="#00BE67",'4'="#00BFC4" ,'5'="#00A9FF",'6'="#C77CFF",'7'="#FF61CC" ,'8'="#E76BF3",'9'="#FF62BC")

cairo_pdf("niche_overlay_spatial_slide_expression_method_k35.pdf", onefile = T)
for( image in Images(seurat)){
  
  print(ImageDimPlot(seurat,fov=image, border.size = NA, group.by = "Expr.Niches_k35",cols=my_cols) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))
  
}
dev.off()

################# Fig 5C. Niche Overlay on spatial (Niche highlight) ######################

setwd("/Niche_Highlight")

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


################# Fig 5D. Celltype Enrichment BarPlot ######################


df <- CellTypeNicheEnrichment(seurat$Niche_ID, seurat$Celltype)
cairo_pdf("niche_cell_type_enrichment_barplots_per_niche_expression_method.pdf", onefile = T)

for( niche in levels( unique(seurat$Niche_ID)) ){
  
  print(PlotCellTypeNicheEnrichmentBarPlot(df, niche=niche))
  
}
dev.off()

################ Fig 5E. Celltype Enrichment Heatmap  ######################

df2<-df[df$Cell.Type %in% c("IMs","Mast Cells","ncMono","cMono","NK Cells","MonoMac","cDC2","Cycling cells","CD8 T cells","CD4 T cells","CD103posCD4 T cells","Tregs","CXCL10+ Mono-Macs","cDC1","Plasma Cells","SPP1 Macs","pDCs","B Cells","Intm AM","FABP4 AM","CD206hi FN1hi AM"),]

niche.results<-df2

niche.results$Estimate <- as.numeric(niche.results$Estimate)
mat <- log2(reshape2::acast(niche.results, `Cell.Type`~ Niche, value.var = "Estimate")+ 0.01)
niche.results$STAR <- "n.s"; niche.results$STAR[niche.results$FDR < 0.05] <- "*"; niche.results$STAR[niche.results$FDR < 0.01] <- "**"; niche.results$STAR[niche.results$FDR < 0.001] <- "***"
mat.labels <-reshape2::acast(niche.results, `Cell.Type`~ Niche, value.var = "STAR")
mat[mat < -4] <- -4; mat[ mat > 4] <- 4
pheatmap::pheatmap(mat, border_color = "black", display_numbers = mat.labels, fontsize_number = 14,
                   color=Seurat:::SpatialColors(n=100), number_color = "black" )


cairo_pdf("niche_cell_type_enrichment_heatmap_Immune_Grp1.pdf")
print(pheatmap::pheatmap(mat, border_color = "black", display_numbers = mat.labels, fontsize_number = 14,
                         color=Seurat:::SpatialColors(n=100), number_color = "black" )
)
dev.off()


################ Fig 5F. T1_IFN_Module Score group overlay  ######################


genes<- c("CXCL10", "CXCL9", "ISG15","MX1","IFIT1","IFIT3","RSAD2","STAT1","STAT2","IRF9","IRF1")

# module score calculation
seurat <- AddModuleScore(
  object = seurat,
  features = list(genes),
  ctrl = 10,
  name = 'T1_IFN_ModuleScore'
)

# Pick cells with upper 10% module scores 

meta<-seurat@meta.data

meta1<-meta[meta$SAMPLE_ID=="S1",]
meta$T1_IFN_ModuleScore2<-ifelse(meta$T1_IFN_ModuleScore1 > quantile(meta1$T1_IFN_ModuleScore1,0.9) & meta$SAMPLE_ID=="S1" ,"High","Low" )

meta1<-meta[meta$SAMPLE_ID=="S2",]
meta$T1_IFN_ModuleScore2<-ifelse(meta$T1_IFN_ModuleScore1 > quantile(meta1$T1_IFN_ModuleScore1,0.9) & meta$SAMPLE_ID=="S2" ,"High",meta$T1_IFN_ModuleScore2 )

meta1<-meta[meta$SAMPLE_ID=="S3",]
meta$T1_IFN_ModuleScore2<-ifelse(meta$T1_IFN_ModuleScore1 > quantile(meta1$T1_IFN_ModuleScore1,0.9) & meta$SAMPLE_ID=="S3" ,"High",meta$T1_IFN_ModuleScore2 )

meta1<-meta[meta$SAMPLE_ID=="S4",]
meta$T1_IFN_ModuleScore2<-ifelse(meta$T1_IFN_ModuleScore1 > quantile(meta1$T1_IFN_ModuleScore1,0.9) & meta$SAMPLE_ID=="S4" ,"High",meta$T1_IFN_ModuleScore2 )

meta1<-meta[meta$SAMPLE_ID=="S5",]
meta$T1_IFN_ModuleScore2<-ifelse(meta$T1_IFN_ModuleScore1 > quantile(meta1$T1_IFN_ModuleScore1,0.9) & meta$SAMPLE_ID=="S5" ,"High",meta$T1_IFN_ModuleScore2 )

seurat$T1_IFN_ModuleScore2<-meta$T1_IFN_ModuleScore2

seurat$T1_IFN_ModuleScore2<-factor(seurat$T1_IFN_ModuleScore2,levels = c("High","Low"))

cairo_pdf("T1_IFN_moduelscore_slides_top10pc.pdf", onefile = T)
for( image in Images(seurat)){
  
  print(ImageDimPlot(seurat, fov=image,border.size = NA,group.by = "T1_IFN_ModuleScore2",cols=c("#FDE725FF","#414487FF")))
  
}
dev.off()

#Highlight on individual niche

xen$T1_IFN_ModuleScore4<-xen$T1_IFN_ModuleScore2
xen$T1_IFN_ModuleScore4<-ifelse(!(xen$Expr.Niches_k35=="0"),NA,xen$T1_IFN_ModuleScore4)
xen$T1_IFN_ModuleScore4<-gsub("1","High",xen$T1_IFN_ModuleScore4)
xen$T1_IFN_ModuleScore4<-gsub("2","Low",xen$T1_IFN_ModuleScore4)

xen$T1_IFN_ModuleScore4<-factor(xen$T1_IFN_ModuleScore4,levels = c("High","Low"))

cairo_pdf("T1_IFN_moduelscore_slides_top10pc_Niche0.pdf", onefile = T)
for( image in Images(xen)){
  
  print(ImageDimPlot(xen, fov=image,border.size = NA,group.by = "T1_IFN_ModuleScore4",cols=c("#FDE725FF","#414487FF"))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) )
  
}
dev.off()




require(Seurat)
require(dplyr)
require(Matrix)
require(magrittr)
library(scales)
library(ggplot2)
library(configr)
library(cowplot)
library(Hmisc)
library(RColorBrewer)



################ Fig 2A. UMAP ##################


BAL<-readRDS("IPF_BAL.rds")

print(head(BAL@meta.data))
DefaultAssay(object = BAL) <- "RNA"
Idents(BAL)<-BAL@meta.data$CellType_v2
DimPlot(BAL, reduction = "umap", group.by="Celltype_v2",cols=c("#d8eab2","grey","grey52","#D62728","#E377C2","#9467BD","#FF9896","#C49C94","#2CA02C" ,"#BCBD22","#8C564B", "#F7B6D2", "#C5B0D5", "#FFBB78", "#AEC7E8" , "#98DF8A", "#1F77B4"),label=T,repel = T)
ggsave('IPF_BAL_UMAP.pdf', dpi=900,width=7, height=5, units="in")


################ Fig 2B. Top10 marker dotplot ##################

BAL1<-BAL[,Idents(BAL) %in% c("FABP4 AM","IGF1+ AM","CD206hi FNhi AM","monoSPP1+ AM","CXCL10+ AMs")]
markers <- FindAllMarkers(BAL1, max.cells.per.ident = 200)


top10<-markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top10_markers<-unique(top10$gene)

Idents(BAL1)<-factor(BAL1$CellType_v2,levels = c("FABP4 AM","IGF1+ AM","CD206hi FNhi AM","CXCL10+ AMs","monoSPP1+ AM"))
DotPlot(object = BAL1, assay="RNA", features  = top5_markers,cols="RdYlBu" ) + scale_y_discrete(limits = rev) +  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
    guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) + guides(size = guide_legend(title = "Percent\nExpressed"),color = guide_colorbar(title = "Average\nExpression"))+ theme(axis.text.x  =element_text( hjust=1, size=10,angle=90))

ggsave('Mac_Marker_DotPlot.pdf', dpi=900,width=10, height=4.1, units="in")

################ Fig 2D. Overlay of merged DEP on BAL data ##################

# This should be the spectra of GEP from 5 merged data
m_Spectra<-read.table("/BatchCorrected/Component_15/BatchCorrected.gene_spectra_score.k_15.dt_0_15.txt", sep='\t', row.names=1, header=TRUE)
m_Spectra_tpm<-read.table("BatchCorrected.gene_spectra_tpm.k_15.dt_0_15.txt", sep='\t', row.names=1, header=TRUE)



exp_mat<-LayerData(BAL, assay="RNA", layer='data')
common_gene<-intersect(colnames(m_Spectra),rownames(exp_mat))

length(common_gene)

exp_mat_1<-exp_mat[rownames(exp_mat) %in% common_gene,]
m_Spectra_1<-m_Spectra[,colnames(m_Spectra) %in% common_gene]
m_Spectra_1_tpm<-m_Spectra_tpm[,colnames(m_Spectra_tpm) %in% common_gene]

all(colnames(m_Spectra_1)==rownames(exp_mat_1))

factor<-t(exp_mat_1)  %*% t(as.matrix(m_Spectra_1))
colnames(factor)<-paste0("Merged_cNMF_",colnames(factor))


BAL<- AddMetaData(
  object = BAL,
  metadata = factor
)


FeaturePlot_scCustom(seurat_object = BAL, features = colnames(factor), num_columns = 2)
ggsave('merged_cNMF_FeaturePlot.pdf', dpi=900,width=8, height=12, units="in")



################ Fig 2E & F. Violin Plot of overlayed GEP ##################
# Insert usage mtx as assay
factor[factor<0]<-0
usage_norm_f<-t(as.matrix(factor))

BAL[["cNMF_mOverlay"]]<-CreateAssayObject(usage_norm_f)


BAL@active.assay<-"cNMF_mOverlay"

features<-rownames(BAL@assays$cNMF_mOverlay@counts)

BAL1<-BAL[,(BAL$CellType_v2 %in% c("mono-SPP1+ mac","CXCL10+ AM","FABP4hi AM","IGF1+ AM","CD206hiFN1hi AM","mono-cDC","mono-pDC"))]
# Set the order after clustering
BAL1$CellType_v2<-factor(BAL1$CellType_v2,levels=c("mono-SPP1+ mac","CXCL10+ AM","FABP4hi AM","IGF1+ AM","CD206hiFN1hi AM","mono-cDC","mono-pDC"))


BAL1@active.assay<-"cNMF_mOverlay"
Idents(BAL1)<-BAL1$CellType_v2



# same order as merged cNMF
features<-c("Usage1","Usage2","Usage4","Usage5","Usage7","Usage6","Usage9","Usage10")

### Save Vertical 8.5 x 8
plot <- VlnPlot(BAL1, features, stack = TRUE, sort = F, flip = TRUE,same.y.lims = T,y.max=1.0,cols=c("#272E6A","#272E6A","#272E6A","#272E6A","#272E6A","#BF0000","#BF0000","#BF0000","#89C75F","#89C75F","#89C75F","#89C75F")) +theme(axis.text.x = element_text(size = 15))
ggsave('merged_cNMF_VlnPlot.pdf', dpi=900,width=8, height=8.5, units="in")


features1<-c("Usage1","Usage2","Usage4","Usage5","Usage7","Usage6","Usage9","Usage10")

plot2 <- VlnPlot(BAL1, features1, stack = TRUE, sort = F, flip = TRUE,same.y.lims = T,split.by ="Disease" , split.plot=T,cols=c("#008494", "#F47D2B")) +theme(axis.text.x = element_text(size = 15))
ggsave('merged_cNMF_VlnPlot_DiseaseSplit.pdf', dpi=900,width=8, height=8.5, units="in")

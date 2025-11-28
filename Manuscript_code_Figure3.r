require(Seurat)
require(dplyr)
require(Matrix)
library(scales)
library(ggplot2)
library(configr)
library(cowplot)
library(Hmisc)
library(RColorBrewer)
library(sccomp)
library(loo)


######## Fig 3B. UMAP (Morse Myeloid subset)#########

Morse<-readRDS("Morse_Myeloid_MAY2024_f.rds")

print(head(Morse@meta.data))
DefaultAssay(object = Morse) <- "RNA"
Idents(Morse)<-Morse@meta.data$CellType_V4
DimPlot(Morse, reduction = "umap", group.by="CellType_V4",label=T,repel = T)
ggsave('Morse_Myeloid_UMAP.pdf', dpi=900,width=7, height=5, units="in")


######## Fig 3C. Overlay of cNMF (from merged) into Morse UMAP #########


# This should be the spectra from cNMF of 5 Merged data
m_Spectra<-read.table("/BatchCorrected/Component_15/BatchCorrected.gene_spectra_score.k_15.dt_0_15.txt", sep='\t', row.names=1, header=TRUE)
m_Spectra_tpm<-read.table("/BatchCorrected/Component_15/BatchCorrected.gene_spectra_tpm.k_15.dt_0_15.txt", sep='\t', row.names=1, header=TRUE)

Morse<-readRDS("../../Morse_Macs_MAR24.rds")

exp_mat<-LayerData(Morse, assay="RNA", layer='data')
common_gene<-intersect(colnames(m_Spectra),rownames(exp_mat))

length(common_gene)

exp_mat_1<-exp_mat[rownames(exp_mat) %in% common_gene,]
m_Spectra_1<-m_Spectra[,colnames(m_Spectra) %in% common_gene]
m_Spectra_1<-m_Spectra[,rownames(exp_mat_1)]
m_Spectra_1_tpm<-m_Spectra_tpm[,colnames(m_Spectra_tpm) %in% common_gene]

all(colnames(m_Spectra_1)==rownames(exp_mat_1))

factor<-t(exp_mat_1)  %*% t(as.matrix(m_Spectra_1))
colnames(factor)<-paste0("Merged_cNMF_",colnames(factor))


Morse<- AddMetaData(
  object = Morse,
  metadata = factor
)


FeaturePlot_scCustom(seurat_object = Morse, features = colnames(factor), num_columns = 4)
ggsave('merged_cNMF_on_Morse.pdf', dpi=900,width=16, height=10, units="in")



######## Fig 3D. Violin Plot of cNMF Score #########


# usage mtx as assay
factor[factor<0]<-0
usage_norm_f<-t(as.matrix(factor))

Morse[["cNMF_mOverlay"]]<-CreateAssayObject(usage_norm_f)


Morse@active.assay<-"cNMF_mOverlay"


Morse1$CellType_v4 <-factor(Morse1$CellType_v4 ,levels=c("DCs","cMonos","ncMonos","HSP macs","IM","SPP1 hi macs","CXCL10+ Macs","Mac UD (transitional mono-mac)","SPP1 mid macs","IntmAM 1","MT Mac","Mac UD (AM_UD1)","IntmAM 2","CD206hi FN1hi AM" ,"FABP4hi AM","IGF-1 AM"))

Morse1@active.assay<-"cNMF_mOverlay"
Idents(Morse1)<-Morse1$CellType_v4



Morse2<-Morse1[,Morse1$Disease_Status_1 %in% c("IPF_LL","IPF_UL")]
Morse2$Disease_Status_1<-factor(Morse2$Disease_Status_1,levels = c("IPF_UL","IPF_LL"))

features1<-c("Merged-cNMF-1","Merged-cNMF-2","Merged-cNMF-4","Merged-cNMF-5","Merged-cNMF-7","Merged-cNMF-6","Merged-cNMF-9","Merged-cNMF-10")

#Save Horizontal 8x 8
plot2 <- VlnPlot(Morse2, features1, stack = TRUE, sort = F, flip = TRUE,same.y.lims = T,split.by ="Disease_Status_1" , split.plot=T,cols=c("#C06CAB", "#D8A767")) +theme(axis.text.x = element_text(size = 15))
ggsave('VlnPlot_merged_cNMF_on_Morse_Disease_Split.pdf', dpi=900,width=8, height=8, units="in")




######## Fig 3E. Differential Abundance test in IPF LL vs UL (sccomp) #########


Morse$Disease_Status_1<-Morse$Disease_Status
Morse$Disease_Status_1<-ifelse(Morse$Disease_Status %in% c("IPF(Lower Lobe)","IPF(Lower Lobe)_Mac_Deplet"),"IPF_LL",Morse$Disease_Status_1)

Morse$Disease_Status_1<-ifelse(Morse$Disease_Status %in% c("IPF(Upper Lobe)"),"IPF_UL",Morse$Disease_Status_1)

Morse$Disease_Status_1<-gsub("[(]","_",Morse$Disease_Status_1)
Morse$Disease_Status_1<-gsub("[)]","",Morse$Disease_Status_1)
Morse$Disease_Status_1<-gsub(" ","",Morse$Disease_Status_1)


res = Morse |>
  sccomp_estimate(
    formula_composition =  ~ 0+Disease,
    .sample = Donor,
    .cell_group = CellType_v4,
    bimodal_mean_variability_association = TRUE,
    cores = 1
  )|>
  sccomp_test(  contrasts =  c("DiseaseIPF_LL - DiseaseIPF_UL"))




plots1 = res |> plot()


plots1$credible_intervals_1D + theme(axis.text.y = element_text( size = 10)) +theme_bw()
ggsave('Whisker Plot.pdf', dpi=900,width=7, height=5, units="in")


################ Fig 3F LRT test between HC vs IPF_UL vs IPF LL #########################
library(DESeq2)
library(EnhancedVolcano)
library(paletteer)


## LRT 
Morse_sub<-Morse[,!(Morse$Donor %in% c("SC31DNOR","SC31NOR","SC45NOR","SC14NOR"))]
table(Morse_sub$Disease_Status)

Morse_sub$Disease_1<-Morse_sub$Disease_Status

Morse_sub$Disease_1<-gsub("IPF[(]Lower Lobe[)]","IPF_L",Morse_sub$Disease_1)
Morse_sub$Disease_1<-gsub("IPF_L_Mac_Deplet","IPF_L",Morse_sub$Disease_1)
Morse_sub$Disease_1<-gsub("IPF[(]Upper Lobe[)]","IPF_U",Morse_sub$Disease_1)
Morse_sub$Disease_1<-ifelse(Morse_sub$Disease=="Normal","Normal",Morse_sub$Disease_1)


for (i in 1:length(Celltype)) {


subs<-Morse_sub[,Morse_sub$CellType_v2 == Celltype[i]]

mtx_aggr<-AggregateExpression(subs,assays = "RNA",slot="count",group.by = "Donor")
cluster_counts<-as.data.frame(mtx_aggr$RNA)


# Metadata generation

 metadata<-as.data.frame(subs@meta.data)

 meta_1<-metadata[,c("Donor","Disease_1")]

 meta_f<-unique(meta_1)

 rownames(meta_f)<-meta_f$Donor

meta_order<-meta_f[order(meta_f$Disease_1),]

 cluster_counts_f<-cluster_counts[,rownames(meta_order)]

meta_final<-meta_order
 # DESeq2 Run 

 meta_final$Sample_Name<-as.factor(meta_final$Donor)
 meta_final$Status<-as.factor(meta_final$Disease_1)

 dds <- DESeqDataSetFromMatrix(cluster_counts_f,
                               colData = meta_final,
                               design = ~ Status)


 #Prefilter
 dds <- dds[rowSums(counts(dds))> 3, ]
 rld <- rlog(dds, blind=TRUE)

  dds <- DESeq(dds, test="LRT", reduced=~1)
  res_LRT <- results(dds)

  head(res_LRT[order(res_LRT$padj),], 10)

  write.table(res_LRT,file =paste0(Celltype[15],"LRT_test_result.txt",sep="\t"))



  #Subset gene with significance
  sig_res_LRT <- res_LRT %>%
                 data.frame() %>%
                 tibble::rownames_to_column(var="gene") %>%
                 as_tibble() %>%
                 filter(padj < 0.05)

dim(sig_res_LRT)

  # Pull out sifnificant genes
  sigLRT_genes <- sig_res_LRT %>%
                  pull(gene)



  # Obtain normalized count values for those significant genes



 DEG_mat <- assay(rld)[ rownames(assay(rld))%in% sig_res_LRT$gene, colnames(assay(rld))%in% colnames(dds)]
 df <- data.frame(Group = SummarizedExperiment::colData(dds)[,c("Status")], row.names = rownames(SummarizedExperiment::colData(dds)))
 

 p<-pheatmap::pheatmap(DEG_mat, cluster_rows=TRUE, show_rownames=TRUE, cluster_cols=TRUE, annotation_col=df,scale="row", clustering_distance_rows = "correlation",clustering_distance_cols = "correlation", fontsize_row = 7.0)


  row_anno <-as.data.frame(cutree(p$tree_row, 4))
  colnames(row_anno) <- "Cluster"
  row_anno$Cluster <- as.factor(row_anno$Cluster)
  dev.off()

pal1<-paletteer_c("ggthemes::Orange-Blue Diverging", 100)

require(lattice)
pdf(paste0(Celltype[i],"_DEG_LRT_Hmap_0.05_1.pdf"),width=8,height=6)
 print(pheatmap::pheatmap(DEG_mat, cluster_rows=TRUE, color =rev(pal1),border_color = NA, show_rownames=TRUE, cluster_cols=TRUE, annotation_col=df,scale="row", clustering_distance_rows = "correlation",clustering_distance_cols = "correlation",annotation_row = row_anno,  fontsize_row = 1.0))
 dev.off()

 write.csv(row_anno,paste0(Celltype[i],"_Cluster_annotation.csv"))
}


########## Fig 5F. Pathway analysis ######################

# (Repeated for each clusters for Celltypes)

library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(ggnewscale)
library(fgsea)
library(dplyr)
library(msigdbr)

Gene_anno<-read.csv("Cluster_annotation.csv")

Gene_GO<-Gene_anno[Gene_anno$Cluster =="1",]
gene_up<-as.vector(Gene_GO$X)


########## GO  #########


ego_up<-enrichGO(gene_up,OrgDb = org.Hs.eg.db,keyType = "SYMBOL",ont="all",pAdjustMethod = "BH",pvalueCutoff = 0.05,qvalueCutoff = 0.05)

head(summary(ego_up))

ego_up

cluster_summary<-data.frame(ego_up)

write.table(cluster_summary,"Cluster1_GO.txt",sep="\t")

dotplot(ego_up, split="ONTOLOGY",label_format=100) + facet_grid(ONTOLOGY~., scale="free") + scale_fill_viridis(direction = -1)



ego_BP<-enrichGO(gene_up, OrgDb= org.Mm.eg.db,keyType = "SYMBOL",ont="BP",pAdjustMethod = "BH",pvalueCutoff = 0.05,qvalueCutoff = 0.05)
dotplot(ego_BP,showCategory=20,label_format=100)+ scale_fill_viridis(direction = -1)

ego2_BP <- simplify(ego_BP)
dotplot(ego2_BP,showCategory=20,label_format=100)+ scale_fill_viridis(direction = -1)

cnetplot(ego2_BP, colorEdge = TRUE)


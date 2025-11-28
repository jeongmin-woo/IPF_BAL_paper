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
library(scCustomize)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(ggnewscale)



######## Fig 1B. UMAP #########

merged<-readRDS("Myeloid_5DataMerge.rds")

print(head(merged@meta.data))
DefaultAssay(object = merged) <- "RNA"
Idents(merged)<-merged@meta.data$Celltype
DimPlot(merged, reduction = "umap", group.by="Celltype",label=T,repel = T)
ggsave('IPF_merged_UMAP.pdf', dpi=900,width=7, height=5, units="in")


######## Fig 1C. GEP Feature Plot #########

merged<-readRDS("Myeloid_5DataMerge.rds")

merged@active.assay<-"cNMF_K15"

FeaturePlot_scCustom(seurat_object = merged, features = rownames(usage_norm_f), num_columns = 4)

######## Fig 1C. Heatmap of top10 Contributing genes #########

usage <- read.table("BatchCorrected/Component_15/BatchCorrected.usages.k_15.dt_0_15.consensus.txt", sep='\t', row.names=1, header=TRUE)
spectra_score <- read.table("BatchCorrected/Component_15/BatchCorrected.gene_spectra_score.k_15.dt_0_15.txt", sep='\t', row.names=1, header=TRUE)
spectra_tpm <- read.table("BatchCorrected/Component_15/BatchCorrected.gene_spectra_tpm.k_15.dt_0_15.txt", sep='\t', row.names=1, header=TRUE)
head(usage)




get_top10_colnames <- function(row) {
  # Orders the values in descending order and gets the names of the top 20
  print(row[1:5])
  top_indices <- order(row, decreasing = TRUE)[1:10]
  return(colnames(spectra_score)[top_indices])
}

top_colnames <- apply(spectra_score, 1, get_top10_colnames)


top_colnames <- as.data.frame(top_colnames)

##### Fetch top10 genes from these GEP programs
gene_1<-as.vector(top_colnames[,1])
gene_2<-as.vector(top_colnames[,2])
gene_3<-as.vector(top_colnames[,4])
gene_4<-as.vector(top_colnames[,5])
gene_5<-as.vector(top_colnames[,7])
gene_6<-as.vector(top_colnames[,6])
gene_7<-as.vector(top_colnames[,9])
gene_8<-as.vector(top_colnames[,10])
gene_9<-as.vector(top_colnames[,3])
gene_10<-as.vector(top_colnames[,15])
gene_11<-as.vector(top_colnames[,14])
gene_12<-as.vector(top_colnames[,8])


genes_u<-Reduce(union, list(gene_1, gene_2, gene_3,gene_4,gene_5,gene_6,gene_7,gene_8,gene_9,gene_10,gene_11,gene_12))

merged1<-merged[,!(merged$Celltype=="Cycling Macs")]
Average_Mat<-AverageExpression(merged1,assays = "RNA",features=genes_u,group.by = "Celltype")
Average_Mac_scaled<- t(scale(t(Average_Mat$RNA), center = T, scale=T))

## Figure size 25 x 6 #
pal1<-paletteer_c("ggthemes::Orange-Blue Diverging", 100)
p1<-pheatmap::pheatmap(Average_Mac_scaled, cluster_cols=TRUE,cluster_rows=FALSE,color =rev(pal1), border_color = "black",fontsize_row=20.0,fontsize_col =20.0,gaps_row  = c(9,18,28,37,47,57,67,77,87,97,107))



######## Fig 1E. GO enrichment barplot  #########


get_top_colnames <- function(row) {
  # Orders the values in descending order and gets the names of the top 20
  print(row[1:5])
  top_indices <- order(row, decreasing = TRUE)[1:100]
  return(colnames(spectra_score)[top_indices])
}

top_colnames <- apply(spectra_score, 1, get_top_colnames)


top_colnames <- as.data.frame(top_colnames)



##### Fetch top100 genes from Usage 10 (repeated with other GEPs)
gene_up<-as.vector(top_colnames[,10])


ego_BP<-enrichGO(gene_up, OrgDb= org.Hs.eg.db,keyType = "SYMBOL",ont="BP",pAdjustMethod = "BH",pvalueCutoff = 0.05,qvalueCutoff = 0.05)
ego2_BP <- clusterProfiler::simplify(ego_BP)
p_bar<-mutate(ego2_BP, qscore = -log(p.adjust, base=10)) %>% 
  barplot(x="qscore",showCategory=10)

egomat<-mutate(ego2_BP, qscore = -log(p.adjust, base=10))
go_data<-as.data.frame(egomat)

#Top10 Bar 
go_data<-go_data %>% arrange(desc(qscore))
go_data<-go_data[1:10,]
go_data<-go_data %>% arrange(qscore)
go_data$Description <- factor(go_data$Description, levels = go_data$Description)

p_bar<-ggplot(data=go_data[1:10,], aes(x=Description, y=qscore,fill=Count)) +
  geom_bar(stat="identity", color="black") +scale_fill_continuous(low="#1C5F9E", high="#CE352E")+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 35)) +
  theme_bw()+coord_flip() +theme(axis.text.x = element_text(size = 16),axis.text.y  = element_text(size = 16),axis.title=element_text(size=16))
p_bar


require(lattice)
pdf("GEP10_GO_BP_BarPlot.pdf",width=10,height=6)
print(p_bar)
dev.off()

######## Fig 1F . UMAP eneration  #########

merged$Celltype-factor(merged$Celltype,levels = c("IGF1+ AM","FABP4hi AM","AM_UD","IntmAM","MT_AM","CD206hi FN1hi AM","ncMono","cMono","CXCL10+ Mono-Mac","MonoMac","IM","MonoDC","cDC2","cDC1","pDC","HSP+ Mac","SPP1+ Mac","ISGhi AM","Cycling Macs"))
DimPlot(merged,group.by="CellType_Integrated_v3",cols=c("#FEA316","#FFE353","#FFCC8A","#FFEFA3","#FFF2B6","#FFE6C4", "#1FAB1F", "#7B9A8D", "#84CB76","#A8DA9C" ,"#BAE2BA","#CBE9C3" ,"#DDF0D7","#49973B","#2E8B57" ,"#6D908B","#FF9896", "#FFAB32",  "#9467BD"),label.size = 5,label=T,repel=T)

                   
######## Fig 1G & H . Violin Plot of GEP usage score  #########


merged1<-merged[,!(merged$Celltype=="Cycling Macs")]
#Re-order based on clustering result
merged1$Celltype<-factor(merged1$Celltype,levels=c( "SPP1+ Mac", "IntmAM","CD206hi FN1hi AM","MonoMac", "AM_UD", "FABP4hi AM", "IGF1+ AM","MT_AM", "CXCL10+ Mono-Mac","ISGhi AM","pDC", "cDC1","cDC2","ncMono","cMono","MonoDC","IM","HSP+ Mac"  ))

merged1@active.assay<-"cNMF_K15"
Idents(merged1)<-merged1$Celltype



features<-c("Usage1","Usage2","Usage4","Usage5","Usage7","Usage6","Usage9","Usage10")

### Save Vertical 8.5 x 8
plot <- VlnPlot(merged1, features, stack = TRUE, sort = F, flip = TRUE,same.y.lims = T,y.max=1.0,cols=c("#272E6A","#272E6A","#272E6A","#272E6A","#272E6A","#BF0000","#BF0000","#BF0000","#89C75F","#89C75F","#89C75F","#89C75F")) +theme(axis.text.x = element_text(size = 15))


features1<-c("Usage1","Usage2","Usage4","Usage5","Usage7","Usage6","Usage9","Usage10")

plot2 <- VlnPlot(merged1, features1, stack = TRUE, sort = F, flip = TRUE,same.y.lims = T,split.by ="Disease_total" , split.plot=T,cols=c("#008494", "#F47D2B")) +theme(axis.text.x = element_text(size = 15))



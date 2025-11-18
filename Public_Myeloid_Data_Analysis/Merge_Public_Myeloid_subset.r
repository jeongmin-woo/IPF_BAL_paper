 library(UpSetR)
 library(ComplexHeatmap)


####### Generat UpSet Plot ################

 gene_adam<-rownames(Adams_sub@assays$RNA@counts)
 gene_HO<-rownames(BAL_sub@assays$RNA@counts)
 gene_Morse<-rownames(Morse_sub@assays$RNA@counts)
 gene_Habermann<-rownames(Haberman@assays$RNA@counts)
 gene_reyfman<-rownames(reyfman_sub@assays$RNA@counts)


lt<-list(Adam=gene_adam,Ho=gene_HO,Habermann=gene_Habermann,Reyfman=gene_reyfman,Morse=gene_Morse)
m = make_comb_mat(lt)


UpSet(m,comb_order = order(comb_size(m)),top_annotation = upset_top_annotation(m, add_numbers = TRUE))


########## only keep common genes from seurat object ###########
common_gene<-Reduce(intersect, lt)

Adams_c<-Adams_sub@assays$RNA@counts[rownames(Adams_sub@assays$RNA@counts) %in% common_gene,]
Adams_1<-CreateSeuratObject(Adams_c,project = "Adams",meta.data =Adams_sub@meta.data )


BAL_c<-BAL_sub@assays$RNA@counts[rownames(BAL_sub@assays$RNA@counts) %in% common_gene,]
HO_1<-CreateSeuratObject(BAL_c,project = "HO",meta.data =BAL_sub@meta.data )

Haberman_c<-Haberman@assays$RNA@counts[rownames(Haberman@assays$RNA@counts) %in% common_gene,]
Haberman_1<-CreateSeuratObject(Haberman_c,project = "Haberman",meta.data =Haberman@meta.data )

reyfman_c<-reyfman_sub@assays$RNA@counts[rownames(reyfman_sub@assays$RNA@counts) %in% common_gene,]
Reyfman_1<-CreateSeuratObject(reyfman_c,project = "Reyfman",meta.data =reyfman_sub@meta.data )

Morse_c<-Morse_sub@assays$RNA@counts[rownames(Morse_sub@assays$RNA@counts) %in% common_gene,]
Morse_1<-CreateSeuratObject(Morse_c,project = "Morse",meta.data =Morse_sub@meta.data )


 obj_list<-list(Reyfman_1,Haberman_1,HO_1,Morse_1)

merged <- merge(Adams_1, obj_list)
Idents(merged)<-merged$Data
merged_s<-subset(merged,downsample=15000)

merged<-merged_s



###prelim clustering########

merged@active.assay <- "RNA"
merged <- NormalizeData(merged)
merged <- FindVariableFeatures(merged)
merged <- ScaleData(merged)
merged <- RunPCA(object = merged, verbose = FALSE)
ElbowPlot(merged, ndims = 50)
merged <- RunUMAP(object = merged, dims = 1:20, verbose = FALSE)
merged <- FindNeighbors(object = merged, dims = 1:20, verbose = FALSE)
merged <- FindClusters(object = merged, verbose = FALSE, resolution = .7)
DimPlot(object = merged, label = TRUE, reduction="umap")

DimPlot(object = merged, group.by="Donor", reduction="umap")
DimPlot(object = merged, group.by="Data", reduction="umap")


merged <- harmony::RunHarmony(merged, group.by.vars = c("Donor","Data"),lambda=NULL)
merged <- RunUMAP(object = merged, dims = 1:20, verbose = FALSE, reduction = "harmony")

merged <- FindNeighbors(object = merged, dims = 1:20, verbose = FALSE,reduction = "harmony")
merged <- FindClusters(object = merged, verbose = FALSE, resolution = .7)
DimPlot(object = merged, label = TRUE, reduction="umap")

DimPlot(object = merged, label = TRUE, reduction="umap",group.by="Data")

merged$CellType_Adams<-ifelse(merged$Data=="Adams",merged$CellType_v2,NA)
merged$CellType_Habermann<-ifelse(merged$Data=="Haberman",merged$CellType_v2,NA)
merged$CellType_HO<-ifelse(merged$Data=="HO",merged$CellType_v2,NA)
merged$CellType_Morse<-ifelse(merged$Data=="Morse",merged$CellType_v2,NA)
merged$CellType_Reyfman<-ifelse(merged$Data=="Reyfman",merged$CellType_v2,NA)

DimPlot(object = merged, reduction="umap",group.by = "Celltype_origin",split.by = "Data",ncol = 2)

DimPlot(object = merged, reduction="umap",group.by = "CellType_v2",split.by = "Data",label=T,repel=T,ncol = 2)

DimPlot(object = merged, reduction="umap",group.by = "CellType_Adams",split.by = "Data",label=T,repel=T,ncol = 2)

ggplot(merged@meta.data, aes(RNA_snn_res.0.7, fill=CellType_v2)) + geom_bar(position="fill", colour="black") +theme_classic()+ labs(y="Proportion") +  theme(axis.text.x = element_text(angle = 45, hjust=1,size=11))
ggplot(merged@meta.data, aes(RNA_snn_res.0.7, fill=CellType_Adams)) + geom_bar(position="fill", colour="black") +theme_classic()+ labs(y="Proportion") +  theme(axis.text.x = element_text( hjust=1,size=11))
ggplot(merged@meta.data, aes(RNA_snn_res.0.7, fill=CellType_Habermann)) + geom_bar(position="fill", colour="black") +theme_classic()+ labs(y="Proportion") +  theme(axis.text.x = element_text( size=11))
ggplot(merged@meta.data, aes(RNA_snn_res.0.7, fill=CellType_HO)) + geom_bar(position="fill", colour="black") +theme_classic()+ labs(y="Proportion") +  theme(axis.text.x = element_text( size=11))
ggplot(merged@meta.data, aes(RNA_snn_res.0.7, fill=CellType_Morse)) + geom_bar(position="fill", colour="black") +theme_classic()+ labs(y="Proportion") +  theme(axis.text.x = element_text( size=11))
ggplot(merged@meta.data, aes(RNA_snn_res.0.7, fill=CellType_Reyfman)) + geom_bar(position="fill", colour="black") +theme_classic()+ labs(y="Proportion") +  theme(axis.text.x = element_text( size=11))
ggplot(merged@meta.data, aes(RNA_snn_res.0.7, fill=Data)) + geom_bar(position="fill", colour="black") +theme_classic()+ labs(y="Proportion") +  theme(axis.text.x = element_text( size=11))


###### After annotating merged Cluster ####

### Frequency_BarPlot 

ggplot(merged@meta.data,aes(x=CellType_Integrated,fill=Data)) + geom_bar(color="black") +theme_classic()+ labs(y="Count") +  theme(axis.text.x = element_text(angle = 45, hjust=1,size=11))

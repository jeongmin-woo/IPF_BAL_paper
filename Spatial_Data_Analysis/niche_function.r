library(Seurat )
library(pheatmap)
library(ggplot2)
library(dplyr)
library(Matrix)


#' Creates a matrix of expression values for each cell's nearest spatial neighbours
#' To be used as an alternative niche assay without needing to specify cell types apriori - can test for enrichment afterwards
#' So this is a more unsupervised approach
#'@param object seurat object
#'@param fov name of the image/fov to compute neighbours from
#'@param neighbors.k nearest neighbours to use
#'@param assay.counts which spatial assay to pull counts from
#'@return matrix of nearest neighbour aggregated counts for each feature, excluding the cell itself
#'
BuildNicheExpressionAssay.PerImage <- function(
    object,
    fov,
    neighbors.k = 30,
    assay.counts = "XENIUM"
) {
  # find neighbors based on tissue position
  coords <- GetTissueCoordinates(object[[fov]], which = "centroids")
  cells <- coords$cell
  rownames(coords) <- cells
  coords <- as.matrix(coords[ , c("x", "y")])
  neighbors <- FindNeighbors(coords, k.param = neighbors.k)
  neighbors$nn <- neighbors$nn[cells, cells]

  diag(neighbors$nn) <- 0 # dont count transcriptome of the cell itself, just neighbours?

  mt <- object[[assay.counts]]@counts[, cells]

  sum.mtx <- as.matrix(neighbors$nn %*% t(mt))

  return(sum.mtx)
}

#' Creates a matrix of expression values for each cell's nearest spatial neighbours
#' To be used as an alternative niche assay without needing to specify cell types apriori - can test for enrichment afterwards
#' So this is a more unsupervised approach
#' This function will just apply 'BuildNicheExpressionAssay.PerImage' to the whole seurat object of multiple merged slides
#'@param object seurat object
#'@param neighbors.k nearest neighbours to use
#'@param assay.counts which spatial assay to pull counts from
#'@return matrix of nearest neighbour aggregated counts for each feature, excluding the cell itself
#'

BuildNicheExpressionAssay <- function(object,
                                      neighbors.k = 30,
                                      assay.counts = "XENIUM"){

  res <- lapply(Images(object), function(x){
    print("Calculating FOV:")
    print(x)
    results <- BuildNicheExpressionAssay.PerImage(object, fov=x,  neighbors.k = neighbors.k, assay.counts = assay.counts)

  })

  res <- do.call(rbind, res)
  res
}

#' Given a vector of niches and a corresponding vector of cell types, will calculate cell type enrichment table in each niche
#' to be used with downstream visualisation functions
#' @niches - vector of niche names per cell
#' @cell.types - vector of cell types per cell, in the same cell order as that of niches
#' @return data frame with enrichment estimate and p value for each  cell type - cell niche combination
CellTypeNicheEnrichment <- function(niches, cell.types){

  df <- list()

  for ( niche in unique(niches)){

    for( cell in unique(cell.types)){

      res <- fisher.test(x=niches == niche, y=cell.types == cell )

      df[[paste0(cell, niche)]] <- c(res$estimate, res$p.value, cell, niche)
    }
  }

  df <- as.data.frame(do.call(rbind, df))
  colnames(df) <- c("Estimate", "p.val", "Cell Type", "Niche")
  df$p.val <- as.numeric(df$p.val)
  df$FDR <- p.adjust(df$p.val)
  df
}

#' given niche results from CellTypeNicheEnrichment, will plot a cell type enrichment barplot for any query niche
#' @niche.results - output table from CellTypeNicheEnrichment
#' @niche - which niche to plot, must be present in the results
#' @return barplot with enrichment results
PlotCellTypeNicheEnrichmentBarPlot <- function(niche.results, niche=0){

  df <- niche.results[niche.results$Niche == niche,]
  df <- df[order(as.numeric(df$Estimate)), ]
  df$`Cell Type` <- factor(df$`Cell Type`, levels=df$`Cell Type`)
  ggplot(df, aes(`Cell Type`, log2(as.numeric(Estimate)),
                 fill=-log10(as.numeric(FDR) + 10^-300))
  ) + geom_bar(stat="identity", colour="black"
  )  + coord_flip() + scale_fill_viridis_c() + theme_classic(base_size = 16) + labs(
    x="Cell Type", y="Relative Depletion <-----> Relative Enrichment", fill="-log10 FDR"
  ) + geom_hline(yintercept = c(-1, 1), lty=2, colour="red")
}


#'Given niche results from cellTypeNicheEnrichment, will plot an overview heatmap of cell type to cell niche enrichment.
#' @niche.results  - output table from CellTypeNicheEnrichment
#' @return heatmap with enrichment results
PlotCellTypeNicheEnrichmentHeatmap <- function(niche.results){

  niche.results$Estimate <- as.numeric(niche.results$Estimate)
  mat <- log2(reshape2::acast(niche.results, `Cell Type`~ Niche, value.var = "Estimate")+ 0.01)
  niche.results$STAR <- "n.s"; niche.results$STAR[niche.results$FDR < 0.05] <- "*"; niche.results$STAR[niche.results$FDR < 0.01] <- "**"; niche.results$STAR[niche.results$FDR < 0.001] <- "***"
  mat.labels <-reshape2::acast(niche.results, `Cell Type`~ Niche, value.var = "STAR")
  mat[mat < -4] <- -4; mat[ mat > 4] <- 4
  pheatmap::pheatmap(mat, border_color = "black", display_numbers = mat.labels, fontsize_number = 14,
                     color=Seurat:::SpatialColors(n=100), number_color = "black" )

}

#' Creates an aggregated expression assay for cell spatial neighbours, but only if they are of a particular cell type/group etc.
#' This is for downstream testing, e.g. whether a T-cell in a neighbourhood of an epithelial cell expresses anything different to when it is in a neighbourhood of say a b-cell
#' @object seurat object
#' @fov field of view/image for which to calculate spatial neighbours
#' @group.by meta data column containing cell type grouping
#' @idents which cell identities in the cell type column to consider
#' @neighbours.k how many nearest neighbours to consider
#' @assay.counts which spatial assay to use
#' @return matrix of gene expression counts for cell spatial neighbours of the specified cell type(s)
#'
BuildCellSpecificNicheExpressionAssay.PerImage <- function(
    object,
    fov,
    group.by,
    idents=c(14, 0, 6),
    neighbors.k = 30,
    assay.counts = "XENIUM"
) {

  #print(object)
  # find neighbors based on tissue position
  coords <- GetTissueCoordinates(object[[fov]], which = "centroids")
  cells <- coords$cell
  rownames(coords) <- cells
  coords <- as.matrix(coords[ , c("x", "y")])
  neighbors <- FindNeighbors(coords, k.param = neighbors.k)
  neighbors$nn <- neighbors$nn[cells, cells]

  diag(neighbors$nn) <- 0 # dont count transcriptome of the cell itself, just neighbours

  mt <- object[[assay.counts]]@counts[, cells]
  mt[, !(object[, cells][[group.by]][, 1] %in% idents)] <- 0 # zero counts from cells that are not selected

  sum.mtx <- as.matrix(neighbors$nn %*% t(mt))

  return(sum.mtx)
}


#' Will find niche clusters given an arbitraty niche matrix - can be expression or cell type counts matrix
#' @object seurat object
#' @niche.matrix any niche matrix that matrches the cells of the seurat
#' @resolution - clustering resolution, this generally should be lower than for single cell
#' @dims - how many dims to use, again normally 10 is good, too many gets weird small clusters
#' @niche.assay.name name of the new assay where teh data will be stored in the seurat
#' @harmonise whether to harmonise the niche matrix between any variables -e.g. different slides
#' @harmony.group.by if harmonise, provide group.by arguement
#' @niche.name what name to give to the niche groups in the object meta data
#'
FindNiches <- function(object, niche.matrix, resolution=0.1, dims=1:10, niche.assay.name="niche",
                       harmonise=F, harmony.group.by=NULL, niche.name="TissueNiche"){


  object[[niche.assay.name]] <- CreateAssayObject(niche.matrix[, Cells(object)])
  object@active.assay <- niche.assay.name
  object <- SCTransform(object, assay = niche.assay.name, clip.range = c(-10, 10))
  object <- RunPCA(object, npcs = 30)
  tmp <- object@active.ident

  if( harmonise){

    object <- harmony::RunHarmony(object, group.by.vars=harmony.group.by, assay=niche.assay.name)
    object <- FindNeighbors(object, reduction = "harmony", dims = dims)
    object <- FindClusters(object, resolution = resolution)
  }
  else{

    object <- FindNeighbors(object, reduction = "pca", dims = dims)
    object <- FindClusters(object, resolution = resolution)
  }

  object <- StashIdent(object, niche.name) # stack the niche cluster identity in the meta data
  object@active.ident <- tmp ## reset back to original active ident so it is not accidentally lost, like I am likely to do
  object
}


FindNiches.NoDimRed <- function(object, niche.matrix, k.centers=15, niche.assay.name="niche",
                        niche.name="TissueNiche"){


  object[[niche.assay.name]] <- CreateAssayObject(niche.matrix[, Cells(object)])
  object@active.assay <- niche.assay.name
  object <- ScaleData(object)

  results <- kmeans(
    x = t(object[[niche.assay.name]]@scale.data),
    centers = k.centers,
    nstart = 30
  )
  object@meta.data[, niche.name] <- results[["cluster"]]
  object
}

#' Preprocess the single-cell raw data using functions in the \code{Seurat} package
#'
#' #' This function provide a simplified-version of Seurat analysis pipeline for single-cell RNA-seq data. It contains the following steps in the pipeline:
#' \itemize{
#'    \item Normalize the count data present in a given assay.
#'    \item Identify the variable features.
#'    \item Scales and centers features in the dataset.
#'    \item Run a PCA dimensionality reduction.
#'    \item Constructs a Shared Nearest Neighbor (SNN) Graph for a given dataset.
#'    \item Runs the Uniform Manifold Approximation and Projection (UMAP) dimensional reduction technique.
#' }
#'
#' @title run seurat preprocessing function
#'
#' @description Single-cell preprocess
#'
#' @details you can use this function to caculate x+1,then return the value of x+1.
#'
#' @param Obejct Obejct is a Seurat object
#' @param normalization.method Method for normalization.
#'   \itemize{
#'   \item LogNormalize: Feature counts for each cell are divided by the total counts for that cell and multiplied by the scale.factor.
#'   This is then natural-log transformed using log1p.
#'   \item CLR: Applies a centered log ratio transformation.
#'   \item RC: Relative counts. Feature counts for each cell are divided by the total counts for that cell and multiplied by the scale.factor.
#'   No log-transformation is applied. For counts per million (CPM) set \code{scale.factor = 1e6}.
#' }
#' @param scale.factor Sets the scale factor for cell-level normalization.
#' @param selection.method How to choose top variable features. Choose one of :
#'   \itemize{
#'   \item vst: First, fits a line to the relationship of log(variance) and log(mean) using local polynomial regression (loess).
#'   Then standardizes the feature values using the observed mean and expected variance (given by the fitted line).
#'   Feature variance is then calculated on the standardized values after clipping to a maximum (see clip.max parameter).
#'   \item mean.var.plot (mvp): First, uses a function to calculate average expression (mean.function) and dispersion (dispersion.function)
#'   for each feature. Next, divides features into num.bin (deafult 20) bins based on their average expression, and calculates
#'   z-scores for dispersion within each bin. The purpose of this is to identify variable features while controlling for the strong
#'   relationship between variability and average expression.
#'   \item dispersion (disp): selects the genes with the highest dispersion values
#'   }
#' @param dims_Neighbors Dimensions of reduction to use as input.
#' @param dims_UMAP Which dimensions to use as input features for UMAP.
#' @param verbose Print output.
#'
#' @return a Seurat object
#' @import Seurat
#' @export


### Single-cell preprocess
run_seurat <- function(Obejct,
                       normalization.method = "LogNormalize",
                       scale.factor = 10000,
                       selection.method = "vst",
                       dims_Neighbors = 1:40,
                       dims_UMAP = 1:10,
                       verbose = TRUE){
  
  Obejct <- Seurat::NormalizeData(object = Obejct, normalization.method = normalization.method, scale.factor = scale.factor, verbose = verbose)
  Obejct <- Seurat::FindVariableFeatures(object = Obejct, nfeatures = 3000,selection.method = selection.method, verbose = verbose)
  Obejct <- Seurat::ScaleData(object = Obejct, verbose = verbose)
  Obejct <- Seurat::RunPCA(object = Obejct, features = VariableFeatures(Obejct), verbose = verbose)
  Obejct <- Seurat::FindNeighbors(object = Obejct, dims = dims_Neighbors, verbose = verbose)
  Obejct <- Seurat::RunUMAP(object = Obejct, dims = dims_UMAP, verbose = verbose)
}



#' The scAB_data Class
#'
#' The scAB_data object is created from a single-cell RNA-seq data, bulk RNA-seq data, and phenotype data.
#' Its input consists of a Seurat object, a digital matrix, and a matrix or vector providing phenotypic information. Among them, the Seurat object should contain Shared Nearest Neighbor (SNN) Graph information.
#' The class provides functions for data preprocessing, integrative analysis, and visualization.
#'
#' The key slots used in the scAB_data object are described below.
#'
#' @slot X Correlation matrix between bulk data and single-cell data (individuals should be in rows and cells in columns)
#' @slot S Diagonal matrix for information of individuals phenotype
#' @slot A Shared Nearest Neighbor (SNN) Graph matrix
#' @slot D Degree Matrix of Shared Neighbor Graph
#' @slot L Laplacian Matrix
#' @slot phenotype Matrix of phenotype information
#'
#' @exportClass scAB_data
#' @useDynLib scAB

setClass("scAB_data",slots=list(X="matrix",
                                S="matrix",
                                L="matrix",
                                D="matrix",
                                A="matrix",
                                phenotype="matrix",
                                mathod="character")
)


###  scAB_data preprocess
#' preprocess the single-cell data, bulk data, and phenotype data.
#'
#' @param Obejct  Seurat object
#' @param bulk_dataset  matrix of bulk data
#' @param phenotype Phenotype data, a matrix with two columns "time" and "state", or a vector
#' @param method method "survival" or "binary"
#'
#' @return a scAB_data
#' @importFrom preprocessCore normalize.quantiles
#' @export
#'
#' @examples
create_scAB <- function(Obejct,bulk_dataset,phenotype,method=c("survival","binary")){
  # cell neighbors
  method=match.arg(method)
  if(is.null(Obejct@ graphs$ RNA_snn)) print("'RNA_snn' is not found, please run FindNeighbors function in Seurat.")
  A <- as.matrix(Obejct@ graphs$ RNA_snn)
  diag(A) <- 0
  A[which(A != 0)] <- 1
  D <- diag(rowSums(A))
  D12 <- diag(1/sqrt(rowSums(A)))
  L <- D12%*%(D-A)%*%D12
  Dhat <- D12%*%(D)%*%D12
  Ahat <- D12%*%(A)%*%D12
  
  # similarity matrix
  sc_exprs <- as.data.frame(Obejct@ assays$ RNA@data)
  common <- intersect(rownames(bulk_dataset), rownames(sc_exprs))
  dataset0 <- cbind(bulk_dataset[common,], sc_exprs[common,])         # Dataset before quantile normalization.
  dataset1 <- preprocessCore::normalize.quantiles(as.matrix(dataset0))                           # Dataset after  quantile normalization.
  rownames(dataset1) <- rownames(dataset0)
  colnames(dataset1) <- colnames(dataset0)
  Expression_bulk <- dataset1[,1:ncol(bulk_dataset)]
  Expression_cell <- dataset1[,(ncol(bulk_dataset) + 1):ncol(dataset1)]
  X <- cor(Expression_bulk, Expression_cell)
  X=X/norm(X,"F")
  
  # phenotype ranking
  if(method=="survival"){
    ss <- guanrank(phenotype[,c("time","status")])
    S <- diag(1-ss[rownames(phenotype),3])
  }
  else{
    S <- diag(1-phenotype)
  }
  
  # return
  obj <- list(X=X,S=S,L=L,D=Dhat,A=Ahat,phenotype=phenotype,method=method)
  class(obj) <- "scAB_data"
  return(obj)
}

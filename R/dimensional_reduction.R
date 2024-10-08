#' RUNPCA
#'
#' @description
#' Perform principal component analysis on the original gene expression matrix
#'
#' @param object MGPfact object
#' @param pc_num Number of principal components
#' @param center a logical value indicating whether the variables should be shifted to be zero centered.
#' Alternately, a vector of length equal the number of columns of x can be supplied. The value is passed to scale.
#' @param scale a logical value indicating whether the variables should be scaled to have unit variance before the analysis takes place.
#'
#' @export
#'
RUNPCA <- function(object,
                   pc_num = 30,
                   center = TRUE,
                   scale = FALSE){

  s = object@assay$data_matrix
  s = as.matrix(s)

  # slower
  # object@DimReducs@PCA = prcomp(s, center = center, scale = scale)

  # faster
  s2 = scale(s, center = center, scale = scale)
  pca.results <- irlba(A = s2, nv = pc_num)
  cell.embeddings <- pca.results$u %*% diag(pca.results$d) # 常规操作
  # cell.embeddings <- pca.results$u
  feature.loadings <- pca.results$v
  sdev <- pca.results$d/sqrt(max(1, ncol(s) - 1))

  command = GetCommand()
  object@Command$DimReducs$PCA = command
  object@DimReducs@PCA = list(cell.embeddings = cell.embeddings,
                              feature.loadings = feature.loadings,
                              sdev = sdev)
  return(object)
}

#' RUNtSNE
#'
#' @description
#' Perform tSNE dimensionality reduction based on the principal component analysis
#'
#' @param object MGPfact object
#' @param npcs which pc
#' @export
#'
RUNtSNE <- function(object,
                    npcs = 1:30){

  pca = object@DimReducs@PCA$cell.embeddings
  tsne = Rtsne(pca[,npcs], pca = FALSE)
  command = GetCommand()
  object@Command$DimReducs$tSNE = command
  object@DimReducs@tSNE = list(cell.embeddings = tsne$Y)
  return(object)
}

#' RUNUMAP
#'
#' @description
#' Perform UMAP dimensionality reduction based on the principal component analysis
#'
#' @param object MGPfact object
#' @param npcs which pc
#' @export
#'
RUNUMAP <- function(object,
                    npcs = 1:30){
  x = object@DimReducs@PCA$cell.embeddings[,npcs]
  umap = umap(x)
  command = GetCommand()
  object@Command$DimReducs$UMAP = command
  object@DimReducs@UMAP = list(cell.embeddings = umap$layout)
  return(object)
}

#' RUNDM
#'
#' @description
#' Perform diffusion map dimensionality reduction based on the expression matrix
#'
#' @param object MGPfact object
#' @param neigen number of components
#' @param method the method of compute dist matrix
#' @param cores the number of threads that can be used to calculate distances
#'
#' @import parallelDist philentropy diffusionMap
#'
#' @export
#'
# RUNDM <- function(object,
#                   neigen = 10,
#                   method = "euclidean",
#                   cores = 1
#                   ){
#   s = object@assay$data_matrix
#   s = as.matrix(s)
#   # mat = distance(s, method = method)
#   mat = parDist(s, method = method, diag = TRUE, upper = TRUE, threads = cores)
#   mat = as.matrix(mat)
#
#   dm = diffuse(D = mat, neigen = neigen)
#   object@DimReducs@DM = dm
#   command = GetCommand()
#   object@Command$DimReducs$DM = command
#   return(object)
# }

#' GetNPC
#' @description
#' Get the number of principal components, default is sd > 1
#'
#' @param pca_result pca result from prcomp
#' @param sd_cutoff Cut-off threshold for principal component variance
#'
#' @export
#'
GetNPC <- function(pca_result = NULL, sd_cutoff = 1){

  # pca = MGPfact_object@DimReduc$PCA$pca_result
  pca = pca_result
  npc = length(which(pca$sdev >= sd_cutoff))
  npc = ifelse(npc == 0, npc+1, npc)
  npc = ifelse(npc == 1, npc+1, npc)
  npc = ifelse(npc > 10, 10, npc)

  return(npc)
}
#' Marchenko-Pastur Significant PCs (urd)
#' @description
#' The Marchenko Pastur Law (MP) predicts the theoretical upper and lower bounds
#' on the null distribution of eigenvalues for an MxN random matrix. We take
#' significant principal components (PCs) as those with eigenvalues greater than
#' the maximum eigenvalue predicted for random data. This function assumes that
#' the data has mean 0 and variance 1 (i.e. that the data has been centered and
#' scaled).
#'
#' @param M (Numeric) Number of rows in input data
#' @param N (Numeric) Number of columns in input data
#' @param pca.sdev (Numeric vector) Standard deviations for each principal component
#' @param factor (Numeric) Factor to multiply eigenvalue null upper bound before determining significance.
#' @param do.print (Logical) Whether to report the results
#' @return Logical vector of whether each PC is significant.
#'
#' @export
#'
pcaMarchenkoPastur <- function(M, N, pca.sdev, factor = 2, do.print=T) {
  pca.eigenvalue <- (pca.sdev)^2
  marchenko.pastur.max <- (1+sqrt(M/N))^2
  pca.sig <- pca.eigenvalue > (marchenko.pastur.max * factor)
  if (do.print) {
    print(paste("Marchenko-Pastur eigenvalue null upper bound:", marchenko.pastur.max))
    if (factor != 1) {
      print(paste(length(which(pca.sig)), "PCs have eigenvalues larger than", factor, "times null upper bound."))
    } else {
      print(paste(length(which(pca.eigenvalue > marchenko.pastur.max)), "PCs have larger eigenvalues."))
    }
  }
  pca.sig = which(pca.sig)
  return(pca.sig)
}

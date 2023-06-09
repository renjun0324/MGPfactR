#' RUNPCA
#'
#' @description
#' Perform principal component analysis on the original gene expression matrix
#'
#' @param object celltrek object
#' @param center a logical value indicating whether the variables should be shifted to be zero centered.
#' Alternately, a vector of length equal the number of columns of x can be supplied. The value is passed to scale.
#' @param scale a logical value indicating whether the variables should be scaled to have unit variance before the analysis takes place.
#'
#' @export
#'
RUNPCA <- function(object,
                   center = TRUE,
                   scale = FALSE){

  s = object@assay$data_matrix
  s = as.matrix(s)
  object@DimReducs@PCA = prcomp(s, center = center, scale = scale)
  command = GetCommand()
  object@Command$DimReducs$PCA = command

  return(object)
}

#' RUNtSNE
#'
#' @description
#' Perform tSNE dimensionality reduction based on the principal component analysis
#'
#' @param object celltrek object
#' @param npc which pc
#' @export
#'
RUNtSNE <- function(object,
                    npc = 1:30){
  require(Rtsne)

  pca = object@DimReducs@PCA$x
  object@DimReducs@tSNE = Rtsne(pca[,npc], pca = FALSE)
  command = GetCommand()
  object@Command$DimReducs$tSNE = command
  return(object)
}

#' RUNUMAP
#'
#' @description
#' Perform UMAP dimensionality reduction based on the principal component analysis
#'
#' @param object celltrek object
#' @param npc which pc
#' @export
#'
RUNUMAP <- function(object,
                    npc = 1:30){
  require(umap)

  pca = object@DimReducs@PCA$x
  object@DimReducs@UMAP = umap::umap(pca[,npc])
  command = GetCommand()
  object@Command$DimReducs$UMAP = command
  return(object)
}

#' RUNDM
#'
#' @description
#' Perform diffusion map dimensionality reduction based on the expression matrix
#'
#' @param object celltrek object
#' @param neigen number of components
#' @param method the method of compute dist matrix
#' @param cores the number of threads that can be used to calculate distances
#'
#' @export
#'
RUNDM <- function(object,
                  neigen = 10,
                  method = "euclidean",
                  cores = 1
                  ){

  require(parallelDist)
  require(philentropy)
  require(diffusionMap)

  s = object@assay$data_matrix
  s = as.matrix(s)
  # mat = distance(s, method = method)
  mat = parDist(s, method = method, diag = TRUE, upper = TRUE, threads = cores)
  mat = as.matrix(mat)

  dm = diffuse(D = mat, neigen = neigen)
  object@DimReducs@DM = dm
  command = GetCommand()
  object@Command$DimReducs$DM = command
  return(object)
}

#' GetNPC
#'
#' Get the number of principal components, default is sd > 1
#'
#' @param pca_result pca result from prcomp
#' @param sd_cutoff Cut-off threshold for principal component variance
#'
#' @export
#'
GetNPC <- function(pca_result = NULL, sd_cutoff = 1){

  # pca = celltrek_object@DimReduc$PCA$pca_result
  pca = pca_result
  npc = length(which(pca$sdev >= sd_cutoff))
  npc = ifelse(npc == 0, npc+1, npc)
  npc = ifelse(npc == 1, npc+1, npc)
  npc = ifelse(npc > 10, 10, npc)

  return(npc)
}

#' Marchenko-Pastur Significant PCs (urd)
#'
#' The Marchenko Pastur Law (MP) predicts the theoretical upper and lower bounds
#' on the null distribution of eigenvalues for an MxN random matrix. We take
#' significant principal components (PCs) as those with eigenvalues greater than
#' the maximum eigenvalue predicted for random data. This function assumes that
#' the data has mean 0 and variance 1 (i.e. that the data has been centered and
#' scaled). This called automatically by \code{\link{calcPCA}} and the results
#' are stored in slot \code{pca.sig}.
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
  return(pca.sig)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                              Class definitions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' The DimReduc-class
#'
#' @slot PCA pca result
#' @slot tSNE tsne result
#' @slot DM diffusion map result
#' @slot UMAP umap result
#'
#' @name DimReduc-class
#' @rdname DimReduc-class
#' @exportClass DimReduc
#'·
DimReducs <- setClass(
  Class = 'DimReducs',
  slots = c(
    PCA = 'ANY',
    tSNE = 'ANY',
    DM = 'ANY',
    UMAP = 'ANY'
  )
)

#' The MGPfact-class
#'
#' @slot assay expression matrix
#' @slot OptimResult julia results
#' @slot MURP murp results
#' @slot DimReducs dimensionality reduction results
#' @slot MetaData meta data
#' @slot Date object creation date
#' @slot BeautGene results related to genes
#' @slot Tree trajectory tree
#' @slot GPR gaussian process results
#' @slot Command execution log
#'
#' @name MGPfact-class
#' @rdname MGPfact-class
#' @exportClass MGPfact
#'
MGPfact <- setClass(
  Class = 'MGPfact',
  slots = c(
    assay = 'list',
    OptimResult = 'ANY',
    MURP = 'ANY',
    DimReducs = 'DimReducs',
    MetaData = 'data.frame',
    Settings = 'ANY',
    Date = 'ANY',
    BeautGene = 'list',
    Tree = 'list',
    GPR = 'list',
    Command = 'list'
  )
)

#' The PseOptimResult-class
#'
#' @slot sample_value sampling results for each iteration
#' @slot dimnames names of sampling value
#' @slot burnin burning period
#' @slot chains chain number
#' @slot iter iteration number
#' @slot hasinits is there an initial value
#' @slot inits initial value
#' @slot iter_chain_cor correlation between pseudotime and
#' sampling time points across different iterations and Markov chains.
#' @slot t_pred pseudotime of murp
#' @slot sortindex_cor (ready to delete)
#' @slot diagnostic convergence judgment results
#' @slot logpdf The objective function value of each optimization
#' @name PseOptimResult-class
#' @rdname PseOptimResult-class
#' @exportClass PseOptimResult
#'
PseOptimResult <- setClass(
  Class = 'PseOptimResult',
  slots = c(
    sample_value = 'array',
    dimnames = 'list',
    burnin = 'numeric',
    chains = 'vector',
    iter = 'numeric',
    hasinits = 'logical',
    inits = 'list',
    iter_chain_cor = 'list',
    t_pred = 'list',
    sortindex_cor = 'list',
    diagnostic = 'list',
    logpdf = 'list'
  )
)

#' The TrackOptimResult-class
#'
#' @slot sample_value sampling results for each iteration
#' @slot dimnames names of sampling value
#' @slot burnin burning period
#' @slot chains chain number
#' @slot iter iteration number
#' @slot hasinits is there an initial value
#' @slot inits inits
#' @name TrackOptimResult-class
#' @rdname TrackOptimResult-class
#' @exportClass TrackOptimResult
TrackOptimResult <- setClass(
  Class = 'TrackOptimResult',
  slots = c(
    sample_value = 'array',
    dimnames = 'list',
    burnin = 'numeric',
    chains = 'vector',
    iter = 'numeric',
    hasinits = 'logical',
    inits = 'list',
    diagnostic = 'list',
    logpdf = 'list',
    gene_weight = 'list'
  )
)

#' The Settings-class
#'
#' @slot sample_value sampling results for each iteration
#' @slot dimnames names of sampling value
#' @slot burnin burning period
#' @slot chains chain number
#' @slot iter iteration number
#' @slot hasinits is there an initial value
#' @slot inits inits
#' @name Settings-class
#' @rdname Settings-class
#' @exportClass Settings
Settings <- setClass(
  Class = 'Settings',
  slots = c(
    datasetTitle = "ANY",
    settings="list"
  )
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                              Create function
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Create PseOptimResult Object
#'
#' @param sample_value sampling results for each iteration
#' @param dimnames names of sampling value
#' @param burnin burning period
#' @param chains chain number,vector
#' @param iter iteration number
#' @param hasinits is there an initial value
#' @param inits inital value
#' @param iter_chain_cor Correlation of pseudo-time results between different Markov chains
#' @param sortindex_cor (ready to delete)
#' @param t_pred pseudotime list
#' @param diagnostic convergence judgment results
#' @param logpdf The objective function value of each optimization
#'
#' @export
#'
CreatePseOptimResultObject <- function(sample_value = NULL,
                                       dimnames = NULL,
                                       burnin = NULL,
                                       chains = NULL,
                                       iter = NULL,
                                       hasinits = FALSE,
                                       inits = list(),
                                       iter_chain_cor = list(),
                                       sortindex_cor = list(),
                                       t_pred = list(),
                                       diagnostic = list(),
                                       logpdf = list()

){

  # check sample_value and dimnames
  if(missing(x = sample_value) && missing(x = dimnames)){
    stop("Must provide either 'sample_value' or 'dimnames'")
  }

  # check hasinits and inits
  if(hasinits){
    if(missing(inits)){
      stop("has inits? if has, please provide it in 'inits = '")
    }
    if(!is.list(inits)){
      stop("inits must be list")
    }
  }

  # check inits
  object.list <- list(
    'Class' = 'PseOptimResult',
    'sample_value' = sample_value,
    'dimnames' = dimnames,
    'chains' = chains,
    'burnin' = burnin,
    'iter' = iter,
    'hasinits' = hasinits,
    'inits' = inits,
    'iter_chain_cor' = iter_chain_cor,
    't_pred' = t_pred,
    'sortindex_cor' = sortindex_cor,
    'diagnostic' = diagnostic,
    'logpdf' = logpdf
  )

  object <- do.call(what = 'new', args = object.list)
  return(object)
}

#' Create TrackOptimResult Object
#' @description
#' create track julia result object
#' @param sample_value sampling results for each iteration
#' @param dimnames names of sampling value
#' @param burnin burning period
#' @param chains chain number,vector
#' @param iter iteration number
#' @param hasinits is there an initial value
#' @param inits inital value
#' @param gene_weight gene weight
#' @param diagnostic convergence judgment results
#' @param logpdf The objective function value of each optimization

#'
#' @export
#'
CreateTrackOptimResultObject <- function(sample_value = NULL,
                                         dimnames = NULL,
                                         burnin = NULL,
                                         chains = NULL,
                                         iter = NULL,
                                         hasinits = FALSE,
                                         inits = list(),
                                         diagnostic = list(),
                                         logpdf = list(),
                                         gene_weight = list()
){

  # check sample_value and dimnames
  if(missing(x = sample_value) && missing(x = dimnames)){
    stop("Must provide either 'sample_value' or 'dimnames'")
  }

  # check hasinits and inits
  if(hasinits){
    if(missing(inits)){
      stop("has inits? if has, please provide it in 'inits = '")
    }
    if(!is.list(inits)){
      stop("inits must be list")
    }
  }

  # check inits
  object.list <- list(
    'Class' = 'TrackOptimResult',
    'sample_value' = sample_value,
    'dimnames' = dimnames,
    'chains' = chains,
    'burnin' = burnin,
    'iter' = iter,
    'hasinits' = hasinits,
    'inits' = inits,
    'diagnostic' = diagnostic,
    'logpdf' = logpdf,
    'gene_weight' = gene_weight
  )

  object <- do.call(what = 'new', args = object.list)
  return(object)
}

#' Create DimReducs Object
#' @description
#' create dimensionality reduction object
#' @param pca pca result, from 'RUNPCA'
#' @param tsne tsne result, from 'RUNtSNE'
#' @param umap umap result, from 'RUNUMAP'
#' @param dm dm result, from 'RUNDM'
#' @export
#'
CreateDimReducObject <- function(pca = NULL,
                                 tsne = NULL,
                                 umap = NULL,
                                 dm = NULL){

  # check sample_value and dimnames
  # if(missing(x = reduc_result)){
  #   stop("Must provide 'reduc_result'")
  # }

  # pca = new('Class' = 'PCAResult')
  # tsne = new('Class' = 'tSNEResult')
  # umap = new('Class' = 'UMAPResult')
  # dm = new('Class' = 'DMResult')

  # pca = NULL
  # tsne = NULL
  # umap = NULL
  # dm = NULL

  dimreduc = new(Class = 'DimReducs',
                 PCA = NULL,
                 tSNE = NULL,
                 DM = NULL,
                 UMAP = NULL)

  return(dimreduc)

}

#' Create DimReducs Object
#' @description
#' create settings object
#' @param datasetTitle project name
#' @param settings parameter settings
#' @export
#'
CreateSettings <- function(datasetTitle = NULL,
                           settings = NULL){

  Settings = new(Class = 'Settings',
                 datasetTitle = datasetTitle,
                 settings = list())
  return(Settings)
}

#' Create MGPfact Object
#' @description
#' create new object of MGPfact
#' @param data_matrix gene expression matrix, maybe the expression matrix after normalization or centralization
#' @param count_matrix original count matrix
#' @param dir path to create folder
#' @param datasetTitle project name
#' @param OptimResult optim result from julia package MGPfact
#' @param MURP MURP result from R package MURP
#' @param DimReducs dimensionality reduction results
#' @param BeautGene results related to genes
#' @param MetaData cell information
#' @param Tree trajectory toplog
#' @param GPR gaussian process regression
#' @param Command execution log
#' @export
#'
CreateMGPfactObject <- function(data_matrix = NULL,
                                count_matrix = NULL,
                                datasetTitle = "project",
                                dir = NULL,
                                OptimResult = NULL,
                                MURP = NULL,
                                DimReducs = CreateDimReducObject(),
                                BeautGene = list(),
                                MetaData = NULL,
                                Tree = NULL,
                                GPR = NULL,
                                Command = NULL){

  if(is.null(dir)){
    dir = "MGPfact_result"
    dir.create(dir)
  }else{
    if(!dir.exists(dir)) dir.create(dir)
  }
  setwd(dir)
  Initialize()
  cat("Create MGPfact Object in ", dir, "\n")

  # check count_matrix and data_matrix
  if(missing(x = data_matrix)){
    stop("Must provide either 'data_matrix'")
  }else if (!missing(x = count_matrix)) {
    if (!inherits(x = count_matrix, what = 'dgCMatrix')) {
      count_matrix <- as(object = as.matrix(x = count_matrix), Class = 'dgCMatrix')
    }
  }

  # save to 0_input
  save(data_matrix, file = "0_input/data_matrix.rda")
  save(MetaData, file = "0_input/MetaData.rda")

  # check metadata
  if(missing(MetaData)){
    MetaData <- data.frame(row.names = rownames(data_matrix),
                           cellname = rownames(data_matrix),
                           stringsAsFactors = FALSE)
  }else {
    if (is.null(x = rownames(x = MetaData))) {
      stop("Row names not set in metadata. Please ensure that rownames of metadata match column names of data matrix")
    }
    if (length(x = setdiff(x = rownames(x = MetaData), y = rownames(x = data_matrix)))) {
      warning("Some cells in meta.data not present in provided data_matrix matrix.")
      MetaData <- MetaData[intersect(x = rownames(x = MetaData), y = rownames(x = data_matrix)), , drop = FALSE]
    }
    if (is.data.frame(x = MetaData)) {
      # save(MetaData, file = "~/x.rda")
      new.MetaData <- data.frame(row.names = rownames(x = data_matrix))
      for (ii in 1:ncol(MetaData)) {
        # cat(ii, "\n")
        # cat(colnames(MetaData)[ii],"\n")
        new.MetaData[rownames(MetaData), colnames(x = MetaData)[ii]] <- MetaData[, ii, drop = FALSE]
      }
      MetaData <- new.MetaData
    }
  }

  # check count
  if (!missing(x = count_matrix)){
    # check colnames
    if (anyDuplicated(colnames(x = count_matrix)) | anyDuplicated(rownames(x = data_matrix))) {
      stop( "Non-unique cell names (colnames) present in the input matrix, making unique" )
    }
    # check rownames
    if (anyDuplicated(rownames(x = count_matrix)) | anyDuplicated(colnames(x = data_matrix))) {
      stop( "Non-unique cell names (colnames) present in the input matrix, making unique" )
    }
    if (nrow(x = count_matrix) > 0 && is.null(x = rownames(x = count_matrix))) {
      stop("No feature names (rownames) names present in the input matrix")
    }
    if (any(rownames(x = count_matrix) == '')) {
      stop("Feature names of count_matrix matrix cannot be empty", call. = FALSE)
    }
    if (is.null(x = colnames(x = count_matrix))) {
      stop("No cell names (colnames) names present in the input matrix")
    }

    # Ensure row- and column-names are vectors, not arrays
    if (!is.vector(x = rownames(x = count_matrix))) {
      rownames(x = count_matrix) <- as.vector(x = rownames(x = count_matrix))
    }
    if (!is.vector(x = colnames(x = count_matrix))) {
      colnames(x = count_matrix) <- as.vector(x = colnames(x = count_matrix))
    }
  }

  # create settings
  settings <- CreateSettings(datasetTitle = datasetTitle,
                             settings = NULL)

  # create MGPfact object
  MGPfact <- new(
    Class = 'MGPfact',
    assay = list(count_matrix = count_matrix, data_matrix = data_matrix),
    OptimResult = OptimResult,
    MURP = MURP,
    DimReducs = DimReducs,
    MetaData = MetaData,
    Date = date(),
    Settings = settings,
    BeautGene = list(),
    Tree = list(),
    GPR = list(),
    Command = list())
  if(!is.null(BeautGene)){
    MGPfact@BeautGene = BeautGene
  }
  if(!is.null(Tree)){
    MGPfact@Tree = Tree
  }
  if(!is.null(GPR)){
    MGPfact@GPR = GPR
  }
  return(MGPfact)
}

#' set parameters of settings
#'
#' @description
#' A folder about MGPfact will be created in the current path
#' @export
Initialize <- function(){
  # create directory in the local
  # dir.create("1_input", showWarnings=FALSE)
  dir.create("0_input", showWarnings=FALSE)
  dir.create("1_murp", showWarnings=FALSE)

  dir.create("2_pseudotime", showWarnings=FALSE)
  dir.create("2_pseudotime/2.1_julia_result", showWarnings=FALSE)
  dir.create("2_pseudotime/2.2_binarytree", showWarnings=FALSE)
  dir.create("2_pseudotime/2.3_tbtree", showWarnings=FALSE)
  dir.create("2_pseudotime/2.4_cov_matrix", showWarnings=FALSE)

  dir.create("3_tracking", showWarnings=FALSE)
  dir.create("3_tracking/3.1_optim_result", showWarnings=FALSE)
  dir.create("3_tracking/3.2_trajectory", showWarnings=FALSE)
  dir.create("3_tracking/3.3_gene_weight", showWarnings=FALSE)

  dir.create("4_differential_genes", showWarnings=FALSE)
  dir.create("5_combine_plot", showWarnings=FALSE)
}

#' set parameters of settings
#'
#' @description
#' Set parameters necessary in optimization
#'
#' @param object MGPfact object
#' @param murp_pc_number the number of principal components of the expression matrix after MURP downsampling
#' @param trajectory_number the number of trajectories to be deconstructed
#' @param pse_optim_iterations the number of optimization iterations to predict pseudotime
#' @param start_murp root point, If not required, set it to a number larger than the number of murps
#' @param chains_number the number of Markov chains
#'
#' @export
#'
SetSettings <- function(object,
                        murp_pc_number = 3,
                        trajectory_number = 3,
                        pse_optim_iterations = 10,
                        start_murp = 1,
                        chains_number = 10){

  setting = list(cell_number = nrow(object@assay$data_matrix),
                 gene_number = ncol(object@assay$data_matrix),
                 murp_number = object@MURP$Recommended_K,
                 murp_pc_number = murp_pc_number,
                 trajectory_number = trajectory_number,
                 pse_optim_iterations = pse_optim_iterations,
                 start_murp = start_murp,
                 chains_number = chains_number)
  new_c = setdiff(names(object@Settings@settings),names(setting))
  x = c(object@Settings@settings[new_c],setting)
  object@Settings@settings = x
  return(object)
}

#' assignSettings
#'
#' @description
#' add a single parameter
#' @param object MGPfact object
#' @param slotName slot name
#' @param values slot value
#' @export
#'
assignSettings <- function(object,
                           slotName = NULL,
                           values = NULL){
  if(slotName %in% names(object@Settings@settings)){
    # object@Settings@settings[[slotName]] = c(object@Settings@settings[[slotName]],values) %>% unique
    object@Settings@settings[[slotName]] = values
  }else{
    len = length(object@Settings@settings)
    object@Settings@settings[[len+1]] = values
    names(object@Settings@settings)[len+1] = slotName
  }
  return(object)
}

#' assignSettingsMulti
#'
#' @description
#' add multiple parameters
#'
#' @param object MGPfact object
#' @param slot a list containing parameters, must be named with the parameter name
#'
#' @export
#'
assignSettingsMulti <- function(object,
                                slot){
  for(i in 1:length(slot) ){
    slotName = names(slot)[i]
    values = slot[[i]]
    object = assignSettings(object, slotName, values)
  }
  return(object)
}

#' writeSettings
#'
#' @description
#' Write out the parameter settings for easy viewing
#'
#' @param object MGPfact object
#' @export
#'
writeSettings <- function(object){
  x = object@Settings@settings
  write.table("\n", file = paste0(getwd(),"/settings.txt"),
              quote = FALSE, append = FALSE, row.names = FALSE, col.names = FALSE, sep = "\n")
  for(i in 1:length(x)){
    write.table(paste0("*", names(x)[i]), file = paste0(getwd(),"/settings.txt"),
                quote = FALSE, append = TRUE, row.names = FALSE, col.names = FALSE, sep = "\n")
    write.table(paste0(x[[i]], collapse = ", "), file = paste0(getwd(),"/settings.txt"),
                quote = FALSE, append = TRUE, row.names = FALSE, col.names = FALSE, sep = "\n\n")
    write.table("\n", file = paste0(getwd(),"/settings.txt"),
                quote = FALSE, append = TRUE, row.names = FALSE, col.names = FALSE, sep = "\n")
  }
}


#' GetCommand
#'
#' @description
#' Get the parameter set of the running command,
#' which can be conveniently stored in "Command"
#' @noRd
GetCommand <- function(){
  time.stamp <- Sys.time()
  command.name = sys.calls()[[1]]
  command.name = strsplit(as.character(command.name[[1]]),"\\(")[[1]]

  argnames <- argnames <- names(x = formals(fun = sys.function(which = sys.parent(n = 1))))
  params <- list()
  p.env <- parent.frame(n = 1)
  argnames <- intersect(x = argnames, y = ls(name = p.env))
  argnames <- setdiff(argnames, c("mamba","murp","metadata","object"))
  for (arg in argnames) {
    param_value <- get(x = arg, envir = p.env)
    # if (inherits(x = param_value)) {
    #   next
    # }
    params[[arg]] <- param_value
  }

  # p = list(name = command.name, time.stamp = time.stamp, argnames = argnames, params = params)
  p = list(name = command.name, time.stamp = time.stamp, params = params)
  return(p)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                                 setMethod
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @export
setMethod("show", "MGPfact",
          function(object) {
            cat("An object of class ", class(object), "\n", sep = "")
            cat("Default data: data_matrix \n")
            cat(" ", ncol(object@assay$data_matrix), " genes across ",
                nrow(object@assay$data_matrix), " samples.\n", sep = "")
            invisible(NULL)
          }
)

#' @export
setMethod("show", "PseOptimResult",
          function(object) {
            P = length(grep("T\\[",object@dimnames[[2]]))
            L = length(grep("Yx",object@dimnames[[2]])) / P
            cat("An object of class ", class(object), "\n", sep = "")
            cat("Pseudo-time optimization results with", P, "MURPs \n")
            cat("Trajectory Number: ", L, "\n")
            cat("Iterations: ",object@iter," \n")
            cat("Chains Number: ", length(object@dimnames[[3]]), " \n")
            invisible(NULL)
          }
)

#' @export
setMethod("show", "TrackOptimResult",
          function(object) {
            # P = getParams(object,"murp_number")
            # L = 3
            cat("An object of class ", class(object), "\n", sep = "")
            # cat("Trajectory optimization results with", P, "MURPs \n")
            # cat("Trajectory Number: ", L, "\n")
            cat("Iterations: ",object@iter," \n")
            invisible(NULL)
          }
)

#' dim
#' @description
#' Returns the expression matrix in the MGPfact object and the dimensions of the downsampling matrix
#' @param x An MGPfact object.
#' @export
#'
setMethod("dim", "MGPfact",
          function(x) {
            d1 = dim(x@assay$data_matrix)
            d2 = dim(x@MURP$Recommended_K_cl$centers)
            cat("expression matrix: ", d1[1], " cells, ", d1[2], " genes ",  "\n", sep = "")
            cat("murp matrix: ", d2[1], " murps, ", d2[2], " genes ",  "\n", sep = "")
            return(d1)
          }
)

#' $
#' @description $
#' @param x An MGPfact object
#' @param name attribute columns
#' @export
setMethod("$", "MGPfact",
          function(x, name) {
            x = x@MetaData[,name]
            return(x)
            # invisible(NULL)
          }
)

#' $-MGPfact-method
#' @description Replace or add a column in the metadata
#' @param x An MGPfact object
#' @param name The column name
#' @param value The column values
setReplaceMethod("$", signature = "MGPfact",
                 function(x, name, value) {
                   x@MetaData = cbind(x@MetaData, value)
                   colnames(x@MetaData)[ncol(x@MetaData)] = name
                   return(x)
                 })

#' Get Params
#' @param x An MGPfact object
#' @param ... other parameters
#' @export
setGeneric(name="getParams",
           def=function(x,...) standardGeneric("getParams"))

#' getParams
#' @param x An MGPfact object
#' @param slotName slot name
#' @export
setMethod("getParams",
          signature="MGPfact",
          definition = function(x, slotName=NULL)
          {
            if(is.null(slotName))
            {
              print(x@Settings@settings)
            }else{
              return(x@Settings@settings[[slotName]])
            }
            invisible(NULL)
          })

#' getPse
#' @param x An MGPfact object
#' @export
setGeneric(name="getPse",
           def=function(x) standardGeneric("getPse"))

#' getPse
#' @param x An MGPfact object
#' @export
setMethod("getPse",
          signature="MGPfact",
          definition = function(x)
          {
            return(x@OptimResult$pse)
            invisible(NULL)
          })

#' getTrack
#' @param x An MGPfact object
#' @export
setGeneric(name="getTrack",
           def=function(x) standardGeneric("getTrack"))

#' getTrack
#' @param x An MGPfact object
#' @export
setMethod("getTrack",
          signature="MGPfact",
          definition = function(x)
          {
            return(x@OptimResult$track)
            invisible(NULL)
          })

# --------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                                   back up
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' The PCAResult-class
#'
#' @slot pca_result pca result
#' @slot MURP MURP result
#' @slot OptimResult mamba result from julia
#'
#' @name PCAResult-class
#' @rdname PCAResult-class
#' @exportClass PCAResult
#'
# PCAResult <- setClass(
#   Class = 'PCAResult',
#   slots = c(
#     pca_result = 'ANY',
#     MURP = 'ANY',
#     MambaResult = 'ANY'
#   )
# )

#' The tSNEResult-class
#'
#' @slot tsne_result tsne result
#' @slot MURP MURP result
#' @slot OptimResult mamba result from julia
#'
#' @name tSNEResult-class
#' @rdname tSNEResult-class
#' @exportClass tSNEResult
#'
# tSNEResult <- setClass(
#   Class = 'tSNEResult',
#   slots = c(
#     tsne_result = 'ANY',
#     MURP = 'ANY',
#     MambaResult = 'ANY'
#   )
# )

#' The DMResult-class
#'
#' @slot dm_result dm result
#' @slot MURP MURP result
#' @slot OptimResult mamba result from julia
#'
#' @name DMResult-class
#' @rdname DMResult-class
#' @exportClass DMResult
#'
# DMResult <- setClass(
#   Class = 'DMResult',
#   slots = c(
#     dm_result = 'ANY',
#     MURP = 'ANY',
#     MambaResult = 'ANY'
#   )
# )

#' The UMAPResult-class
#'
#' @slot umap_result umap result
#' @slot MURP MURP result
#' @slot OptimResult mamba result from julia
#'
#' @name UMAPResult-class
#' @rdname UMAPResult-class
#' @exportClass UMAPResult
#'
# UMAPResult <- setClass(
#   Class = 'UMAPResult',
#   slots = c(
#     umap_result = 'ANY',
#     MURP = 'ANY',
#     MambaResult = 'ANY'
#   )
# )

#' The OptimeResult-class
#'
#' @slot sample_value sampling results for each iteration
#' @slot dimnames names of sampling value
#' @slot burnin
#' @slot chains
#' @slot iter iteration number
#' @slot hasinits
#' @slot inits initial values
#' @slot iter_chain_cor
#' @slot t_pred
#' @slot sortindex_cor
#' @slot diagnostic
#' @slot logpdf
#'
#' @exportClass OptimResult
#'
# AllOptimResult <- setClass(
#   Class = 'AllOptimResult',
#   slots = c(
#     sample_value = 'array',
#     dimnames = 'list',
#     burnin = 'numeric',
#     chains = 'vector',
#     iter = 'numeric',
#     hasinits = 'logical',
#     inits = 'list',
#     iter_chain_cor = 'list',
#     t_pred = 'list',
#     sortindex_cor = 'list',
#     diagnostic = 'list',
#     logpdf = 'list'
#   )
# )

#' Create PseMambaResult Object
#'
#' @param sample_value sampling results for each iteration
#' @param dimnames names of sampling value
#' @param burnin burning period
#' @param chains chain number,vector
#' @param iter iteration number
#' @param hasinits is there an initial value
#' @param inits inits
#'
#' @export
#'
# CreateAllOptimResultObject <- function(sample_value = NULL,
#                                        dimnames = NULL,
#                                        burnin = NULL,
#                                        chains = NULL,
#                                        iter = NULL,
#                                        hasinits = FALSE,
#                                        inits = list(),
#                                        iter_chain_cor = list(),
#                                        sortindex_cor = list(),
#                                        t_pred = list(),
#                                        diagnostic = list(),
#                                        logpdf = list()
# ){
#
#   # check sample_value and dimnames
#   if(missing(x = sample_value) && missing(x = dimnames)){
#     stop("Must provide either 'sample_value' or 'dimnames'")
#   }
#
#   # check hasinits and inits
#   if(hasinits){
#     if(missing(inits)){
#       stop("has inits? if has, please provide it in 'inits = '")
#     }
#     if(!is.list(inits)){
#       stop("inits must be list")
#     }
#   }
#
#   # check inits
#   object.list <- list(
#     'Class' = 'AllOptimResult',
#     'sample_value' = sample_value,
#     'dimnames' = dimnames,
#     'chains' = chains,
#     'burnin' = burnin,
#     'iter' = iter,
#     'hasinits' = hasinits,
#     'inits' = inits,
#     'iter_chain_cor' = iter_chain_cor,
#     't_pred' = t_pred,
#     'sortindex_cor' = sortindex_cor,
#     'diagnostic' = diagnostic,
#     'logpdf' = logpdf
#   )
#
#   object <- do.call(what = 'new', args = object.list)
#   return(object)
# }

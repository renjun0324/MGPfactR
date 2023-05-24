

## runSCENIC_2_createRegulons, mingene = 5
## runSCENIC_3_scoreCells, ncores = 1; keey regulon>3 stay; if regulon<2, rm hclust;
## runSCENIC_4_aucell_binarize
library(data.table)
library(foreach)

# Help files will be automatically generated from the coments starting with #'
# (https://cran.r-project.org/web/packages/roxygen2/vignettes/rd.html)

#' @importFrom stats sd setNames
#'
#' @title Add motif annotation
#' @description Select significant motifs and/or annotate motifs to
#' genes or transcription factors.
#' The motifs are considered significantly enriched if they pass the the
#' Normalized Enrichment Score (NES) threshold.
#' @param auc Output from calcAUC.
#' @param nesThreshold Numeric. NES threshold to calculate the motif significant
#' (3.0 by default). The NES is calculated -for each motif- based on the AUC
#' distribution of all the motifs for the gene-set [(x-mean)/sd].
#' @param digits Integer. Number of digits for the AUC and NES in the
#' output table.
#' @param motifAnnot Motif annotation database containing the 
#' annotations of the motif to transcription factors.
#' The names should match the ranking column names.
#' @param motifAnnot_highConfCat Categories considered as source for 
#' 'high confidence' annotations. By default, 
#' "directAnnotation" (annotated in the source database), and 
#' "inferredBy_Orthology" (the motif is annotated to an homologous/ortologous 
#' gene).
#' @param motifAnnot_lowConfCat Categories considered 
#' 'lower confidence' source for annotations. By default, the annotations 
#' inferred based on motif similarity ("inferredBy_MotifSimilarity", 
#' "inferredBy_MotifSimilarity_n_Orthology").
#' @param idColumn Annotation column containing the ID (e.g. motif, accession)
#' @param highlightTFs Character. If a list of transcription factors is
#' provided, the column TFinDB in the otuput table will indicate whether any
#' of those TFs are included within the 'high-confidence' annotation 
#' (two asterisks, **)
#' or 'low-confidence' annotation (one asterisk, *) of the motif.
#' The vector can be named to indicate which TF to highlight for each gene-set.
#' Otherwise, all TFs will be used for all geneSets.
#' @param keepAnnotationCategory Include annotation type in the TF information?
#' @return \code{\link[data.table]{data.table}} with the folowing columns:
#' \itemize{
#' \item geneSet: Name of the gene set
#' \item motif: ID of the motif
#' (colnames of the ranking, it might be other kind of feature)
#' \item NES: Normalized enrichment score of the motif in the gene-set
#' \item AUC: Area Under the Curve (used to calculate the NES)
#' \item TFinDB: Indicates whether the highlightedTFs are included within the
#' high-confidence annotation (two asterisks, **)
#' or lower-confidence annotation (one asterisk, *)
#' \item TF_highConf: Transcription factors annotated to the motif 
#' based on high-confidence annotations.
#' \item TF_lowConf: Transcription factors annotated to the motif according to
#' based on lower-confidence annotations.
#' }
#' @seealso Next step in the workflow: \code{\link{addSignificantGenes}}.
#'
#' Previous step in the workflow: \code{\link{calcAUC}}.
#'
#' See the package vignette for examples and more details:
#' \code{vignette("RcisTarget")}
#' @example inst/examples/example_addMotifAnnotation.R
#' @export
addMotifAnnotation <- function(auc, nesThreshold=3.0, digits=3,
                               motifAnnot=NULL, 
                               motifAnnot_highConfCat=c("directAnnotation", "inferredBy_Orthology"), 
                               motifAnnot_lowConfCat=c("inferredBy_MotifSimilarity", 
                                                       "inferredBy_MotifSimilarity_n_Orthology"), 
                               idColumn="motif",
                               highlightTFs=NULL,
                               keepAnnotationCategory=TRUE)
{
  auc <- getAUC(auc)
  #### Check inputs
  if(!is.null(highlightTFs))
  {
    if(is.null(motifAnnot))
      stop("To hightlight TFs, please provide a motif-TF annotation.")
    if(is.null(names(highlightTFs))) {
      warning("The input TFs are not named, ",
              "all TFs will be used with all Gene Sets.")
      highlightTFs <- setNames(rep(list(highlightTFs), nrow(auc)),
                               rownames(auc))
    }
    
    if(!all(names(highlightTFs) %in% rownames(auc))) warning("TFs 1")
    if(!all(rownames(auc) %in% names(highlightTFs))) warning("TFs 2")
  }
  
  if(!is.null(motifAnnot))
  {
    if(!is.data.table(motifAnnot))
      stop("motifAnnot should be a data.table")
    # if(!is.null(motifAnnot_highConfCat) && 
    #    any(!motifAnnot_highConfCat %in% levels(motifAnnot$annotationSource)))
    #   warning("'motifAnnot_highConfCat' 
    #        should be a value in the column 'annotationSource'.") 
    # 
    # if(!is.null(motifAnnot_lowConfCat) && 
    #    any(!motifAnnot_lowConfCat %in% levels(motifAnnot$annotationSource)))
    #   warning("'motifAnnot_lowConfCat' 
    #        should be a value in the column 'annotationSource'.") 
    
    commonCat <- intersect(motifAnnot_highConfCat, motifAnnot_lowConfCat)
    if(length(commonCat)>0)
      warning("The following annotation types are both in", 
              "'motifAnnot_highConfCat' and 'motifAnnot_lowConfCat': ", commonCat)
  }
  
  #### Runs "auc.asTable" on each signature/geneset
  applyFun <- lapply
  if((nrow(auc)>100000) && ("BiocParallel" %in% installed.packages())) 
  {
    applyFun <- BiocParallel::bplapply
    message("Using BiocParallel...")
  }
  
  ret <- applyFun(rownames(auc), function(geneSet) {
    tfs <- highlightTFs[[geneSet]]
    aucTable <- .auc.asTable(auc[geneSet,],
                             nesThreshold=nesThreshold, digits=digits, idColumn=idColumn)
    if(nrow(aucTable)>0)
    {
      aucTable <- .addTfs(aucTable,
                          motifAnnot=motifAnnot,
                          TFs=tfs, 
                          motifAnnot_highConfCat=motifAnnot_highConfCat,
                          motifAnnot_lowConfCat=motifAnnot_lowConfCat,
                          idColumn=idColumn,
                          keepAnnotationCategory=keepAnnotationCategory)
      aucTable <- data.table::data.table(geneSet=geneSet, aucTable)
    }else{
      aucTable <- NULL
    }
    aucTable
  })
  
  ## Merge the results from each signature/geneSet/regionSet into a single dt
  # ret <- do.call(rbind, unname(ret))  # Slower?
  # library(data.table)
  ret <- data.table::rbindlist(ret)
  
  if(nrow(ret)>0)
    colnames(ret)[which(colnames(ret) == "ranking")] <- "motif"
  return(ret)
}


############ PRIVATE
.calcNES <- function(AUC)
{
  meanAUC <- mean(AUC)
  sdAUC <- sd(AUC)
  
  # NES = (AUC-mean)/sd
  NES <- vapply(AUC, function(x) (x-meanAUC)/sdAUC,
                FUN.VALUE=numeric(1))
  return(NES)
}

.auc.asTable <- function(auc, nesThreshold=3.0, digits=3, idColumn="motif")
{
  nes <- .calcNES(auc)
  nes <- sort(nes, decreasing=TRUE)
  
  signifRankings <- names(nes)[which(nes >= nesThreshold)]
  aucTable <- data.table::data.table(motif=signifRankings,
                                     NES=signif(nes[signifRankings], digits=digits),
                                     AUC=signif(auc[signifRankings],digits=digits))
  colnames(aucTable)[1] <- idColumn
  aucTable
}

.addTfs <- function(aucTable,
                    motifAnnot=NULL,
                    TFs=NULL,
                    motifAnnot_highConfCat=NULL,
                    motifAnnot_lowConfCat=NULL,
                    idColumn="motif",
                    keepAnnotationCategory=keepAnnotationCategory)
{
  if(!is.null(TFs))
  {
    aucTable <- data.table::data.table(aucTable,
                                       highlightedTFs=paste(TFs, collapse=", ") ,
                                       TFinDB="")
    
    if(!is.null(motifAnnot)) {
      motifAnnot_subset <- motifAnnot[(motifAnnot[[idColumn]] %in% aucTable[[idColumn]]) 
                                      & (motifAnnot$TF %in% TFs), ][,c(idColumn, "TF", "annotationSource"),with=F]
      motifAnnot_subset <- split(motifAnnot_subset, motifAnnot_subset[[idColumn]])
      for(motifName in names(motifAnnot_subset))
      {
        if(any(as.character(motifAnnot_subset[[motifName]]$annotationSource) 
               %in% motifAnnot_lowConfCat))
          aucTable[aucTable[[idColumn]]==motifName,"TFinDB"] <- "*"
        
        # overrides lowConf
        if(any(as.character(motifAnnot_subset[[motifName]]$annotationSource) 
               %in% motifAnnot_highConfCat))
          aucTable[aucTable[[idColumn]]==motifName,"TFinDB"] <- "**"
      }
    }
  }
  
  if(!is.null(motifAnnot))
  {
    if(!is.null(motifAnnot_highConfCat))
    {
      TF_highConf <- .formatTfs(motifs=aucTable[[idColumn]], 
                                motifAnnot=motifAnnot,
                                annotCats=motifAnnot_highConfCat,
                                idColumn=idColumn,
                                keepAnnotationCategory=keepAnnotationCategory)
      
      aucTable <- data.table::data.table(aucTable, TF_highConf=TF_highConf)
    }
    
    if(!is.null(motifAnnot_lowConfCat))
    {
      TF_lowConf <- .formatTfs(motifs=aucTable[[idColumn]], 
                               motifAnnot=motifAnnot,
                               annotCats=motifAnnot_lowConfCat,
                               idColumn=idColumn, 
                               keepAnnotationCategory=keepAnnotationCategory)
      
      aucTable <- data.table::data.table(aucTable, TF_lowConf=TF_lowConf)
    }
  }
  
  aucTable
}


## 26 apr 2019
# Replaced input: aucTable by motifs. In calls:  .formatTfs(motifs=aucTable[[idColumn]]
# aucTable$motif --> motifs
# nrow(aucTable) --> length(motifs)
.formatTfs <- function(motifs, motifAnnot, annotCats, idColumn, keepAnnotationCategory)
{
  motifAnnot_subset <- motifAnnot[motifAnnot$annotationSource %in% annotCats, ] 
  motifAnnot_subset <- motifAnnot_subset[motifAnnot_subset[[idColumn]] %in% motifs, ] 
  if(keepAnnotationCategory)
  {
    motifAnnot_Cats <- vapply(split(motifAnnot_subset, motifAnnot_subset[[idColumn]]), 
                              function(mat){
                                mat <- split(mat$TF, factor(mat$annotationSource))
                                tfsByCat <- vapply(names(mat),
                                                   function(x) paste(paste(unlist(mat[[x]]),
                                                                           collapse="; "),
                                                                     " (",x,"). ",
                                                                     sep=""), "")
                                paste(tfsByCat, collapse="")
                              }, FUN.VALUE="")
  }else{
    motifAnnot_Cats <- vapply(split(motifAnnot_subset, motifAnnot_subset[[idColumn]]), 
                              function(mat){
                                paste(unique(unlist(mat$TF)), collapse="; ")
                              }, FUN.VALUE="")
  }
  
  ret <- setNames(rep("", length(motifs)), motifs)
  ret[names(motifAnnot_Cats)] <- motifAnnot_Cats
  return(ret)
}

#' @title Get motif annotation
#' @description Get the genes/transcription factors annotated to the given motifs

#' @param motifs Motif IDs 
#' @param motifAnnot Motif annotation database containing the 
#' annotations of the motif to genes or transcription factors.
#' @param annotCats Annotation categories to be considered:  
#' "directAnnotation" (annotated in the source database),  
#' "inferredBy_Orthology" (the motif is annotated to an homologous/ortologous 
#' gene), or inferred based on motif similarity ("inferredBy_MotifSimilarity", 
#' "inferredBy_MotifSimilarity_n_Orthology").
#' @param idColumn Annotation column containing the ID (e.g. motif, accession)
#' @param returnFormat Determines the output format. Choose one of the following values:
#' @param keepAnnotationCategory Include annotation type in the TF information?
#' \itemize{
#' \item \code{asCharacter}: Named vector with the genes or TFs annotated to the given motifs (in the same order, including empty and duplicated values).
#' \item \code{subset}: Subset of the annotation table (list split by motif)
#' \item \code{list}: List of TF names (unique values), duplicated motifs or motifs without annotation are not returned.
#' }
#' @return See argument \code{returnFormat}
#' @seealso \code{\link{addMotifAnnotation}} add the annotation directly to the motif enrichment results.
#' 
#' See the package vignette for examples and more details:
#' \code{vignette("RcisTarget")}
#' @example inst/examples/example_addMotifAnnotation.R
#' @export
getMotifAnnotation <- function(motifs, 
                               motifAnnot, 
                               annotCats=c("directAnnotation",
                                           "inferredBy_MotifSimilarity",
                                           "inferredBy_Orthology",
                                           "inferredBy_MotifSimilarity_n_Orthology"),
                               idColumn="motif",
                               returnFormat=c("asCharacter","subset","list")[1],
                               keepAnnotationCategory=TRUE)
{
  ## Check inputs:
  returnFormat <- tolower(returnFormat)
  if(!returnFormat %in% c("ascharacter","subset","list")) stop("returnFormat should be eihter 'asCharacter', 'subset', or 'list'.")
  if(length(returnFormat)>1) stop("Please, choose ONE returnFormat.")
  
  ## Run:
  if(returnFormat=="ascharacter"){
    ret <- .formatTfs(motifs=motifs, 
                      motifAnnot=motifAnnot,
                      annotCats=annotCats,
                      idColumn=idColumn,
                      keepAnnotationCategory=keepAnnotationCategory) 
  }else{
    ret <- .getTfs(motifs=motifs, 
                   motifAnnot=motifAnnot,
                   annotCats=annotCats,
                   idColumn=idColumn,
                   returnFormat=returnFormat) 
  }
  return(ret)
}

.getTfs <- function(motifs, motifAnnot, annotCats, idColumn, returnFormat)
{
  motifAnnot_subset <- motifAnnot[motifAnnot$annotationSource %in% annotCats, ] 
  motifAnnot_subset <- motifAnnot_subset[motifAnnot_subset[[idColumn]] %in% motifs, ] 
  # motifAnnot_subset <- split(motifAnnot_subset[,c("TF", "directAnnotation", "inferred_Orthology", "inferred_MotifSimil","annotationSource")], motifAnnot_subset[[idColumn]])
  motifAnnot_subset <- split(motifAnnot_subset, motifAnnot_subset[[idColumn]])
  
  ret <- motifAnnot_subset # returnFormat=="subset"
  
  if(returnFormat=="list")
    ret <- lapply(ret, function(x) sort(unique(x[["TF"]])))
  
  return(ret)
}



################################################################################
# Step 2. Identifying regulons (direct TF targets) based on DNA motif enrichment
################################################################################


#' @title runSCENIC_2_createRegulons
#' @description Step 2: RcisTarget (prune co-expression modules using TF-motif enrichment analysis)
#' @param scenicOptions Fields used: TODO
#' @param minGenes Minimum size of co-expression gene set (default: 20 genes)
#' @param signifGenesMethod Method for 'addSignificantGenes'
#' @param coexMethods Allows to select the method(s) used to generate the co-expression modules
#' @param minJakkardInd Merge overlapping modules (with Jakkard index >=minJakkardInd; reduces running time).
#' @param onlyPositiveCorr Whether to include only positive-correlated targets in the regulons (default: TRUE).
#' @param dbIndexCol Column containing the feature name in the database (default: 'features')
#' @return The output is written in the folders 'int' and 'ouput'
#' @details See the detailed vignette explaining the internal steps.
#' @examples 
#' scenicOptions <- readRDS("int/scenicOptions.Rds")
#' # In case any settings need to be modified:
#' scenicOptions@settings$nCores <- 20
#' scenicOptions@inputDatasetInfo$org <- "mgi" 
#' 
#' runSCENIC_2_createRegulons(scenicOptions)
#' @export
runSCENIC_2_createRegulons <- function(scenicOptions, 
                                       minGenes=5, 
                                       coexMethods=NULL, 
                                       minJakkardInd=0.8,
                                       signifGenesMethod="aprox", 
                                       onlyPositiveCorr=TRUE,
                                       onlyBestGsPerMotif=TRUE,
                                       dbIndexCol='features'
)
{
  nCores <- getSettings(scenicOptions, "nCores")
  
  tfModules_asDF <- tryCatch(loadInt(scenicOptions, "tfModules_asDF"),
                             error = function(e) {
                               if(getStatus(scenicOptions, asID=TRUE) < 2) 
                                 e$message <- paste0("It seems the co-expression modules have not been built yet. Please, run runSCENIC_1_coexNetwork2modules() first.\n", 
                                                     e$message)
                               stop(e)
                             })
  if(!is.null(coexMethods)) tfModules_asDF <- tfModules_asDF[which(tfModules_asDF$method %in% coexMethods),]
  if(!is.null(minJakkardInd)) tfModules_asDF <- mergeOverlappingModules(tfModules_asDF, minJakkardInd=minJakkardInd) # New
  if(nrow(tfModules_asDF)==0) stop("The co-expression modules are empty.")
  
  # Set cores for RcisTarget::addMotifAnnotation(). The other functions use foreach package.
  if("BiocParallel" %in% installed.packages() && (nCores>1)) {
    library(BiocParallel)
    register(MulticoreParam(nCores), default=TRUE) 
  }
  
  msg <- paste0(format(Sys.time(), "%H:%M"), "\tStep 2. Identifying regulons")
  if(getSettings(scenicOptions, "verbose")) message(msg)
  
  ### Check load DBs
  library(AUCell)
  library(RcisTarget)
  motifAnnot <- getDbAnnotations(scenicOptions)
  
  if(is.null(names(getSettings(scenicOptions, "dbs")))) 
  {
    names(scenicOptions@settings$"dbs") <- scenicOptions@settings$"dbs"
    tmp <- sapply(strsplit(getSettings(scenicOptions, "dbs"),"-", fixed=T), function(x) x[grep("bp|kb",x)])
    if(all(lengths(tmp)>0)) names(scenicOptions@settings$"dbs") <- tmp
  }
  
  loadAttempt <- sapply(getDatabases(scenicOptions), dbLoadingAttempt, indexCol=dbIndexCol)
  if(any(!loadAttempt)) stop("It is not possible to load the following databses: \n",
                             paste(dbs[which(!loadAttempt)], collapse="\n"))
  
  genesInDb <- unique(unlist(lapply(getDatabases(scenicOptions), function(dbFilePath) {
    rf <- arrow::ReadableFile$create(dbFilePath)
    fr <- arrow::FeatherReader$create(rf)
    genesInDb <- names(fr)
    rnktype <- "features"        #TODO: add as option for custom dbs
    genesInDb <- genesInDb[genesInDb != rnktype]
  })))
  
  ## Check if annotation and rankings (potentially) match:
  featuresWithAnnot <- checkAnnots(scenicOptions, motifAnnot)
  if(any(featuresWithAnnot == 0)) warning("Missing annotations\n", names(which(rankingsInDb==0)))
  
  ### Filter & format co-expression modules
  # Remove genes missing from RcisTarget databases
  #  (In case the input matrix wasn't already filtered)
  tfModules_asDF$TF <- as.character(tfModules_asDF$TF)
  tfModules_asDF$Target <- as.character(tfModules_asDF$Target)
  allTFs <- getDbTfs(scenicOptions)
  tfModules_asDF <- tfModules_asDF[which(tfModules_asDF$TF %in% allTFs),]
  geneInDb <- tfModules_asDF$Target %in% genesInDb
  missingGene <- sort(unique(tfModules_asDF[which(!geneInDb),"Target"]))
  if(length(missingGene)>0) 
    warning(paste0("Genes in co-expression modules not available in RcisTargetDatabases: ", 
                   paste(missingGene, collapse=", ")))
  tfModules_asDF <- tfModules_asDF[which(geneInDb),]
  
  ######
  # Targets with positive correlation
  if(all(is.na(tfModules_asDF$corr)))
  {
    warning("no correlation info available") 
    tfModules_Selected <- tfModules_asDF
    tfModules_Selected$geneSetName <- paste(tfModules_Selected$TF, tfModules_Selected$method, sep="_")
  }else{
    tfModules_Selected <- tfModules_asDF[which(tfModules_asDF$corr==1),]
    tfModules_Selected$geneSetName <- paste(tfModules_Selected$TF, tfModules_Selected$method, sep="_")
    
    if(!onlyPositiveCorr)
    {
      tfModules_IgnCorr <- tfModules_asDF[which(tfModules_asDF$corr!=1),]
      tfModules_IgnCorr$geneSetName <- paste0(tfModules_IgnCorr$TF,"_", tfModules_IgnCorr$method)
      
      # Include positive corrs for these geneSets: 
      # gplots::venn(list(pos=unique(tfModules_SelectedgeneSetName),ign=unique(tfModulesIgnCorrgeneSetName)))
      posCorr <- tfModules_Selected[which(tfModules_Selected$geneSetName %in% unique(tfModules_IgnCorr$geneSetName)),]
      
      tfModules_IgnCorr <- rbind(tfModules_IgnCorr, posCorr)
      tfModules_IgnCorr$geneSetName <- paste0(tfModules_IgnCorr$geneSetName, "IgnCorr")
      
      tfModules_Selected <- rbind(tfModules_Selected, tfModules_IgnCorr)
    }
  }
  
  tfModules_Selected$geneSetName <- factor(as.character(tfModules_Selected$geneSetName))
  # head(tfModules_Selected)
  allGenes <- unique(tfModules_Selected$Target)
  
  
  #####
  # Split into tfModules (TF-modules, with several methods)
  tfModules <- split(tfModules_Selected$Target, tfModules_Selected$geneSetName)
  
  # Add TF to the gene set (used in the following steps, careful if editing)
  tfModules <- setNames(lapply(names(tfModules), function(gsn) {
    tf <- strsplit(gsn, "_")[[1]][1]
    unique(c(tf, tfModules[[gsn]]))
  }), names(tfModules))
  
  # Keep gene sets with at least 'minGenes' genes
  tfModules <- tfModules[which(lengths(tfModules)>=minGenes)]
  saveRDS(tfModules, file=getIntName(scenicOptions, "tfModules_forEnrichment")) #TODO as geneset? & previous step?
  
  if(getSettings(scenicOptions, "verbose")) {
    tfModulesSummary <- t(sapply(strsplit(names(tfModules), "_"), function(x) x[1:2]))
    message("tfModulesSummary:")
    print(cbind(sort(table(tfModulesSummary[,2]))))
  }
  
  ################################################################
  ### 1. Calculate motif enrichment for each TF-module (Run RcisTarget)
  
  ### 1.1 Calculate enrichment
  msg <- paste0(format(Sys.time(), "%H:%M"), "\tRcisTarget: Calculating AUC")
  if(getSettings(scenicOptions, "verbose")) message(msg)
  
  motifs_AUC <- lapply(getDatabases(scenicOptions), function(rnkName) {
    ranking <- importRankings(rnkName, columns=allGenes)
    message("Scoring database: ", ranking@description)
    RcisTarget::calcAUC(tfModules, ranking, aucMaxRank=0.03*getNumColsInDB(ranking), nCores=nCores, verbose=FALSE)})
  saveRDS(motifs_AUC, file=getIntName(scenicOptions, "motifs_AUC"))
  
  ### 1.2 Convert to table, filter by NES & add the TFs to which the motif is annotated
  # (For each database...)
  msg <- paste0(format(Sys.time(), "%H:%M"), "\tRcisTarget: Adding motif annotation")
  message(msg)
  motifEnrichment <- lapply(motifs_AUC, function(aucOutput)
  {
    # Extract the TF of the gene-set name (i.e. MITF_w001):
    tf <- sapply(setNames(strsplit(rownames(aucOutput), "_"), rownames(aucOutput)), function(x) x[[1]])
    
    # Calculate NES and add motif annotation (provide tf in 'highlightTFs'):
    addMotifAnnotation(aucOutput, 
                       nesThreshold=3, digits=3, 
                       motifAnnot=motifAnnot,
                       motifAnnot_highConfCat=c("directAnnotation", "inferredBy_Orthology"),
                       motifAnnot_lowConfCat=c("inferredBy_MotifSimilarity",
                                               "inferredBy_MotifSimilarity_n_Orthology"), 
                       highlightTFs=tf)
  })
  
  # Merge both tables, adding a column that contains the 'motifDb'
  motifEnrichment <- do.call(rbind, lapply(names(motifEnrichment), function(dbName){
    cbind(motifDb=dbName, motifEnrichment[[dbName]])
  }))
  saveRDS(motifEnrichment, file=getIntName(scenicOptions, "motifEnrichment_full"))
  msg <- paste0("Number of motifs in the initial enrichment: ", nrow(motifEnrichment))
  if(getSettings(scenicOptions, "verbose")) message(msg)
  
  ### 1.3 Keep only the motifs annotated to the initial TF
  motifEnrichment_selfMotifs <- motifEnrichment[which(motifEnrichment$TFinDB != ""),, drop=FALSE]
  msg <- paste0("Number of motifs annotated to the matching TF: ", nrow(motifEnrichment_selfMotifs))
  if(getSettings(scenicOptions, "verbose")) message(msg)
  rm(motifEnrichment)
  
  if(nrow(motifEnrichment_selfMotifs)==0) 
    stop("None of the co-expression modules present enrichment of the TF motif: There are no regulons.")
  
  ####
  if(onlyBestGsPerMotif)
  {
    met_byDb <- split(motifEnrichment_selfMotifs, motifEnrichment_selfMotifs$motifDb)
    for(db in names(met_byDb))
    {
      met <- met_byDb[[db]]
      met <- split(met, factor(met$highlightedTFs))
      met <- lapply(met, function(x){
        rbindlist(lapply(split(x, x$motif), function(y) y[which.max(y$NES),]))
      })
      # sapply(met, nrow)
      met_byDb[[db]] <- rbindlist(met)
    }
    motifEnrichment_selfMotifs <- rbindlist(met_byDb)
    rm(met_byDb); rm(met)
  }
  ####
  
  
  ################################################################
  # 2. Prune targets
  msg <- paste0(format(Sys.time(), "%H:%M"), "\tRcisTarget: Pruning targets")
  if(getSettings(scenicOptions, "verbose")) message(msg)
  
  dbNames <- getDatabases(scenicOptions)
  motifEnrichment_selfMotifs_wGenes <- lapply(names(dbNames), function(motifDbName){
    ranking <- importRankings(dbNames[motifDbName], columns=allGenes)
    addSignificantGenes(resultsTable=motifEnrichment_selfMotifs[motifEnrichment_selfMotifs$motifDb==motifDbName,],
                        geneSets=tfModules,
                        rankings=ranking, 
                        plotCurve = FALSE,
                        maxRank=5000, 
                        method=signifGenesMethod, 
                        nMean=100,
                        nCores=nCores)
  })
  
  suppressPackageStartupMessages(library(data.table))
  motifEnrichment_selfMotifs_wGenes <- rbindlist(motifEnrichment_selfMotifs_wGenes)
  saveRDS(motifEnrichment_selfMotifs_wGenes, file=getIntName(scenicOptions, "motifEnrichment_selfMotifs_wGenes"))
  
  if(getSettings(scenicOptions, "verbose")) 
  {
    # TODO messages/print
    message(format(Sys.time(), "%H:%M"), "\tNumber of motifs that support the regulons: ", nrow(motifEnrichment_selfMotifs_wGenes))
    motifEnrichment_selfMotifs_wGenes[order(motifEnrichment_selfMotifs_wGenes$NES,decreasing=TRUE),][1:5,(1:ncol(motifEnrichment_selfMotifs_wGenes)-1), with=F] 
  }
  
  # Save as text:
  if(!file.exists("output")) dir.create("output") 
  write.table(motifEnrichment_selfMotifs_wGenes, file=getOutName(scenicOptions, "s2_motifEnrichment"),
              sep="\t", quote=FALSE, row.names=FALSE)
  if("DT" %in% installed.packages() && nrow(motifEnrichment_selfMotifs_wGenes)>0)
  {
    nvm <- tryCatch({
      colsToShow <- c("motifDb", "logo", "NES", "geneSet", "TF_highConf", "TF_lowConf")
      motifEnrichment_2html <- viewMotifs(motifEnrichment_selfMotifs_wGenes, colsToShow=colsToShow, options=list(pageLength=100))
      
      fileName <- getOutName(scenicOptions, "s2_motifEnrichmentHtml")
      
      dirName <- dirname(fileName)
      fileName <- basename(fileName)
      suppressWarnings(DT::saveWidget(motifEnrichment_2html, fileName))
      file.rename(fileName, file.path(dirName, fileName))
      if(getSettings(scenicOptions, "verbose")) message("\tPreview of motif enrichment saved as: ", file.path(dirName, fileName))
    }, error = function(e) print(e$message))
  }
  
  ################################################################
  # Format regulons & save
  motifEnrichment.asIncidList <- apply(motifEnrichment_selfMotifs_wGenes, 1, function(oneMotifRow) {
    genes <- strsplit(oneMotifRow["enrichedGenes"], ";")[[1]]
    oneMotifRow <- data.frame(rbind(oneMotifRow), stringsAsFactors=FALSE)
    data.frame(oneMotifRow[rep(1, length(genes)),c("NES", "motif", "highlightedTFs", "TFinDB", "geneSet", "motifDb")], genes, stringsAsFactors = FALSE)
  })
  motifEnrichment.asIncidList <- rbindlist(motifEnrichment.asIncidList)
  # colnames(motifEnrichment.asIncidList) <- c("NES", "motif", "TF", "annot", "gene", "motifDb", "geneSet")
  colnames(motifEnrichment.asIncidList)[which(colnames(motifEnrichment.asIncidList)=="highlightedTFs")] <- "TF"
  colnames(motifEnrichment.asIncidList)[which(colnames(motifEnrichment.asIncidList)=="TFinDB")] <- "annot"
  colnames(motifEnrichment.asIncidList)[which(colnames(motifEnrichment.asIncidList)=="genes")] <- "gene"
  
  motifEnrichment.asIncidList <- data.frame(motifEnrichment.asIncidList, stringsAsFactors = FALSE)
  
  # Get targets for each TF, but keep info about best motif/enrichment
  # (directly annotated motifs are considered better)
  regulonTargetsInfo <- lapply(split(motifEnrichment.asIncidList, motifEnrichment.asIncidList$TF), function(tfTargets){
    # print(unique(tfTargets$TF))
    tfTable <- as.data.frame(do.call(rbind, lapply(split(tfTargets, tfTargets$gene), function(enrOneGene){
      highConfAnnot <- "**" %in% enrOneGene$annot
      enrOneGeneByAnnot <- enrOneGene
      if(highConfAnnot) enrOneGeneByAnnot <- enrOneGeneByAnnot[which(enrOneGene$annot == "**"),]
      bestMotif <- which.max(enrOneGeneByAnnot$NES)
      
      tf <- unique(enrOneGene$TF)
      cbind(TF=tf, 
            gene=unique(enrOneGene$gene), 
            highConfAnnot=highConfAnnot,
            nMotifs=nrow(enrOneGene),
            bestMotif=as.character(enrOneGeneByAnnot[bestMotif,"motif"]), NES=as.numeric(enrOneGeneByAnnot[bestMotif,"NES"]), 
            motifDb=as.character(enrOneGeneByAnnot[bestMotif,"motifDb"]), coexModule=gsub(paste0(tf,"_"), "", as.character(enrOneGeneByAnnot[bestMotif,"geneSet"]), fixed=TRUE)
      )
    })), stringsAsFactors=FALSE)
    tfTable[order(tfTable$NES, decreasing = TRUE),]
  })
  rm(motifEnrichment.asIncidList)
  regulonTargetsInfo <- rbindlist(regulonTargetsInfo)
  
  
  # Optional: Add correlation
  corrMat <- loadInt(scenicOptions, "corrMat", ifNotExists="null")
  if(!is.null(corrMat))
  {
    regulonTargetsInfo$spearCor <- NA_real_
    for(tf in unique(regulonTargetsInfo$TF))
    {
      regulonTargetsInfo[which(regulonTargetsInfo$TF==tf),"spearCor"] <- corrMat[tf, unlist(regulonTargetsInfo[which(regulonTargetsInfo$TF==tf),"gene"])]
    }
  }else warning("It was not possible to add the correlation to the regulonTargetsInfo table.")
  
  
  # Optional: Add Genie3 score
  linkList <- loadInt(scenicOptions, "genie3ll", ifNotExists="null")
  if(!is.null(linkList) & ("weight" %in% colnames(linkList)))
  {
    if(is.data.table(linkList)) linkList <- as.data.frame(linkList)
    
    uniquePairs <- nrow(unique(linkList[,c("TF", "Target")]))
    if(uniquePairs == nrow(linkList)) {
      linkList <- linkList[which(linkList$weight>=getSettings(scenicOptions, "modules/weightThreshold")),]  # TODO: Will not work with GRNBOOST!
      rownames(linkList) <- paste(linkList$TF, linkList$Target,sep="__")
      regulonTargetsInfo <- cbind(regulonTargetsInfo, CoexWeight=linkList[paste(regulonTargetsInfo$TF, regulonTargetsInfo$gene,sep="__"),"weight"])
    }else {
      warning("There are duplicated regulator-target (gene id/name) pairs in the co-expression link list.",
              "\nThe co-expression weight was not added to the regulonTargetsInfo table.")
    }
  }else warning("It was not possible to add the weight to the regulonTargetsInfo table.")
  
  saveRDS(regulonTargetsInfo, file=getIntName(scenicOptions, "regulonTargetsInfo"))
  
  write.table(regulonTargetsInfo, file=getOutName(scenicOptions, "s2_regulonTargetsInfo"),
              sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
  rm(linkList)
  
  # Split into regulons... (output: list TF --> targets)
  regulonTargetsInfo_splitByAnnot <- split(regulonTargetsInfo, regulonTargetsInfo$highConfAnnot)
  regulons <- NULL
  if(!is.null(regulonTargetsInfo_splitByAnnot[["TRUE"]]))
  {
    regulons <- lapply(split(regulonTargetsInfo_splitByAnnot[["TRUE"]], regulonTargetsInfo_splitByAnnot[["TRUE"]][,"TF"]), function(x) sort(as.character(unlist(x[,"gene"]))))
  }
  regulons_extended <- NULL
  if(!is.null(regulonTargetsInfo_splitByAnnot[["FALSE"]]))
  {
    regulons_extended <- lapply(split(regulonTargetsInfo_splitByAnnot[["FALSE"]],regulonTargetsInfo_splitByAnnot[["FALSE"]][,"TF"]), function(x) unname(unlist(x[,"gene"])))
    regulons_extended <- setNames(lapply(names(regulons_extended), function(tf) sort(unique(c(regulons[[tf]], unlist(regulons_extended[[tf]]))))), names(regulons_extended))
    names(regulons_extended) <- paste(names(regulons_extended), "_extended", sep="")
  }
  regulons <- c(regulons, regulons_extended)
  saveRDS(regulons, file=getIntName(scenicOptions, "regulons"))
  
  # Save as incidence matrix (i.e. network)
  incidList <- reshape2::melt(regulons)
  incidMat <- table(incidList[,2], incidList[,1])
  saveRDS(incidMat, file=getIntName(scenicOptions, "regulons_incidMat"))
  rm(incidMat)
  #TODO NMF::aheatmap(incidMat)
  
  if(getSettings(scenicOptions, "verbose")) 
  {
    # Number of regulons and summary of sizes:
    length(regulons) # TODO
    summary(lengths(regulons))
  }
  
  # Finished. Update status.
  scenicOptions@status$current <- 2
  invisible(scenicOptions)
}

################################################################################
# Step 3. Analyzing the network activity in each individual cell
################################################################################
#' @title runSCENIC_3_scoreCells
#' @description Step 3: AUCell (scoring the regulons on the individual cells) 
#' @param scenicOptions Fields used: TODO
#' @param exprMat Expression matrix
#' @param skipBinaryThresholds Whether to skip the automatic binarization step
#' @param skipHeatmap Whether to plot the AUC heatmap
#' @param skipTsne Whether to plot the t-SNE
#' @return The output is written in the folders 'int' and 'ouput'
#' @details See the detailed vignette explaining the internal steps.
#' @examples 
#' runSCENIC_3_scoreCells(scenicOptions)
#' @export
runSCENIC_3_scoreCells <- function(scenicOptions, exprMat, 
                                   skipBinaryThresholds=FALSE, skipHeatmap=FALSE, skipTsne=FALSE)
{
  nCores <- getSettings(scenicOptions, "nCores")
  
  if(is.data.frame(exprMat)) 
  {
    supportedClasses <- paste(gsub("AUCell_buildRankings,", "", methods("AUCell_buildRankings")), collapse=", ")
    supportedClasses <- gsub("-method", "", supportedClasses)
    
    stop("'exprMat' should be one of the following classes: ", supportedClasses, 
         "\n(data.frames are not supported. Please, convert the expression matrix to one of these classes.)")
  }
  
  ################################################################
  ## Prepare regulons
  regulons <- tryCatch(loadInt(scenicOptions, "regulons"),
                       error = function(e) {
                         if(getStatus(scenicOptions, asID=TRUE) < 2) 
                           e$message <- paste0("It seems the regulons have not been built yet. Please, run runSCENIC_2_createRegulons() first.\n", 
                                               e$message)
                         stop(e)
                       })
  regulons <- regulons[order(lengths(regulons), decreasing=TRUE)]
  ### rj debug
  regulons <- regulons[lengths(regulons)>=3]
  # if(length(regulons) <1)  stop("Not enough regulons with at least 10 genes.")
  
  # Add the TF to the regulon (keeping it only once) & rename regulon
  regulons <- setNames(lapply(names(regulons), function(tf) sort(unique(c(gsub("_extended", "", tf), regulons[[tf]])))), names(regulons))
  names(regulons) <- paste(names(regulons), " (",lengths(regulons), "g)", sep="")
  saveRDS(regulons, file=getIntName(scenicOptions, "aucell_regulons"))
  
  msg <- paste0(format(Sys.time(), "%H:%M"), "\tStep 3. Analyzing the network activity in each individual cell")
  if(getSettings(scenicOptions, "verbose")) message(msg)
  
  biggestRegulons <- grep("_extended",names(regulons),invert = T, value = T)
  biggestRegulons <- biggestRegulons[1:min(length(biggestRegulons),10)]
  msg <- paste0("\tNumber of regulons to evaluate on cells: ", length(regulons),
                "\nBiggest (non-extended) regulons: \n",
                paste("\t", biggestRegulons, collapse="\n")) # TODO maxlen?
  if(getSettings(scenicOptions, "verbose")) message(msg)
  
  ################################################################
  # AUCell
  library(AUCell)
  # 1. Create rankings
  set.seed(getSettings(scenicOptions,"seed"))
  tryCatch({
    .openDev(fileName=getIntName(scenicOptions, "aucell_genesStatsPlot"),
             devType=getSettings(scenicOptions, "devType"))
    aucellRankings <- AUCell_buildRankings(exprMat, 
                                           plotStats=TRUE, verbose=getSettings(scenicOptions, "verbose"))
    abline(v=aucellRankings@nGenesDetected["1%"], col="skyblue3", lwd=5, lty=3)
    dev.off()
  },error = function(e) {
    message("Catched error in AUCell_buildRankings() or in the histogram plot: ", e$message)
  })
  saveRDS(aucellRankings, file=getIntName(scenicOptions, "aucell_rankings"))
  
  # 2. Calculate AUC
  regulonAUC <- AUCell_calcAUC(regulons, aucellRankings, 
                               # aucMaxRank=aucellRankings@nGenesDetected["1%"],
                               # ceiling(0.05 * nrow(aucellRankings))
                               nCores=nCores) 
  
  if(nrow(regulonAUC)>=2){
    # Order the modules by similarity, for easier exploration in the upcoming steps & save
    variableRegulons <- names(which(apply(getAUC(regulonAUC), 1, sd) > 0))
    reguDist <- as.dist(1-cor(t(getAUC(regulonAUC)[variableRegulons,]), method="spear"))
    reguClust <- hclust(reguDist, method="ward.D2")
    regulonClusters <- setNames(dynamicTreeCut::cutreeDynamic(reguClust, distM=as.matrix(reguDist), verbose = FALSE), reguClust$labels)
    regulonOrder <- reguClust$labels[reguClust$order]
    regulonOrder <- regulonOrder[order(regulonClusters[regulonOrder], decreasing = TRUE)]
    regulonAUC <- regulonAUC[regulonOrder,]
  }
  
  saveRDS(regulonAUC, file=getIntName(scenicOptions, "aucell_regulonAUC"))
  
  # 3. Default thresholds
  cells_AUCellThresholds <- NULL
  if(!skipBinaryThresholds)
  {
    cells_AUCellThresholds <- AUCell_exploreThresholds(regulonAUC, 
                                                       smallestPopPercent=getSettings(scenicOptions,"aucell/smallestPopPercent"),
                                                       assignCells=TRUE, plotHist=FALSE, 
                                                       verbose=FALSE, nCores=nCores)
    saveRDS(cells_AUCellThresholds, file=getIntName(scenicOptions, "aucell_thresholds"))
    
    # Get cells assigned to each regulon
    regulonsCells <- getAssignments(cells_AUCellThresholds)
    
    ### Save threshold info as text (e.g. to edit/modify...)
    trhAssignment <- getThresholdSelected(cells_AUCellThresholds)
    trhAssignment <- signif(trhAssignment, 3) # TODO why is it sometimes a list? https://github.com/aertslab/AUCell/issues/3
    commentsThresholds <- sapply(cells_AUCellThresholds, function(x) unname(x$aucThr$comment))
    
    table2edit <- cbind(regulon=names(cells_AUCellThresholds),
                        threshold=trhAssignment[names(cells_AUCellThresholds)],
                        nCellsAssigned=lengths(regulonsCells)[names(cells_AUCellThresholds)],
                        AUCellComment=commentsThresholds[names(cells_AUCellThresholds)],
                        nGenes=gsub("[\\(g\\)]", "", regmatches(names(cells_AUCellThresholds), gregexpr("\\(.*?\\)", names(cells_AUCellThresholds)))),
                        clusteringOrder=1:length(cells_AUCellThresholds),
                        # clusterGroup=regulonClusters[names(cells_AUCellThresholds)],
                        clusterGroup = ifelse(nrow(regulonAUC)>=2, regulonClusters[names(cells_AUCellThresholds)], "none"),
                        onlyNonDuplicatedExtended=(names(cells_AUCellThresholds) %in% onlyNonDuplicatedExtended(names(cells_AUCellThresholds))),
                        personalNotes="")
    write.table(table2edit, file=getIntName(scenicOptions, "aucell_thresholdsTxt"), row.names=F, quote=F, sep="\t")
    rm(trhAssignment)
  }
  ####################################################################
  # Plots
  msg <- paste0(format(Sys.time(), "%H:%M"), "\tFinished running AUCell.")
  if(getSettings(scenicOptions, "verbose")) message(msg)
  
  if(!skipHeatmap){
    msg <- paste0(format(Sys.time(), "%H:%M"), "\tPlotting heatmap...")
    if(getSettings(scenicOptions, "verbose")) message(msg)
    
    nCellsHeatmap <- min(500, ncol(regulonAUC))
    cells2plot <- sample(colnames(regulonAUC), nCellsHeatmap)
    
    cellInfo <- loadFile(scenicOptions, getDatasetInfo(scenicOptions, "cellInfo"), ifNotExists="null")   #TODO check if exists, if not... create/ignore?
    if(!is.null(cellInfo)) cellInfo <- data.frame(cellInfo)[cells2plot,,drop=F]
    colVars <- loadFile(scenicOptions, getDatasetInfo(scenicOptions, "colVars"), ifNotExists="null")
    
    fileName <- getOutName(scenicOptions, "s3_AUCheatmap")
    
    fileName <- .openDevHeatmap(fileName=fileName, devType=getSettings(scenicOptions, "devType"))
    NMF::aheatmap(getAUC(regulonAUC)[,cells2plot,drop=FALSE],
                  annCol=cellInfo,
                  annColor=colVars,
                  Rowv = ifelse(nrow(regulonAUC)>=2,TRUE, NA),
                  main="AUC",
                  sub=paste("Subset of",nCellsHeatmap," random cells"),
                  filename=fileName)
    .closeDevHeatmap(devType=getSettings(scenicOptions, "devType"))
  }
  ####################################################################
  # Plots
  if(!skipTsne){
    msg <- paste0(format(Sys.time(), "%H:%M"), "\tPlotting t-SNEs...")
    if(getSettings(scenicOptions, "verbose")) message(msg)
    
    tSNE_fileName <- tsneAUC(scenicOptions, aucType="AUC", onlyHighConf=FALSE) # default: nPcs, perpl, seed, tsne prefix
    tSNE <- readRDS(tSNE_fileName)
    
    # AUCell (activity) plots with the default tsne, as html: 
    fileName <- getOutName(scenicOptions, "s3_AUCtSNE_colAct")
    plotTsne_AUCellHtml(scenicOptions, exprMat, fileName, tSNE) #open the resulting html locally
    
    # Plot cell properties:
    sub <- ""; if("type" %in% names(tSNE)) sub <- paste0("t-SNE on ", tSNE$type)
    cellInfo <- loadFile(scenicOptions, getDatasetInfo(scenicOptions, "cellInfo"), ifNotExists="null") 
    colVars <- loadFile(scenicOptions, getDatasetInfo(scenicOptions, "colVars"), ifNotExists="null")
    pdf(paste0(getOutName(scenicOptions, "s3_AUCtSNE_colProps"),".pdf"))
    plotTsne_cellProps(tSNE$Y, cellInfo=cellInfo, colVars=colVars, cex=1, sub=sub)
    dev.off()
  }
  
  # Finished. Update status.
  scenicOptions@status$current <- 3
  invisible(scenicOptions)
}
#' @title runSCENIC_4_aucell_binarize
#' @description Step 4: Binarize the AUC (and, optional: re-cluster)
#' @param scenicOptions Fields used: TODO
#' @param skipBoxplot Whether to plot the boxplots
#' @param skipHeatmaps Whether to plot the Binary heatmaps
#' @param skipTsne Whether to calculate the binary t-SNE (and optionally, generate a HTML preview)
#' @param exprMat If !skipTsne, the expression matrix can be provided to plot the TF expression on the t-SNE (plotTsne_AUCellHtml)
#' @return The output is written in the folders 'int' and 'ouput'
#' @details See the detailed vignette explaining the internal steps.
#' @examples 
#' runSCENIC_4_aucell_binarize(scenicOptions)
#' @export
runSCENIC_4_aucell_binarize <- function(scenicOptions, 
                                        skipBoxplot=FALSE, 
                                        skipHeatmaps=FALSE, 
                                        skipTsne=FALSE, 
                                        exprMat=NULL)
{
  nCores <- getSettings(scenicOptions, "nCores")
  regulonAUC <- tryCatch(loadInt(scenicOptions, "aucell_regulonAUC"),
                         error = function(e) {
                           if(getStatus(scenicOptions, asID=TRUE) < 3) 
                             e$message <- paste0("It seems the regulons have not been scored on the cells yet. Please, run runSCENIC_3_scoreCells() first.\n", 
                                                 e$message)
                           stop(e)
                         })
  thresholds <- loadInt(scenicOptions, "aucell_thresholds")
  thresholds <- AUCell::getThresholdSelected(thresholds)
  
  # Assign cells
  regulonsCells <- setNames(lapply(names(thresholds), 
                                   function(x) {
                                     trh <- thresholds[x]
                                     names(which(getAUC(regulonAUC)[x,]>trh))
                                   }),names(thresholds))
  ### Convert to matrix (regulons with zero assigned cells are lost)
  regulonActivity <- reshape2::melt(regulonsCells)
  binaryRegulonActivity <- t(table(regulonActivity[,1], regulonActivity[,2]))
  class(binaryRegulonActivity) <- "matrix"
  saveRDS(binaryRegulonActivity, file=getIntName(scenicOptions, "aucell_binary_full"))
  if(nrow(binaryRegulonActivity)==0) stop("No cells passed the binarization.")
  
  # Keep only non-duplicated thresholds
  # (e.g. only "extended" regulons if there is not a regulon based on direct annotation)
  binaryRegulonActivity_nonDupl <- binaryRegulonActivity[which(rownames(binaryRegulonActivity) %in% onlyNonDuplicatedExtended(rownames(binaryRegulonActivity))),,drop=FALSE]
  saveRDS(binaryRegulonActivity_nonDupl, file=getIntName(scenicOptions, "aucell_binary_nonDupl"))
  
  minCells <- ncol(binaryRegulonActivity) * .01
  msg <- paste0("Binary regulon activity: ",
                nrow(binaryRegulonActivity_nonDupl), " TF regulons x ",
                ncol(binaryRegulonActivity), " cells.\n(",
                nrow(binaryRegulonActivity), " regulons including 'extended' versions)\n",
                sum(rowSums(binaryRegulonActivity_nonDupl)>minCells),
                " regulons are active in more than 1% (", minCells, ") cells.")
  if(getSettings(scenicOptions, "verbose")) message(msg)
  
  
  if(!skipBoxplot)
  {
    cat("/////boxplot")
    .openDev(fileName=getOutName(scenicOptions, "s4_boxplotBinaryActivity"),
             devType=getSettings(scenicOptions, "devType"))
    par(mfrow=c(1,2))
    boxplot(rowSums(binaryRegulonActivity_nonDupl), main="nCells per regulon",
            sub='number of cells \nthat have the regulon active',
            col="darkolivegreen1", border="#001100", lwd=2, frame=FALSE)
    boxplot(colSums(binaryRegulonActivity_nonDupl), main="nRegulons per Cell",
            sub='number of regulons \nactive per cell',
            col="darkolivegreen1", border="#001100", lwd=2, frame=FALSE)
    dev.off()
  }
   
  ################################################################################
  # Binary activity heatmap
  if(!skipHeatmaps)
  {
    cat("/////-----prepare\n")
    regulonSelection <- loadInt(scenicOptions, "aucell_regulonSelection", ifNotExists="null", verbose=FALSE)
    if(is.null(regulonSelection)) {
      cat("//////-----regulonselection\n")
      regulonSelection <- regulonSelections(binaryRegulonActivity, binaryRegulonActivity_nonDupl, minCells)
    }

    cat("/////-----1\n")
    cellInfo <- loadFile(scenicOptions, getDatasetInfo(scenicOptions, "cellInfo"), ifNotExists="null")
    cellInfo <- data.frame(cellInfo)
    colVars <- loadFile(scenicOptions, getDatasetInfo(scenicOptions, "colVars"), ifNotExists="null")
    
    
    ### Plot heatmap:
    cat("/////-----2\n")
    for(selRegs in names(regulonSelection$labels))
    {
      cat(selRegs,"\n")
      if(length(regulonSelection[[selRegs]])>0)
      {
        regulonSelection[[selRegs]] <- regulonSelection[[selRegs]][which(regulonSelection[[selRegs]] %in% rownames(binaryRegulonActivity))]
        binaryMat <- binaryRegulonActivity[regulonSelection[[selRegs]],,drop=FALSE]
        
        if(nrow(binaryMat)>0) 
        {
          fileName <- paste0(getOutName(scenicOptions, "s4_binaryActivityHeatmap"),selRegs)
          fileName <- .openDevHeatmap(fileName=fileName, devType=getSettings(scenicOptions, "devType"))
          
          rowv <- ifelse(nrow(binaryMat) >= 2, T, NA)
          colv <- ifelse(ncol(binaryMat) >= 2, T, NA)

          if(nrow(binaryMat)==1){
            NMF::aheatmap(binaryMat, 
                          scale="none", 
                          breaks = c(0,1),
                          revC=TRUE, 
                          main=selRegs,   
                          annCol=cellInfo[colnames(binaryMat),, drop=FALSE],
                          annColor=colVars,
                          Rowv=rowv,
                          Colv=colv,
                          color = c("white", "black"),
                          filename=fileName)

          }else{
            NMF::aheatmap(binaryMat, 
                          scale="none", 
                          breaks = c(0,1),
                          revC=TRUE, 
                          main=selRegs,   
                          annCol=cellInfo[colnames(binaryMat),, drop=FALSE],
                          annColor=colVars,
                          Rowv=rowv,
                          Colv=colv,
                          color = c("white", "black"),
                          filename=fileName)
          }

          if(getSettings(scenicOptions, "devType")!="pdf") dev.off()
        }else{
          if(getSettings(scenicOptions, "verbose")) message(paste0("No regulons to plot for regulon selection '", selRegs, "'. Skipping."))
        }
      }
    }
  }
  
  ################################################################################
  # Tsne - on binary activity
  if(!skipTsne)
  {
    tSNE_fileName <- tsneAUC(scenicOptions, aucType="Binary", filePrefix=getIntName(scenicOptions, "tsne_prefix"), onlyHighConf=FALSE) # default: nPcs, perpl, seed
    if(!is.null(tSNE_fileName))
    {
      tSNE <- readRDS(tSNE_fileName)
      
      # AUCell (activity) as html: 
      fileName <- getOutName(scenicOptions, "s4_binarytSNE_colAct")
      plotTsne_AUCellHtml(scenicOptions, exprMat, fileName, tSNE) #open the resulting html locally
      
      # Plot cell properties:
      sub <- ""; if("type" %in% names(tSNE)) sub <- paste0("t-SNE on ", tSNE$type)
      cellInfo <- loadFile(scenicOptions, getDatasetInfo(scenicOptions, "cellInfo"), ifNotExists="null")
      colVars <- loadFile(scenicOptions, getDatasetInfo(scenicOptions, "colVars"), ifNotExists="null")
      pdf(paste0(getOutName(scenicOptions, "s4_binarytSNE_colProps"),".pdf"))
      plotTsne_cellProps(tSNE$Y, cellInfo=cellInfo, colVars=colVars, cex=1, sub=sub)
      dev.off()
    }
  }
  
  # Finished. Update status.
  scenicOptions@status$current <- 4
  invisible(scenicOptions)
}


################################################################################
# Regulon orders/selection for plots
#' @export
regulonSelections <- function(binaryRegulonActivity, binaryRegulonActivity_nonDupl, minCells)
{
  #binaryRegulonActivity <- loadInt(scenicOptions, "aucell_binary_full")
  #binaryRegulonActivity_nonDupl <- loadInt(scenicOptions, "aucell_binary_nonDupl")
  
  ### Select regulons:
  regulonSelection <- list(labels=c(all="All regulons \n (including duplicated regulons)",
                                    corr="Regulons with any other regulon correlated\n with abs(cor)>0.30 \n(and active in at least 1% of cells)",
                                    onePercent="Regulons active in more than 1% of cells",
                                    notCorr="Regulons with no other regulons correlated\n abs(cor)>0.30 \n or active in fewer than 1% of cells"))
  
  # All regulons.
  regulonSelection[["all"]] <- rownames(binaryRegulonActivity)
  
  # Active in > 1% cells
  regMinCells <- names(which(rowSums(binaryRegulonActivity_nonDupl) > minCells))
  regulonSelection[["onePercent"]] <- regMinCells
  
  # Correlation across regulons (based on binary cell activity)
  reguCor <- cor(t(binaryRegulonActivity_nonDupl[regMinCells,]))
  reguCor[which(is.na(reguCor))] <- 0
  diag(reguCor) <- 0
  
  # Regulons that co-ocurr in similar cells. If a regulon is relevant by itself it will not be shown, also check the regulons ignored.
  corrRegs <- names(which(rowSums(abs(reguCor) > 0.30) > 0))
  regulonSelection[["corr"]]  <- corrRegs
  missingRegs <- rownames(binaryRegulonActivity_nonDupl)[which(!rownames(binaryRegulonActivity_nonDupl) %in% corrRegs)]
  regulonSelection[["notCorr"]]  <- missingRegs
  cat("////save1")
  saveRDS(regulonSelection, file=getIntName(scenicOptions, "aucell_regulonSelection"))
  
  ## Set regulon order (only plotting most correlated regulons)
  reguCor_dist <- as.dist(1-reguCor[corrRegs,corrRegs])
  if(length(reguCor_dist) >= 2) 
  {
    binaryRegulonOrder <- hclust(reguCor_dist)
    binaryRegulonOrder <- binaryRegulonOrder$labels[binaryRegulonOrder$order]
  } else 
  {
    binaryRegulonOrder <- labels(reguCor_dist)
  }
  saveRDS(binaryRegulonOrder, file=getIntName(scenicOptions, "aucell_binaryRegulonOrder"))
  
  return(regulonSelection)
}

.openDev <- function(fileName, devType, ...)
{
  if(devType=="pdf")
    pdf(paste0(fileName, ".pdf"), ...)
  
  if(devType=="png")
    png(paste0(fileName, ".png", type="cairo"), ...)
  
  if(devType=="cairo_pfd") # similar to Cairo::CairoPDF?
    grDevices::cairo_pdf(paste0(fileName, ".pdf"), ...)
}

.openDevHeatmap <- function(fileName, devType)
{
  if(devType!="pdf") 
  {
    if(devType=="png") .openDev(fileName=fileName, devType=devType, width=1200,height=1200)
    if(devType!="png") .openDev(fileName=fileName, devType=devType)
    fileName <- NA
  }else{
    fileName <- paste0(fileName,".pdf")
  }
  return(fileName)
}

.closeDevHeatmap <- function(devType)
{
  if(devType!="pdf") 
  {
    dev.off()
  }
}



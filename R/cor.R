#' GetIterCor
#'
#' @description
#' use this function, obtain correlations under different
#' iterations of different chains
#'
#'
#' @param mamba mamba result object
#' @param murp murp object
#' @param metadata attribute of cell
#' @param iteration_list a list, list(c(1, 2000),c(1,4000))
#' @param adjust TRUE
#'
GetIterCor <- function(object,
                       iteration_list = NULL,
                       chains = NULL,
                       adjust = TRUE){

  # if(missing(iterations_vector) | missing(chains) | missing(T_pred) | missing(cellnum) | missing(truetime) | missing(cluster_label)){
  #   stop("missing something")
  # }
  # if(!is.numeric(truetime)){
  #   truetime <- as.numeric(as.character(truetime))
  # }

  require(dplyr)
  mamba = getPse(object)
  murp = object@MURP
  metadata = object@MetaData

  if(is.null(chains)){
    chains <- mamba@chains
  }
  if(is.null(iteration_list)){
    seq_tmp <- c(seq(0,mamba@iter,2000))[-1]
    iteration_list <- lapply(seq_tmp, function(x) c(1, x))
  }

  svalue <- mamba@sample_value
  parnames <-  mamba@dimnames[[2]]
  iter_chain_cor <- mamba@iter_chain_cor

  cellnum <- nrow(metadata)
  #centers_num <- murpObject$Recommended_K
  cluster_label <- murp$Recommended_K_cl$cluster

  # get T_pred
  T_pred <- lapply(iteration_list, function(iter){
    # cat(iter)
    c = lapply(chains, function(ch){
      if(iter[1]==iter[2]){
        svalue[iter[1]:iter[2], grep(pattern="T\\[", parnames),ch]
      }else{
        apply(svalue[iter[1]:iter[2], grep(pattern="T\\[", parnames),ch], 2, mean)
      }

    })
    names(c) = paste0("chain",chains)
    c
  })
  names(T_pred) <- paste0("iter_", lapply(iteration_list, paste0, collapse="_"))

  # get pseudotime
  pseudotime <- lapply(1:length(iteration_list), function(i){
    tmp = T_pred[[i]]
    pt = lapply(1:length(chains), function(ch){
      GetPseudoTt(CellNum = cellnum, pred_t = tmp[[ch]], Cluster_Label = cluster_label)
    })
    names(pt) = paste0("chain",chains)
    pt
  })
  names(pseudotime) <- paste0("iter_",iteration_list)

  # get cor
  if("truetime" %in% colnames(metadata)){
    truetime <- as.numeric(as.character(metadata$truetime))
    ct_spearman <- matrix(0, nc = length(iteration_list), nr = length(chains))
    colnames(ct_spearman) <- sapply(iteration_list, function(c){paste0("iter_",c[1],"_",c[2])})
    rownames(ct_spearman) <- paste0('chain', chains)
    ct_pearson <- ct_spearman
    for (k in 1:length(iteration_list)){
      for (i in 1:length(chains) ){
        ct_spearman[i, k] <- cor(truetime, pseudotime[[k]][[i]], method = 'spearman')
        ct_pearson[i, k] <- cor(truetime, pseudotime[[k]][[i]], method = 'pearson')
      }
    }

  }else{
    ct_spearman = ct_pearson = NULL
  }


  ## prepare return
  # mi = min(unlist(iteration_list))
  # ma = max(unlist(iteration_list))
  # se = ma/length(iteration_list)
  ## add command
  # command = GetCommand()
  # if(length(iter_chain_cor)==0){
  #   iter_chain_cor[[1]] = list(t_pred = T_pred,
  #                              pseudot = pseudotime,
  #                              spearman = ct_spearman,
  #                              pearson = ct_pearson,
  #                              command = command)
  #   names(iter_chain_cor)[1] =  paste0("iter_", mi, "_", ma, "_seq_", se)
  #
  # }else{
  #   n = length(iter_chain_cor)
  #   iter_chain_cor[[n+1]] = list(t_pred = T_pred,
  #                                pseudot = pseudotime,
  #                                spearman = ct_spearman,
  #                                pearson = ct_pearson,
  #                                command = command)
  #   names(iter_chain_cor)[n+1] <- paste0("iter_", mi, "_", ma, "_seq_", se)
  # }
  iter_chain_cor = list(t_pred = T_pred,
                        pseudot = pseudotime,
                        spearman = ct_spearman,
                        pearson = ct_pearson)

  # 是否要矫正 adjust
  if(adjust){

    # get adjust t_pred
    # T_pred=mamba@iter_chain_cor[[1]]$t_pred
    chain_t = do.call(rbind, tail(T_pred,1)[[1]])
    chain_cor = cor(t(chain_t))
    c = which.max(apply(chain_cor,1,mean))
    # 根据所有的chain的均值得到要反转的chain是哪几个
    if(length(which(chain_cor[c,] > 0)) < (nrow(chain_t)/2) ){
      ch_rev <- which(chain_cor[c,] > 0)
    } else{
      ch_rev <- which(chain_cor[c,] < 0)
    }
    cat(ch_rev,"\n")

    # reverse
    adjust_T_pred <- lapply(1:length(T_pred), function(i){
      chain_t = do.call(rbind, T_pred[[i]])
      chain_t[ch_rev,] <- 1 - chain_t[ch_rev,]
      chain_t
    })
    names(adjust_T_pred) = paste0("iter_", lapply(iteration_list, paste0, collapse="_"))

    # get adjust pseudotime
    adjust_pseudotime <- lapply(1:length(adjust_T_pred), function(i){
      tmp = adjust_T_pred[[i]]
      pt = lapply(1:length(chains), function(ch){
        GetPseudoTt(CellNum = cellnum, pred_t = tmp[ch,], Cluster_Label = cluster_label)
      })
      names(pt) = paste0("chain",chains)
      pt
    })
    names(adjust_pseudotime) <- paste0("iter_", lapply(iteration_list, paste0, collapse="_"))

    # get adjust cor
    adjust_ct_spearman <- matrix(0, nc = length(iteration_list), nr = length(chains))
    colnames(adjust_ct_spearman) <- sapply(iteration_list, function(c){paste0("iter_",c[1],"_",c[2])})
    rownames(adjust_ct_spearman) <- paste0('chain', chains)
    adjust_ct_pearson <- adjust_ct_spearman

    if("truetime" %in% colnames(metadata)){
      for (k in 1:length(iteration_list)){
        for (i in 1:length(chains) ){
          adjust_ct_spearman[i, k] <- cor(truetime, adjust_pseudotime[[k]][[i]], method = 'spearman')
          adjust_ct_pearson[i, k] <- cor(truetime, adjust_pseudotime[[k]][[i]], method = 'pearson')
        }
      }
    }else{
      adjust_ct_spearman = adjust_ct_pearson = NULL
    }


    # prepare
    if(length(ch_rev)==0){
      adjust_chain = NULL
    }else{
      adjust_chain = paste0("chain",ch_rev)
    }

    iter_chain_cor = list(t_pred = T_pred,
                          pseudot = pseudotime,
                          spearman = ct_spearman,
                          pearson = ct_pearson,
                          adjust_chain = adjust_chain,
                          adjust_t_pred = adjust_T_pred,
                          adjust_pseudot = adjust_pseudotime,
                          adjust_spearman = adjust_ct_spearman,
                          adjust_pearson = adjust_ct_pearson )
    # iter_chain_cor = mamba@iter_chain_cor
    # if(length(iter_chain_cor)==0){
    #   iter_chain_cor[[1]] = list(t_pred = T_pred,
    #                              pseudot = pseudotime,
    #                              spearman = ct_spearman,
    #                              pearson = ct_pearson,
    #                              adjust_chain = adjust_chain,
    #                              adjust_t_pred = adjust_T_pred,
    #                              adjust_pseudot = adjust_pseudotime,
    #                              adjust_spearman = adjust_ct_spearman,
    #                              adjust_pearson = adjust_ct_pearson,
    #                              command = command)
    #   names(iter_chain_cor)[1] =  paste0("iter_", mi, "_", ma, "_seq_", se)
    #
    # }else{
    #   n = length(iter_chain_cor)
    #   iter_chain_cor[[n+1]] = list(t_pred = T_pred,
    #                                pseudot = pseudotime,
    #                                spearman = ct_spearman,
    #                                pearson = ct_pearson,
    #                                adjust_chain = adjust_chain,
    #                                adjust_t_pred = adjust_T_pred,
    #                                adjust_pseudot = adjust_pseudotime,
    #                                adjust_spearman = adjust_ct_spearman,
    #                                adjust_pearson = adjust_ct_pearson,
    #                                command = command)
    #   names(iter_chain_cor)[n+1] <- paste0("iter_", mi, "_", ma, "_seq_", se)
    # }

  }

  mamba@iter_chain_cor = iter_chain_cor

  object@OptimResult$pse = mamba
  command = GetCommand()
  object@Command$pse$GetIterCor = command

  return(object)

}

#' GetPredT
#'
#' @description
#' use this function, obtain correlations under different
#' iterations of different chains
#'
#' @param mamba mamba result object
#' @param murp murp object
#' @param metadata attribute of cell
#' @param iteration_list a list, list(c(1, 2000),c(1,4000))
#' @param chains c(1,2,3,4)
#' @param adjust TRUE
#' @param iter_chain_cor_index index of iter_chain_cor from function GetIterCor
#'
GetPredT <- function(object,
                     # iter_chain_cor_index = 1,
                     chains = c(1:10),
                     filter_chain = TRUE,
                     mean_th = 0.45,
                     adjust = TRUE){

  # mamba = ct@OptimResult$pse
  # murp = ct@MURP
  # metadata = ct@MetaData
  # iter_chain_cor_index = 1
  # chains = 1:5
  # adjust = TRUE
  # filter_chain = FALSE
  # mean_th = 0.3

  require(dplyr)
  mamba = getPse(object)
  murp = object@MURP
  metadata = object@MetaData

  t_pred_list = mamba@t_pred
  iter_chain_cor = mamba@iter_chain_cor
  # iter_chain_cor = mamba@iter_chain_cor[[iter_chain_cor_index]]

  # get data from iter_chain_cor
  if(adjust){
    t_pred = iter_chain_cor$adjust_t_pred
    pseudot = iter_chain_cor$adjust_pseudot
    ch_rev = iter_chain_cor$adjust_chain
  }else{
    t_pred = iter_chain_cor$t_pred
    pseudot = iter_chain_cor$pseudot
    ch_rev = NULL
  }
  pseudot_r = lapply(pseudot, function(pse){ do.call(rbind,pse) })
  names(pseudot_r) = names(pseudot)

  # filter chains
  chain_cor = cor(t( tail(pseudot_r,1)[[1]][chains,,drop=FALSE] ))
  if(filter_chain){
    chain_cor_mean = apply(chain_cor, 2, mean)
    filter_chain = rownames(chain_cor)[which(chain_cor_mean > mean_th)]
  }else{
    filter_chain = paste0("chain", chains)
  }
  # if(length(filter_chain)==0){
  #   filter_chain = chains
  # }

  # get mean pseudot
  mean_pseudot_r = lapply(pseudot_r, function(df) apply(df[filter_chain, ,drop=FALSE], 2, mean) )

  # get cor
  if("truetime" %in% colnames(metadata)){
    truetime = as.numeric(as.character(metadata$truetime))
    ct_cor <- matrix(0, nc = length(names(pseudot_r)), nr = 2)
    colnames(ct_cor) <- names(pseudot_r)
    rownames(ct_cor) <- c("spearman", "pearson")
    for (k in 1:length(pseudot_r)){
      ct_cor["spearman", k] <- cor(truetime, mean_pseudot_r[[k]], method = 'spearman')
      ct_cor["pearson", k] <- cor(truetime, mean_pseudot_r[[k]], method = 'pearson')
    }
  }else{
    ct_cor = NULL
  }

  # add command
  # command = GetCommand()
  # if(length(t_pred_list)==0){
  #   t_pred_list[[1]] = list(t_pred = t_pred,
  #                           pseudot = pseudot,
  #                           mean_pseudot = mean_pseudot_r,
  #                           cor = ct_cor,
  #                           filter_chain = filter_chain,
  #                           adjust = adjust,
  #                           adjust_chain = ch_rev,
  #                           command = command)
  #   names(t_pred_list)[1] =  paste0(names(mamba@iter_chain_cor)[iter_chain_cor_index], "_filter_",
  #                                   paste0(filter_chain, collapse = "_"),
  #                                   "_adjust_", adjust)
  # }else{
  #   n = length(t_pred_list)
  #   t_pred_list[[n+1]] = list(t_pred = t_pred,
  #                             pseudot = pseudot,
  #                             mean_pseudot = mean_pseudot_r,
  #                             cor = ct_cor,
  #                             filter_chain = filter_chain,
  #                             adjust = adjust,
  #                             adjust_chain = ch_rev,
  #                             command = command)
  #   names(t_pred_list)[n+1] = paste0(names(mamba@iter_chain_cor)[iter_chain_cor_index], "_filter_",
  #                                    paste0(filter_chain, collapse = "_"),
  #                                    "_adjust_", adjust)
  # }

  t_pred_list = list(t_pred = t_pred,
                     pseudot = pseudot,
                     mean_pseudot = mean_pseudot_r,
                     cor = ct_cor,
                     filter_chain = filter_chain,
                     adjust = adjust,
                     adjust_chain = ch_rev)

  mamba@t_pred = t_pred_list
  object@OptimResult$pse = mamba

  command = GetCommand()
  object@Command$pse$GetPredT = command
  return(object)
}

#' MapMURPLabelToAll
#'
#' @description
#' mapping MURP label to all cells
#'
#' @param vecc
#' @param orig
#'
MapMURPLabelToAll <- function(vecc,
                              orig){

  vecc_names = paste0("T[",1:length(vecc),"]")
  murp_t = data.frame(row.names = vecc_names,
                      names = vecc_names,
                      value = vecc)

  if(is.null(names(orig))){
    orig_names = paste0("c",1:length(orig))
  }else{
    orig_names = names(orig)
  }
  cell_t = data.frame(row.names = orig_names,
                      cell_id = orig_names,
                      names =  paste0("T[",orig,"]"))
  cell_t = merge(murp_t, cell_t, by = "names")
  rownames(cell_t) = cell_t$cell_id
  cell_t = cell_t[orig_names,]

  orig_value = cell_t$value
  names(orig_value) = cell_t$cell_id
  return(orig_value)
}

#' GetPseudoTt
#'
#' @description
#' mapping MURP label to all cells
#'
#' @param CellNum
#' @param pred_t
#' @param Cluster_Label
#' @param CellName
#'
GetPseudoTt <- function(CellNum = NULL,
                        pred_t = NULL,
                        Cluster_Label = NULL,
                        CellName = NULL) {
  names(pred_t) = names( table(Cluster_Label) )
  if(is.null(CellNum)){
    CellNum = length(Cluster_Label)
  }
  Cell_Time = rep(0,CellNum)
  # a <- as.numeric( rownames(as.matrix(table(Cluster_Label))) )
  a <- rownames(as.matrix(table(Cluster_Label)))
  for(i in a) {
    Cell_Time[which(Cluster_Label == i)] = pred_t[as.character(i)]
  }

  if (!is.null(CellName)) {
    names(Cell_Time) <- CellName
  }

  return(Cell_Time)
}

#' GetCorMatrix
#'
#' @description
#' calculate the correlation between two continuous vector
#'
#' @param m matrix, row is sample
#'
# GetCorMatrix <- function(m = NULL){
#   require(clusterCrit)
#   require(aricode)
#   require(infotheo)
#
#   x = matrix(0, nrow(m), nrow(m))
#   for(i in 1:nrow(m)){
#     for(j in 1:nrow(m)){
#       c1 = as.integer(m[i,])
#       c2 =  as.integer(m[j,])
#       # P = extCriteria(c1, c2, crit = "Precision")
#       # R = extCriteria(c1, c2, crit = "Recall")
#       # F1 = 2 * R[[1]] * P[[1]] / (R[[1]] + P[[1]])
#
#       # RI = RI(c1, c2)
#       # x[i,j] = aricode::ARI(c1, c2)
#       x[i,j] = length(which(c1==c2))
#     }
#   }
#   return(x)
# }

#' GetContMI
#'
#' @description
#' calculate the Mutual Information of between two continuous vector
#'
#' @param x a vector
#' @param y a vector
#' @param nB1 5
#' @param nB2 5
#'
# GetContMI <- function(x,
#                       y,
#                       nB1 = 3,
#                       nB2 = 3,
#                       na.rm = TRUE,
#                       ...) {
#   if(na.rm){
#     x = x[which(is.na(x)==FALSE)]
#     y = y[which(is.na(y)==FALSE)]
#   }
#   if(length(x)!=length(y)){
#     return(paste0("Two vectors are not equal in length"))
#   }
#   y2d <- entropy::discretize2d(x, y, numBins1=nB1, numBins2=nB1)
#   mi <- entropy::mi.empirical(y2d)
#
#   return(mi)
# }

#' CorSortTC
#'
#' @Description:
#' collect top chains calculation correlations according
#' converged number of T or C
#'
#'
# GetSortCor <- function(mamba = NULL,
#                        murp = NULL,
#                        diagnostic = NULL,
#                        metadata = NULL,
#                        rm_chains = NULL,
#                        sortindex = 'T'){
#
#   truetime <- as.numeric(as.character(metadata$truetime))
#   cellnum <- nrow(murp$rawdata)
#   cluster_label <- murp$Recommended_K_cl$cluster
#   chains <- mamba@chains
#   chain_num <- length(chains)
#   sortindex_cor <- mamba@sortindex_cor
#
#   svalue <- mamba@sample_value
#   parnames <-  mamba@dimnames[[2]]
#
#   if(is.null(dimnames(svalue))){
#     dimnames(svalue) = mamba@dimnames
#   }
#
#
#   # rm chains
#   diagdf <- diagnostic$converged_param_numbers[2:(chain_num+1),]
#   if(is.vector(diagdf))
#     diagdf = t(as.data.frame(diagdf))
#
#   diagdf2 <- diagdf[setdiff(1:chain_num, rm_chains),]
#   if(is.vector(diagdf2)) {
#     diagdf2 = t(as.data.frame(diagdf2))
#     rownames(diagdf2) = "chain1"
#   }
#
#
#   # create cormat
#   cormat <- matrix(0, ncol = nrow(diagdf2), nrow = 3)
#   colnames(cormat) <- paste0(sortindex, '_top', 1:nrow(diagdf2))
#   rownames(cormat) <- c('chains', "spearman", "pearson")
#
#   # compute cor
#   t_pred <- list()
#   pseudot <- list()
#   index_order <- order(diagdf2[, sortindex], decreasing=TRUE)
#   for(i in 1:nrow(diagdf2)){
#     ch = rownames(diagdf2)[index_order[1:i]]
#     ch = as.numeric(substr(ch,6,nchar(ch)))
#
#     t_pred[[i]] = apply(svalue[, grep(pattern="T\\[", parnames),ch], 2, mean)
#     pseudot[[i]] = GetPseudoTt(CellNum = cellnum, pred_t = t_pred[[i]], Cluster_Label = cluster_label)
#
#     cormat[1, i] = paste(ch, collapse=",")
#     cormat[2, i] = cor(truetime, pseudot[[i]], method = "spearman")
#     cormat[3, i] = cor(truetime, pseudot[[i]], method = "pearson")
#   }
#
#   # add command
#   command = GetCommand()
#   if(length(sortindex_cor)==0){
#     sortindex_cor[[1]] = list(t_pred = t_pred,
#                               pseudot = pseudot,
#                               cormat = cormat,
#                               rm_chains = rm_chains,
#                               command = command)
#     names(sortindex_cor)[1] =  paste0("diagnostic_", diagnostic$diag_name, "_sortindex_", sortindex)
#   }else{
#     n = length(sortindex_cor)
#     sortindex_cor[[n+1]] = list(t_pred = t_pred,
#                                 pseudot = pseudot,
#                                 cormat = cormat,
#                                 rm_chains = rm_chains,
#                                 command = command)
#     names(sortindex_cor)[n+1] <- paste0("diagnostic_", diagnostic$diag_name, "_sortindex_", sortindex)
#   }
#
#   mamba@sortindex_cor = sortindex_cor
#
#   return(mamba)
# }

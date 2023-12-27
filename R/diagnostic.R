
#' GetParamMeanChain
#'
#' @description
#' Get the average value of each chain for any parameter
#'
#' @param mamba mamba object
#' @param murp murp result
#' @param param parameter name
#' @param whichPos 'all' / c(1,2,3)
#' @param iter_range Iteration range, default is c(1, 3000)
#' @param cores threads
#' @export
#'
GetParamMeanChain <- function(mamba = NULL,
                              murp = NULL,
                              param = NULL,
                              whichPos = 'all',
                              iter_range = c(1,3000),
                              cores = 1
){
  # importFrom("ggplot2","ggplot", "ggsave","aes", "labs", "geom_line","facet_wrap")
  # require(reshape)
  # require(stringr)
  # require(pbmcapply)
  # require(doParallel)

  # if(is.null(dimnames(mamba_object$forcast_time@sample_value))){
  #   dimnames(mamba_object$forcast_time@sample_value) <- mamba_object$forcast_time@dimnames
  # }
  chains <- mamba@chains
  K <- length(grep("T",mamba@dimnames[[2]]))
  all_params <- mamba@dimnames[[2]]
  sample_value <- mamba@sample_value
  param_number <- as.data.frame.numeric(table(str_split_fixed(all_params,"\\[",2)[,1]))

  dimnames(sample_value) <- list(paste0("iter",1:dim(sample_value)[1]),
                                 all_params,
                                 paste0("chain",1:dim(sample_value)[3]))

  # check param number
  if(param %in% rownames(param_number) ){

    if(param_number[param,]==1){

      # t = sample_value[,which(all_params==param),]
      t = sample_value[,grep(param,all_params),,drop=FALSE]
      # if(is.vector(t)){
      #   t = as.data.frame(t)
      # }
      df = mclapply(chains, function(ch){
        if(is.numeric(iter_range) & length(iter_range)==2){
          if(iter_range[2]>=iter_range[1] ){
            tmp = t[iter_range[1]:iter_range[2],,ch,drop=F]
          }
        }else{
          tmp = t[,ch,drop=F]
        }
        mean(tmp)
      }, mc.cores = cores)
      df = do.call(rbind, df)
      rownames(df) = paste0("chain", chains)

    }else{
      if(whichPos == 'all'){
        t = sample_value[,grep(pattern=paste0(param, "\\["),all_params),,drop=FALSE]
        # if(is.vector(t)){
        #   t = as.data.frame(t)
        # }
        dimnames(t)[[2]] = all_params[grep(pattern=paste0(param, "\\["),all_params)]
      }else{
        if(param == "rho" | param == "Y"){
          index_name = unlist(lapply(whichPos, function(i) paste0(param,"[",i[1],", ",i[2],"]")) )
          t = sample_value[, index_name, drop=F]
          dimnames(t)[[2]] = index_name
        }else if(param == "pi"){
          index_name = unlist(lapply(whichPos, function(i) paste0(param,"[",i[1],", ",i[2],", ",i[3],"]")) )
          t = sample_value[, index_name, drop=F]
          dimnames(t)[[2]] = index_name
        }else{
          t = sample_value[, paste0(param, "[",whichPos,"]"), ]
          dimnames(t)[[2]] = paste0(param, "[",whichPos,"]")
        }
      }

      # compute mean value
      df <- mclapply(chains, function(ch){
        if(is.numeric(iter_range) & length(iter_range)==2){
          if(iter_range[2]>=iter_range[1] ){
            if(length(chains)==1 & length(dim(t))==2){
              tmp = t[iter_range[1]:iter_range[2],,drop=F]
            }else{
              tmp = t[iter_range[1]:iter_range[2],,ch,drop=F]
            }
          }
        }else{
          tmp = t[,,ch,drop=F]
        }
        if(is.vector(tmp)){
          tmp
        }else{
          apply(tmp, 2, mean)
        }

      }, mc.cores = cores)
      df = do.call(rbind, df)
      rownames(df) = paste0("chain", chains)

    }
  }
  return(df)
}

#' GetAllParamMeanChain
#'
#' @description
#' Get the average value of all parameters on each chain
#'

#' @param object mgpfact object
#' @param iter_range Iteration range, default is c(1, 3000)
#' @param cores threads
#' @param aspect pse or track sub object
#' @param save whether save to object
#'
#' @export
#'
GetAllParamMeanChain <- function(object,
                                 iter_range = NULL,
                                 cores = 1,
                                 aspect = "pse",
                                 save = F){

  # aspect = "pse"
  # iter_range = rep(getParams(ct, "pse_optim_iterations"),2)
  # whichPos = 'all'
  # iter_range = iter_range
  # cores = 1
  # param = "m_t"

  mamba = object@OptimResult[[aspect]]
  murp = object@MURP
  all_params <- mamba@dimnames[[2]]
  param_names <- rownames(as.data.frame.numeric(table(str_split_fixed(all_params,"\\[",2)[,1])))
  param_number <- as.data.frame.numeric(table(str_split_fixed(all_params,"\\[",2)[,1]))

  param_meanvalue <- lapply(param_names, function(param){
    print(param)
    GetParamMeanChain(mamba = mamba,
                      murp = murp,
                      param = param,
                      whichPos = 'all',
                      iter_range = iter_range,
                      cores = 1)
  })
  names(param_meanvalue) <- param_names

  if(save){
    if(aspect=="pse"){
      save(param_meanvalue,
           file = paste0("2_pseudotime/param_meanvalue_", iter_range[1],"_",iter_range[2], ".rda"))
    }
    if(aspect=="track"){
      save(param_meanvalue,
           file = paste0("3_tracking/param_meanvalue_", iter_range[1],"_",iter_range[2], ".rda"))
    }
  }
  return(param_meanvalue)
}

#' PlotParamSampling
#'
#' @description
#' Get the average value of all parameters on each chain
#'
#' @param object mgpfact object
#' @param iter_range Iteration range, default is c(1, 3000)
#' @param params parameter name, can be a character vector.
#' @export
#'
PlotParamSampling <- function(object,
                              iter_range = NULL,
                              params = "T"){

  mamba = getPse(object)
  murp = object@MURP
  # iter_range = rep(getParams(ct, "pse_optim_iterations"),2)
  # param = "m_t"

  if(is.null(iter_range)){
    iter_range = c(1, mamba@iter)
  }

  all_params <- mamba@dimnames[[2]]
  param_names <- rownames(as.data.frame.numeric(table(str_split_fixed(all_params,"\\[",2)[,1])))
  param_number <- as.data.frame.numeric(table(str_split_fixed(all_params,"\\[",2)[,1]))
  colnames(param_number) <- "number"
  params_names = mamba@dimnames[[2]]

  gtag = paste0(params,collapse="|")
  gtag = gsub("\\[","\\\\[",gtag)
  gtag = gsub("\\]","\\\\]",gtag)
  pind = grep(gtag, params_names)
  pnames = params_names[pind]

  tmp = mamba@sample_value[iter_range[1]:iter_range[2],pind,,drop=F]
  dimnames(tmp) = list(paste0(1:dim(tmp)[1]), pnames, paste0("chains ",1:dim(tmp)[3]))
  df = reshape2::melt(tmp)
  colnames(df) = c("iter", "param", "chains", "value")

  ggplot(df,aes(x = .data$iter, y = .data$value, color = .data$chains)) +
    # geom_point(size = 0.5, shape = 21) +
    geom_line(size = 0.3, alpha = 0.8) +
    scale_color_simpsons() +
    facet_wrap( ~ param, ncol = 5, scales = "free") +
    labs(y = "Sampling Value") +
    rj.ftheme -> p
  ge = ifelse(length(params)>5,5,length(params))
  gef = paste0(params[1:ge], collapse = "_")
  if(ge<5){
    ggsave(paste0("2_pseudotime/2.1_julia_result/params_",gef,".pdf"),
           p, width = ge*4+0.57, height = 4)
  }else{
    ggsave(paste0("2_pseudotime/2.1_julia_result/params_",gef,".pdf"),
           p, width = 20, height = ceiling(length(params)/5)*4)
  }
}


#' PlotParamLogpdf
#'
#' @description
#' Get the average value of all parameters on each chain
#'
#' @param object mgpfact object
#' @param iter_range Iteration range, default is c(1, 3000)
#'
#' @export
#'
PlotParamLogpdf <- function(object = NULL,
                            iter_range = NULL){

  mamba = getPse(object)
  murp = object@MURP
  logpdf = mamba@logpdf$logpdf
  if(is.null(iter_range)){
    iter_range = c(1, mamba@iter)
  }
  ## params
  tmp = map_df(seq_along(logpdf), function(i){
    x = data.frame(logpdf[[i]][iter_range[1]:iter_range[2],,drop=F])
    dimnames(x) = list(1:nrow(x), paste0("chain ",1:ncol(x)))
    x$params = names(logpdf)[i]
    x$iter = 1:nrow(x)
    x
  })
  df = melt(tmp, id.vars = c("params","iter"))

  ggplot(df,aes(x = .data$iter, y = .data$value, color = .data$variable)) +
    # geom_point(size = 0.5, shape = 21) +
    geom_line(size = 0.3, alpha = 0.8) +
    scale_color_simpsons() +
    facet_wrap( ~ params, ncol = 6, scales = "free") +
    labs(x = "Iterations", y = "logpdf") +
    rj.ftheme -> p
  ggsave(paste0("2_pseudotime/2.1_julia_result/params_logdf_params.pdf"),
         p, width = 24, height = 8)

  ## l_all
  tmp = data.frame(logpdf$l_all)
  dimnames(tmp) = list(1:nrow(tmp), paste0("chain ",1:ncol(tmp)))
  tmp$iter = 1:nrow(tmp)
  df = reshape2::melt(tmp, id.vars = c("iter"))
  ggplot(df,aes(x = .data$iter, y = .data$value, color = .data$variable)) +
    # geom_point(size = 0.5, shape = 21) +
    geom_line(size = 0.7, alpha = 0.8) +
    scale_color_simpsons() +
    labs(x = "Iterations", y = "logpdf") +
    rj.ftheme -> p
  ggsave(paste0("2_pseudotime/2.1_julia_result/params_logdf_all.pdf"),
         p, width = 7, height = 5)
}

#' CovHeatmap
#'
#' @description:
#' get the covariance matrix
#'
#' @param object MGPfact object
#' @export
#'
PlotCovHeatmap <- function(object){

  pl = list()
  L = object@Settings@settings$trajectory_number
  ord = order(GetMURPInfo(object)$T)
  s = object@GPR$murp_cov$cov
  p = Heatmap(s[ord,ord],
              use_raster = TRUE,
              # col = c("#00468bff","black", "yellow"),
              col = c("yellow", "grey10"),
              name = paste0("cov value"),
              cluster_columns = F,
              cluster_rows = F,
              show_row_names = FALSE,
              show_column_names = FALSE,
              row_title = NULL,
              column_title = paste0("All"),
              show_heatmap_legend = T)
  pl = append(pl, list(p))

  s_l = object@GPR$murp_cov$cov_l
  for(l in 1:L){
    s = s_l[[l]]
    p = Heatmap(s[ord,ord],
                use_raster = TRUE,
                # col = c("#00468bff","black", "yellow"),
                col = c("yellow", "grey10"),
                name = paste0("cov value"),
                cluster_columns = F,
                cluster_rows = F,
                show_row_names = FALSE,
                show_column_names = FALSE,
                row_title = NULL,
                column_title = paste0("Trajectory ", l),
                show_heatmap_legend = T)
    pl = append(pl, list(p))
  }

  p = pl[[1]] + pl[[2]] + pl[[3]] + pl[[4]]
  pdf(paste0("2_pseudotime/2.4_cov_matrix/cov_heatmap.pdf"),
      width = 12, height = 3.2)
  ht = draw(p,
            legend_grouping = "original",
            gap = unit(8, "mm"),
            merge_legends = FALSE,
            legend_gap = unit(10, "mm"),
            heatmap_legend_side = "right")
  dev.off()
}

#------------
#' GetParamIter
#'
#' @description
#' Mamba parameter iteration result
#'
#' @param mamba mamba object
#' @param murp murp result
#' @param param parameter name
#' @param whichPos 'all' / c(1,2,3)
#' @param cores threads
#'
#' @export
# GetParamIter <- function(mamba = NULL,
#                          murp = NULL,
#                          param = NULL,
#                          whichPos = 'all',
#                          cores = 2
# ){
#
#   require(reshape)
#   require(stringr)
#   require(pbmcapply)
#   require(doParallel)
#   # if(is.null(dimnames(mamba_object$forcast_time@sample_value))){
#   #   dimnames(mamba_object$forcast_time@sample_value) <- mamba_object$forcast_time@dimnames
#   # }
#
#   K <- murp$Recommended_K
#   all_params <- mamba@dimnames[[2]]
#   sample_value <- mamba@sample_value
#   param_number <- as.data.frame.numeric(table(str_split_fixed(all_params,"\\[",2)[,1]))
#   L <- sqrt(param_number['R',])
#
#   dimnames(sample_value) <- list(paste0("iter", 1:dim(sample_value)[1]),
#                                  all_params,
#                                  paste0("chain", 1:dim(sample_value)[3]))
#
#   # check param number
#   if(param %in% rownames(param_number) ){
#
#     if(param_number[param,]==1){
#
#       t = sample_value[,which(all_params==param),]
#       df = reshape::melt(t)
#       df = data.frame(apply(df,2,as.character),
#                       stringsAsFactors = FALSE)
#       colnames(df) <- c('iter','chain','value')
#       df$iter <- unlist(pbmclapply(df$iter,
#                                    function(i){substr(x = i,start = 5,stop = nchar(i))},
#                                    mc.cores = cores))
#       df$param <- as.factor(param)
#     }else{
#       if(whichPos == 'all'){
#         t = sample_value[,grep(pattern=paste0(param, "\\["),all_params),]
#       }else{
#         if(param == "Y" | param == "Z" | param == "pi" | param == "R" | param == "C"){
#           index_name = unlist(lapply(whichPos, function(i) paste0(param,"[",i[1],", ",i[2],"]")) )
#           t = sample_value[, index_name, ]
#         }else if(param == "xx"){
#           #index_name = unlist(lapply(whichPos, function(i) paste0(param,"[",i[1],", ",i[2],"]")) )
#           index_name = unlist(lapply(whichPos, function(i) paste0(param,"[",i[1],", ",i[2],", ",i[3],"]")) )
#           t = sample_value[, index_name, ]
#         }else{
#           t = sample_value[, paste0(param, "[",whichPos,"]"), ]
#         }
#       }
#
#       # create df
#       df <- reshape::melt(t)
#       df <- data.frame(apply(df,2,as.character),
#                        stringsAsFactors = FALSE)
#       colnames(df) <- c('iter','param','chain','value')
#       df$iter <- unlist(mclapply(df$iter,
#                                  function(i){substr(x = i,start = 5,stop = nchar(i))},
#                                  mc.cores = cores))
#     }
#   }
#   # set attributes
#   if(param %in% rownames(param_number)){
#     if(param_number[param,]!=1){
#
#       if(param == "Y" | param == 'Z' | param == "pi" | param == 'C'){                   #---Y
#         if(whichPos=='all'){
#           level <- paste0(param, "[", 1:(length(table(df$param))/K) )
#           level <- unlist(lapply(level, function(a) paste0(a, ", ", 1:K, "]")) )
#           df$param <- factor(df$param, levels = level)
#         }else{
#           df$param <- factor(df$param)
#         }
#       }else if(param == 'R'){             #---R
#         if(whichPos=='all'){
#           level <- paste0(param, "[", 1:(length(table(df$param))/L) )
#           level <- unlist(lapply(level, function(a) paste0(a, ", ", 1:L, "]")) )
#           df$param <- factor(df$param, levels = level)
#         }else{
#           df$param <- factor(df$param)
#         }
#       }else if(param == "xx"){            #---xx
#         if(whichPos=='all'){
#           level <- paste0(param, "[",1:(length(table(df$param))/2) )
#           level <- unlist(lapply(level, function(a) paste0(a, ", ", 1:length(grep(pattern="T\\[",all_params)), ", ")) )
#           level <- unlist(lapply(level, function(a) paste0(a, 1:2, "]")) )
#           df$param <- factor(df$param, levels = level)
#         }else{
#           df$param <- factor(df$param)
#         }
#       }else{
#         if(whichPos=='all'){
#           df$param <- factor(df$param, levels = paste0(param, "[",1:length(table(df$param)),"]"))
#         }else{
#           df$param <- factor(df$param, levels = paste0(param, "[",whichPos,"]"))
#         }
#       }
#
#     }
#   }
#   df$chain <- factor(df$chain, levels = paste0("chain",1:dim(sample_value)[3]))
#
#   return(df)
# }

#' GetParamMeanIter
#'
#' @description
#' Get the mean value of the parameters under the corresponding iteration
#'
#' @param mamba mamba object
#' @param murp murp result
#' @param param parameter name
#' @param whichPos 'all' / c(1,2,3)
#' @param cores threads
#'
#' @export
# GetParamMeanIter <- function(mamba = NULL,
#                              murp = NULL,
#                              param = NULL,
#                              whichPos = 'all',
#                              cores = 2
# ){
#
#   require(reshape)
#   require(stringr)
#   require(pbmcapply)
#   require(doParallel)
#
#   # if(is.null(dimnames(mamba_object$forcast_time@sample_value))){
#   #   dimnames(mamba_object$forcast_time@sample_value) <- mamba_object$forcast_time@dimnames
#   # }
#
#   K <- murp$Recommended_K
#   all_params <- mamba@dimnames[[2]]
#   sample_value <- mamba@sample_value
#   param_number <- as.data.frame.numeric(table(str_split_fixed(all_params,"\\[",2)[,1]))
#   L <- sqrt(param_number['R',])
#
#   dimnames(sample_value) <- list(paste0("iter", 1:dim(sample_value)[1]),
#                                  all_params,
#                                  paste0("chain", 1:dim(sample_value)[3]))
#
#   # check param number
#   if(param %in% rownames(param_number) ){
#
#     if(param_number[param,]==1){
#
#       t = sample_value[,which(all_params==param),]
#       t2 = t
#       for(i in 1:dim(t2)[1]){
#         for(j in 1:dim(t2)[2]){
#           t2[i, j] =  mean(t2[1:i, j])
#         }
#       }
#       df = reshape::melt(t2)
#
#       df = data.frame(apply(df,2,as.character),
#                       stringsAsFactors = FALSE)
#       colnames(df) <- c('iter','chain','value')
#       df$iter <- unlist(pbmclapply(df$iter,
#                                    function(i){substr(x = i,start = 5,stop = nchar(i))},
#                                    mc.cores = cores))
#       df$param <- as.factor(param)
#     }else{
#       if(whichPos == 'all'){
#         t = sample_value[,grep(pattern=paste0(param, "\\["),all_params),,drop = FALSE]
#       }else{
#         if(param == "Y" | param == 'Z' | param == "pi" | param == 'R' | param == 'C'){
#           #index_name = unlist(lapply(whichPos, function(i) paste0(param,"[",i[1],", ",i[2],", ",i[3],"]")) )
#           index_name = unlist(lapply(whichPos, function(i) paste0(param,"[",i[1],", ",i[2],"]")) )
#           t = sample_value[, index_name, ]
#         }else{
#           t = sample_value[, paste0(param, "[",whichPos,"]"), ]
#         }
#       }
#
#       # create df
#       t2 <- t
#       for(i in 1:dim(t2)[1]){
#         for(j in 1:dim(t2)[2]){
#           for(k in 1:dim(t2)[3]){
#             t2[i, j, k] = mean(t2[1:i, j, k])
#           }
#         }
#       }
#       df <- reshape::melt(t2)
#       df <- data.frame(apply(df,2,as.character),
#                        stringsAsFactors = FALSE)
#       colnames(df) <- c('iter','param','chain','value')
#       df$iter <- unlist(mclapply(df$iter,
#                                  function(i){substr(x = i,start = 5,stop = nchar(i))},
#                                  mc.cores = cores))
#
#     }
#   }
#
#   # set attributes
#   if(param %in% rownames(param_number)){
#     if(param_number[param,]!=1){
#
#       if(param == "Y" | param == 'Z' | param == "pi" | param == 'C'){                   #---Y
#         if(whichPos=='all'){
#           level <- paste0(param, "[", 1:(length(table(df$param))/K) )
#           level <- unlist(lapply(level, function(a) paste0(a, ", ", 1:K, "]")) )
#           df$param <- factor(df$param, levels = level)
#         }else{
#           df$param <- factor(df$param)
#         }
#       }else if(param == 'R'){             #---R
#         if(whichPos=='all'){
#           level <- paste0(param, "[", 1:(length(table(df$param))/L) )
#           level <- unlist(lapply(level, function(a) paste0(a, ", ", 1:L, "]")) )
#           df$param <- factor(df$param, levels = level)
#         }else{
#           df$param <- factor(df$param)
#         }
#       }else{
#         if(whichPos=='all'){
#           df$param <- factor(df$param, levels = paste0(param, "[",1:length(table(df$param)),"]"))
#         }else{
#           df$param <- factor(df$param, levels = paste0(param, "[",whichPos,"]"))
#         }
#       }
#       # end
#     }
#   }
#   df$chain <- factor(df$chain, levels = paste0("chain",1:dim(sample_value)[3]))
#
#   return(df)
# }

#' GetGewDiag
#'
#' @description
#' use this function with julia, obtain the number of convergences
#' for each parameter under different iterations of different chains
#'
#' @param mamba mamba result object
#' @param gvalue gew_value
#' @param cnames gew_colnames
#' @param rnames gew_rownames
#' @param chains chains
#' @param diagnosticMethods "gelmandiag", "gewekediag",
#' "heideldiag", "rafterydiag"
#' @export
# GetGewDiag <- function(mamba = NULL,
#                        gvalue = NULL,
#                        cnames = NULL,
#                        rnames = NULL,
#                        diagnosticMethods = NULL){
#
#   require(stringr)
#
#   chains <- mamba@chains
#   diag_object <- mamba@diagnostic
#
#   if(length(diag_object) > 0){
#     if(length(grep(pattern = diagnosticMethods, names(diag_object))) > 0 ){
#       stop("You have already diagnosed in this way.")
#     }
#   }
#
#   diag <- lapply(1:dim(gvalue)[3], function(x){
#     tmp = gvalue[,,x]
#     rownames(tmp) = rnames
#     colnames(tmp) = cnames
#     tmp
#   })
#
#   cs_number <- t(as.data.frame.numeric(table(str_split_fixed(rnames,"\\[",2)[,1])))
#
#   cs <- matrix(0, ncol = ncol(cs_number), nrow = length(chains) + 1)
#   colnames(cs) <- colnames(cs_number)
#   cs[1,] <- cs_number[1,]
#   rownames(cs) <- c('pal_number', paste0('chain',chains))
#
#   for(i in chains){
#     chain = paste0("chain",i)
#     for (j in colnames(cs)){
#       if(cs[1, j]==1){
#         cs[chain, j] <- length(which(diag[[i]][which(rnames==j), 2] > 0.05))
#       }else{
#         cs[chain, j] <- length(which(diag[[i]][grep(j,rnames), 2] > 0.05))
#       }
#     }
#     print("next")
#   }
#
#   # set
#   if(length(diag_object)==0){
#     diag_object[[1]] = list(diag_frame = diag, converged_param_numbers = cs,diag_name = diagnosticMethods)
#     names(diag_object)[1] =  diagnosticMethods
#   }else{
#     n = length(diag_object)
#     diag_object[[n+1]] = list(diag_frame = diag, converged_param_numbers = cs,diag_name = diagnosticMethods)
#     names(diag_object)[n+1] = diagnosticMethods
#   }
#
#   mamba@diagnostic = diag_object
#
#   return(mamba)
# }

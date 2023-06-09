
#' GetMURPMapLabel
#' @description
#' map the existing attribute labels about cells to murp
#'
#' Arguments:
#' @param object celltrek object
#' @param labels
#'
#' @export
#'
GetMURPMapLabel <- function(object, labels = NULL){

  ind = which(labels %in% colnames(object@MetaData)) %>% length
  if(ind==0) labels = NULL
  len = length(labels)

  ######### 1. get murp_cellinfo
  if(!is.null(labels)){
    tmp <- lapply(labels, function(lab){
      cat(lab, "\n")
      tab = table(object@MetaData[,c(lab,"murp_cluster")]) %>% as.matrix
      tab_prop = apply(tab, 2, function(x){ x/sum(x) })
      if(is.null(dim(tab_prop))) {
        tmp = rep(rownames(tab),length(tab_prop))
        names(tmp) = names(tab_prop)
        tmp
      }else{
        apply(tab_prop, 2, function(x){rownames(tab_prop)[which.max(x)]})
      }

    })
    tmp <- do.call(cbind, tmp) %>% as.data.frame
    colnames(tmp) <- d2_label
    object = AddMURPMetadata(object, tmp)
    # murp_cellinfo <- data.frame(row.names = paste0("T[",1:nrow(object@MURP$Recommended_K_cl$centers),"]"),
    #                             tmp)
    # ori = GetMURPInfo(object)
    # if(is.null(ori)){
    #   object@MURP$murp_cellinfo <- murp_cellinfo
    # }else{
    #
    # }

  }
  object = assignSettings(object,"label",labels)
  return(object)
}

#' GetMURPInfo
#'
#' @param object
#'
GetMURPInfo <- function(object){
  r = object@MURP$murp_cellinfo
  return(r)
}

#' AddMURPMetadata
#'
#' @param object
#' @param df
#'
AddMURPMetadata <- function(object, df){

  meta = GetMURPInfo(object)
  if(is.null(meta)){
    meta = data.frame(row.names = paste0("T[",1:nrow(object@MURP$Recommended_K_cl$centers),"]"),
                      murp_cluster = 1:object@MURP$Recommended_K)
  }
  i = which(colnames(meta) %nin% colnames(df))
  cols = colnames(meta)[i]
  meta = cbind(meta[,cols,drop=F], df)
  object@MURP$murp_cellinfo = meta
  return(object)
}

#' AddMURPMetadata
#'
#' @param object
#' @param df
#'
AddMetadata <- function(object, df){

  meta = object@MetaData
  i = which(colnames(meta) %nin% colnames(df))
  cols = colnames(meta)[i]
  meta = cbind(meta[,cols,drop=F], df)

  object@MetaData = meta
  return(object)
}

#' GetPseSdf
#'
#' @description
#' get all param info from optim result
#'
#' @param object
#' @param param_meanvalue
#' @param unified_direction
#' @param rm_adjust_chain
#' @param rev
#'
GetPseSdf <- function(object,
                      param_meanvalue = NULL,
                      unified_direction = FALSE,
                      rm_adjust_chain = FALSE,
                      rev = FALSE){

  # param_meanvalue = GetAllParamMeanChain(object = ct, aspect = "pse",
  #                                        iter_range = rep(getParams(ct, "pse_optim_iterations"),2))
  # object = ct
  # unified_direction = FALSE
  # rm_adjust_chain = TRUE
  # rev = FALSE
  # unified_direction = TRUE
  # rm_adjust_chain = FALSE
  # rev = FALSE

  if(is.null(param_meanvalue)){
    iter = getParams(ct,"pse_optim_iterations")
    # load(paste0("2_pseudotime/param_meanvalue_",iter,"_",iter,".rda"))
    param_meanvalue <- GetAllParamMeanChain(object = object, aspect = "pse",
                                            iter_range = rep(getParams(object, "pse_optim_iterations"),2))
  }
  if(!is.null(labels)){
    labels = getParams(object,"label")
  }

  mamba = object@OptimResult$pse
  murp = object@MURP
  metadata = object@MetaData
  chains = mamba@chains

  P <- getParams(object, "murp_number")
  L <- getParams(object, "trajectory_number")
  Q <- getParams(object, "murp_pc_number")

  ## 哪些chain是最终保留下来的
  t_pred =mamba@t_pred
  if(object@Settings@settings$chains_number>1){
    if(t_pred$adjust){
      rev_ch <- t_pred$adjust_chain
    }else{
      rev_ch <- NULL
    }
    if(rm_adjust_chain){
      filter_ch <- setdiff(t_pred$keep_chain,rev_ch)
    }else{
      filter_ch <- t_pred$keep_chain
    }
  }else{
    filter_ch <- t_pred$keep_chain
  }
  cat("filter_chain", filter_ch, "\n")

  param_tb <- param_meanvalue$Tb[filter_ch,,drop=FALSE]
  param_T <- param_meanvalue$T[filter_ch,,drop=FALSE]
  param_X <- param_meanvalue$X[filter_ch,,drop=FALSE]
  param_lambda1 <- param_meanvalue$lambda1[filter_ch,, drop=FALSE]
  param_lambda2 <- param_meanvalue$lambda2[filter_ch,, drop=FALSE]
  param_lambda3 <- param_meanvalue$lambda3[filter_ch,, drop=FALSE]
  param_mt <- param_meanvalue$m_t[filter_ch,, drop=FALSE] %>% mean
  param_s2t <- param_meanvalue$s2_t[filter_ch,, drop=FALSE] %>% mean
  param_s2x <- param_meanvalue$s2_x[filter_ch,, drop=FALSE] %>% mean
  param_rho <- param_meanvalue$rho[filter_ch,, drop=FALSE] %>% mean

  ### 用P确定C
  x <- param_meanvalue$C[filter_ch, , drop=FALSE]
  param_p <- param_meanvalue$p[filter_ch, , drop=FALSE]
  param_C <- ifelse(param_p>0.5, 1, 2)
  colnames(param_C) <- colnames(x)

  ### 方向的统一性，需要同时翻转T和Tb
  if(unified_direction & !rm_adjust_chain){
    if(length(which(rev_ch %in% filter_ch))!=0 ){
      inter = intersect(rev_ch, filter_ch)
      param_tb[inter, ] = 1 - param_tb[inter, ]
      param_T[inter, ] = 1 - param_T[inter, ]
    }
  }

  ### 方向翻转
  if(rev){
    param_T <- 1- param_T
    param_tb <- 1 - param_tb
  }

  ### 取均值
  T <- apply(param_T, 2, mean)
  X <- apply(param_X, 2, mean)

  ## 判断是不是单链
  if(nrow(param_tb)==1){

    tb_list = as.list( param_tb )
    l1_list = param_lambda1[1,,drop=F]
    l2_list = param_lambda2[1,,drop=F]
    l3_list = param_lambda3[1,,drop=F]
    tmp = matrix(param_C[1,], nr=L)

    c_combine = tmp
    rownames(c_combine) = paste0("L_", 1:L)

    c0_combine <- c_combine
    for(i in 1:L){
      ind = which(T < mean(tb_list[[i]]))
      if(length(ind)!=0)
        c0_combine[i,ind] = 0
    }

  }else{

    ### 1. 对10个chain的Tb向量化并聚类，得到均值Tb
    param_tb_2 <- param_tb
    colnames(param_tb_2) <- paste0("L",1:L)
    tb_chain <- as.vector(param_tb_2)
    names(tb_chain) <- unlist(lapply(colnames(param_tb_2),
                                     function(x) paste0(rownames(param_tb_2),"_",x)))
    tb_pam <- cluster::pam(tb_chain, L, metric = "manhattan")
    pdf("2_pseudotime/2.2_binarytree/tb_cluster.pdf")
    plot(x = tb_chain, y = tb_chain, col = tb_pam$clustering)
    dev.off()
    # Tb <- tapply(tb_chain, tb_pam$clustering, mean)
    # names(Tb) <- paste0("Tb", 1:L)

    ### 2. 准备一个漂亮的C
    C_list <- lapply(filter_ch, function(ch){
      tmp = matrix(param_C[ch, ], nr=L)
      rownames(tmp) = c(paste0(ch,"_L",1:L))
      tmp
    })
    names(C_list) <- c(filter_ch)
    c <- do.call(rbind, C_list)

    ### 3. 根据聚类对C | Tb | lambda1 | lambda2 | lambda3分组
    cluster <- tb_pam$clustering
    split_c <- split(names(cluster), as.factor(cluster))
    group_list <- lapply(split_c, function(x){
      if(length(x)==1){
        tmp = matrix(c[x,],nr=1)
        dimnames(tmp) = list(x,colnames(c))
        tmp
      }else{
        c[x,]
      }
    })

    param_tb_2 <- param_tb
    colnames(param_tb_2) <- paste0("L",1:L)
    tb <- melt(param_tb_2)
    tb <- data.frame(row.names = paste0(tb[,1],"_",tb[,2]), tb = tb[,3])
    tb_list <- split(tb[names(cluster),], as.factor(cluster))

    param_l1_2 <- param_lambda1
    colnames(param_l1_2) <- paste0("L",1:L)
    l1 <- melt(param_l1_2)
    l1 <- data.frame(row.names = paste0(l1[,1],"_",l1[,2]), l1 = l1[,3])
    l1_list <- split(l1[names(cluster),], as.factor(cluster))

    param_l2_2 <- param_lambda2
    colnames(param_l2_2) <- paste0("L",1:L)
    l2 <- melt(param_l2_2)
    l2 <- data.frame(row.names = paste0(l2[,1],"_",l2[,2]), l2 = l2[,3])
    l2_list <- split(l2[names(cluster),], as.factor(cluster))

    param_l3_2 <- param_lambda3
    colnames(param_l3_2) <- paste0("L",1:L)
    l3 <- melt(param_l3_2)
    l3 <- data.frame(row.names = paste0(l3[,1],"_",l3[,2]), l3 = l3[,3])
    l3_list <- split(l3[names(cluster),], as.factor(cluster))

    ### 5. 翻转
    group_rev_list <- lapply(group_list, function(xt){
      xi = as.matrix(xt)
      x_ind = which(apply(xi, 1, sd) != 0)
      x = xi[x_ind, ]
      if(is.null(dim(x))){
        x = as.matrix(x) %>% t
        max_i = which.max(apply(cor(t(x)),1,mean, na.rm = TRUE))
        x1 = x[max_i,]
      }else if(nrow(x)==0){
        x1 = NULL
      }else{
        max_i = which.max(apply(cor(t(x)),1,mean, na.rm = TRUE))
        x1 = x[max_i,]
      }
      if(is.null(x1)){
        x_rev = xi
      }else{
        x_rev = apply(x, 1, function(y){
          y1 = y;
          y2 = y; y2[which(y==1)]=2; y2[which(y==2)]=1
          c1 = cor(x1, y1); c2 = cor(x1, y2)
          if(c2>c1) y2 else y1
        }) %>% t
      }
      return(x_rev)
    })

    ### 6. 合并
    c_combine = lapply(group_rev_list, function(tb_C){
      apply(tb_C, 2, function(x){
        t = table(x)
        ind = which.max(t)
        as.numeric(names(t)[ind])
      })
    })
    c_combine = do.call(rbind, c_combine)
    rownames(c_combine) = paste0("L_", 1:L)

    ### 7. 生成合并的c0_combine
    c0_combine <- c_combine
    for(i in 1:L){
      ind = which(T < mean(tb_list[[i]]))
      if(length(ind)!=0)
        c0_combine[i,ind] = 0
    }

  }


  ###############################################
  #####             create sdf              #####
  ###############################################
  ### 1. 准备sdf
  ord <- 1:L
  Tb <- sapply(tb_list, mean)
  lambda1 <- sapply(l1_list, mean)
  lambda2 <- sapply(l2_list, mean)
  lambda3 <- sapply(l3_list, mean)
  mt <- param_mt
  s2t <- param_s2t
  s2x <- param_s2x
  rho <- param_rho

  C <- t(c_combine)
  C0 <- t(c0_combine)

  ord <- order(Tb, decreasing = FALSE)
  Tb <- Tb[ord]
  C <- C[,ord]
  C0 <- C0[,ord]
  lambda1 <- lambda1[ord]
  lambda2 <- lambda2[ord]
  lambda3 <- lambda3[ord]

  sdf <- data.frame(T = T,
                    C0,
                    C,
                    matrix(rep(Tb, each=P),nr=P),
                    object@MURP$centersPCA$x[,1:getParams(ct, "murp_pc_number")],
                    X = X,
                    matrix(rep(lambda1, each=P),nr=P),
                    matrix(rep(lambda2, each=P),nr=P),
                    matrix(rep(lambda3, each=P),nr=P),
                    mt = rep(mt, P),
                    s2t = rep(s2t, P),
                    s2x = rep(s2x, P),
                    rho = rep(rho, P))
  colnames(sdf) = c("T",
                    paste0("C0_", 1:L),
                    paste0("C_", 1:L),
                    paste0("Tb_", 1:L),
                    paste0("PC_", 1:Q),
                    "X",
                    paste0("lambda1_", 1:L),
                    paste0("lambda2_", 1:L),
                    paste0("lambda3_", 1:L),
                    "mt","s2t","s2x","rho")

  sdf = sdf[1:nrow(sdf),1:ncol(sdf)]
  sdf$unified_direction = unified_direction

  ### 2. 准备sdf_orig
  orig_sdf <- data.frame(row.names = rownames(object@assay$data_matrix),
                         names = rownames(object@assay$data_matrix),
                         T = MapMURPLabelToAll(vecc = sdf$T, orig = murp$Recommended_K_cl$cluster) )
  for(l in 1:L){
    c0 = paste0("C0_",l)
    c = paste0("C_",l)
   orig_sdf$a = MapMURPLabelToAll(vecc = sdf[,c0], orig = murp$Recommended_K_cl$cluster)
    orig_sdf$b = MapMURPLabelToAll(vecc = sdf[,c], orig = murp$Recommended_K_cl$cluster)
    colnames(orig_sdf)[(ncol(orig_sdf)-1):ncol(orig_sdf)] = c(c0,c)
  }

  ### 3. save
  # save(sdf, file = "2_pseudotime/sdf.rda")
  # save(orig_sdf, file = "2_pseudotime/orig_sdf.rda")

  ### 4. save to object
  # save(sdf, file = "~/sdf.rda")
  object = AddMURPMetadata(object, sdf)
  object = AddMetadata(object, orig_sdf)

  command = GetCommand()
  object@Command$pse$GetPseSdf = command

  return(object)
}

#' PlotPCPseBif
#'
#' @description
#' Visualize bifurcation effects on PC and pseudotime
#'
#' @param object
#' @param pse_sdf
#' @param lm
#'
PlotPCPseBif <- function(object,
                         lm = FALSE){


  pointColorValue <- colorRampPalette(pal_d3("category10")(10))(10)[c(8,1,4)]
  names(pointColorValue) <- c(0,1,2)

  sdf = GetMURPInfo(object)
  L = getParams(object, "trajectory_number")
  Q = getParams(object, "murp_pc_number")

  ## plot
  plist <- lapply(1:L, function(l){
    lapply(1:Q, function(q){

      tb = paste0("Tb_", l)
      c = paste0("C0_", l)
      w = paste0("PC_", q)
      t = "T"
      df = sdf[, c(t, w, c, tb)]
      df[,c] = as.factor(df[,c])

      df$C_1 = NA
      df$C_2 = NA
      df$C_1[which(df[,c]==0 | df[,c]==1)] = 1
      df$C_2[which(df[,c]==0 | df[,c]==2)] = 2
      df$C_1 = as.factor(df$C_1)
      df$C_2 = as.factor(df$C_2)

      df$x_1 = as.numeric(as.character(df$C_1))
      df$x_2 = as.numeric(as.character(df$C_2))
      df[which(df$x_1==1), "x_1"] = df[which(df$x_1==1), w]
      df[which(df$x_2==2), "x_2"] = df[which(df$x_2==2), w]

      colnames(df)[5:8] = c(paste0(c,"_",1:2), paste0(w,"_",1:2))

      if(lm){
        ggplot(df) +
          geom_point(aes_string(x = t, y = w, colour = c), alpha = 0.4, size = 4) +
          geom_vline(xintercept = df[,tb][1], colour = "#990000", linetype = "dashed") +
          geom_smooth(alpha = 0.7, size = 1, method = loess, formula = y ~ x, se = FALSE, na.rm = TRUE,
                      aes_string(x = t, y = paste0(w,"_1"), colour = paste0(c,"_1")) ) +
          geom_smooth(alpha = 0.7, size = 1, method = loess, formula = y ~ x, se = FALSE, na.rm = TRUE,
                      aes_string(x = t, y = paste0(w,"_2"), colour = paste0(c,"_2")) ) +
          geom_rug(alpha = 0.2, aes_string(x = t, y = paste0(w,"_1"), colour = paste0(c,"_1")) ) +
          geom_rug(alpha = 0.2, aes_string(x = t, y = paste0(w,"_2"), colour = paste0(c,"_2")) ) +
          scale_colour_manual(values = pointColorValue) +
          labs(title = paste0("L ",l),
               # subtitle = paste0("liklyhood: ", round(pcl[l, q]),2),
               x = "PseudoT", y = w, color = "Bif") +
          # guides(colour = "none") +
          rj.ftheme
      } else{
        ggplot(df) +
          geom_point(aes_string(x = t, y = w, colour = c), alpha = 0.4, size = 4) +
          geom_vline(xintercept = df[,tb][1], colour = "#990000", linetype = "dashed") +
          # geom_smooth(alpha = 0.7, size = 1, method = loess, formula = y ~ x, se = FALSE, na.rm = TRUE,
          #             aes_string(x = t, y = paste0(w,"_1"), colour = paste0(c,"_1")) ) +
          # geom_smooth(alpha = 0.7, size = 1, method = loess, formula = y ~ x, se = FALSE, na.rm = TRUE,
          #             aes_string(x = t, y = paste0(w,"_2"), colour = paste0(c,"_2")) ) +
          geom_rug(alpha = 0.2, aes_string(x = t, y = paste0(w,"_1"), colour = paste0(c,"_1")) ) +
          geom_rug(alpha = 0.2, aes_string(x = t, y = paste0(w,"_2"), colour = paste0(c,"_2")) ) +
          scale_colour_manual(values = pointColorValue) +
          labs(title = paste0("L ",l),
               # subtitle = paste0("liklyhood: ", round(pcl[l, q]),2),
               x = "PseudoT", y = w, color = "Bif") +
          # guides(colour = "none") +
          rj.ftheme
      }


    })
  })

  ## merge
  p = c()
  for(l in 1:L){ p = c(p,plist[[l]]) }
  ggsave(paste0("2_pseudotime/2.1_julia_result/T_PC_C.pdf"),
         patchwork::wrap_plots(p, ncol = Q, nrow = L, guides = "collect"),
         width = Q*3, height = L*3.3)
}

#' PlotPCPseLabel
#'
#' @description
#' Visualize different labels on PC and pseudotime with
#'
#' @param object
#' @param pse_sdf
#' @param labels
#' @param lm
#'
PlotPCPseLabel <- function(object,
                           labels,
                           lm = FALSE){

  sdf = GetMURPInfo(object)
  L = getParams(object, "trajectory_number")
  Q = getParams(object, "murp_pc_number")

  if(!is.null(labels)){

    for(i in 1:length(labels)){
      lab = labels[i]
      cat(lab, "\n")

      ## plot
      plist <- lapply(1:L, function(l){
        lapply(1:Q, function(q){

          tb = paste0("Tb_", l)
          c = paste0("C0_", l)
          w = paste0("PC_", q)
          t = "T"
          df = sdf[, c(t, w, c, tb)]
          df[,c] = as.factor(df[,c])

          df$C_1 = NA
          df$C_2 = NA
          df$C_1[which(df[,c]==0 | df[,c]==1)] = 1
          df$C_2[which(df[,c]==0 | df[,c]==2)] = 2
          df$C_1 = as.factor(df$C_1)
          df$C_2 = as.factor(df$C_2)

          df$x_1 = as.numeric(as.character(df$C_1))
          df$x_2 = as.numeric(as.character(df$C_2))
          df[which(df$x_1==1), "x_1"] = df[which(df$x_1==1), w]
          df[which(df$x_2==2), "x_2"] = df[which(df$x_2==2), w]

          colnames(df)[5:8] = c(paste0(c,"_",1:2), paste0(w,"_",1:2))
          df = as.data.frame(df)
          df$celltype = as.factor(as.character(sdf[,lab]))

          if(lm){
            ggplot(df, aes_string(x = t, y = w, colour = "celltype")) +
              geom_point( alpha = 0.4, size = 4) +
              geom_vline(xintercept = df[,tb][1], colour = "#990000", linetype = "dashed") +
              geom_smooth(alpha = 0.7, size = 1, method = loess, formula = y ~ x, se = FALSE, na.rm = TRUE,
                          aes_string(x = t, y = paste0(w,"_1"), colour = paste0(c,"_1")) ) +
              geom_smooth(alpha = 0.7, size = 1, method = loess, formula = y ~ x, se = FALSE, na.rm = TRUE,
                          aes_string(x = t, y = paste0(w,"_2"), colour = paste0(c,"_2")) ) +
              # geom_rug(alpha = 0.2, aes_string(x = t, y = paste0(w,"_1"), colour = paste0(c,"_1")) ) +
              # geom_rug(alpha = 0.2, aes_string(x = t, y = paste0(w,"_2"), colour = paste0(c,"_2")) ) +
              # scale_colour_manual(values = pointColorValue, guide="none") +
              scale_colour_d3("category20") +
              labs(title = paste0("L ",l),
                   # subtitle = paste0("liklyhood: ", round(pcl[l, q]),2),
                   x = "PseudoT", y = w, colour = lab) +
              # guides(colour = "none") +
              rj.ftheme
          }else{
            ggplot(df, aes_string(x = t, y = w, colour = "celltype")) +
              geom_point( alpha = 0.4, size = 4) +
              geom_vline(xintercept = df[,tb][1], colour = "#990000", linetype = "dashed") +
              scale_colour_d3("category20") +
              labs(title = paste0("L ",l),
                   # subtitle = paste0("liklyhood: ", round(pcl[l, q]),2),
                   x = "PseudoT", y = w, colour = lab) +
              # guides(colour = "none") +
              rj.ftheme
          }

        })
      })

      ## merge
      p = c()
      for(l in 1:L){ p = c(p,plist[[l]]) }
      ggsave(paste0("2_pseudotime/2.1_julia_result/T_",lab,".pdf"),
             patchwork::wrap_plots(p, ncol = Q, nrow = L, guides = "collect"),
             width =Q*4, height = L*4)

    }

  }

}

#' CorViolinPlot
#'
#' @description
#' draw the correlation between truetime and pseudotime in geom_violin
#'
#' @param mamba_object
#' @param dm.match
#' @param pdfName
#'
CorViolinPlot <- function(truetime = NULL,
                          pseudotime = NULL,
                          title = "TrueTime ~ PseudoTime",
                          ftheme = NULL){

  require(ggplot2)
  rj.ftheme <- theme(panel.background = element_rect(fill='transparent', color="black"),
                     strip.text = element_text(size = 10),
                     panel.grid.minor=element_blank(),
                     panel.grid.major=element_blank(),
                     panel.border = element_rect(fill='transparent', color='black'),
                     plot.title = element_text(size = 15, hjust = 0), # title
                     plot.subtitle = element_text(size = 11.5, hjust = 0), # subtitle
                     legend.key = element_rect( fill = "white"),
                     axis.title.x = element_text(vjust = -1.5, size = 15, colour = 'black'), # face = "bold"
                     axis.title.y = element_text(vjust = 1.5, size = 15, colour = 'black'), # face = "bold"
                     axis.ticks = element_blank(),
                     axis.text.x = element_text(vjust = -0.5, size = 12, colour = 'black'),
                     axis.text.y = element_text(vjust = 0.5, size = 12, colour = 'black'),
                     legend.text = element_text(vjust = 0.4, size = 15, colour = 'black'),
                     legend.title = element_text(vjust = 0.4, size = 15, colour = 'black'),
                     legend.key.size = unit(0.9, "cm") )

  match = data.frame(truetime = truetime,
                     pseudotime = pseudotime,
                     stringsAsFactors = FALSE)

  p <- ggplot(match, aes(x = as.factor(truetime), y = pseudotime)) +
    geom_jitter(alpha = .1) +
    geom_violin(alpha = .75) +
    labs(title = title,
         subtitle = paste0('spearman cor = ',
                           round(cor(match$pseudotime,
                                     as.numeric(as.character(match$truetime)),
                                     method = 'spearman'),2), "\n",
                           'pearson cor = ',
                           round(cor(match$pseudotime,
                                     as.numeric(as.character(match$truetime)),
                                     method = 'pearson'),2)),
         x = "TrueTime", y = "PseudoTime")
  # geom_text(aes(x = 1, y= 1.0, colour=I('black'),
  #               label = paste0('spearman cor = ',
  #                              round(cor(match$pseudotime,
  #                                        as.numeric(as.character(match$truetime)),
  #                                        method = 'spearman'),2), "\n",
  #                              'pearson cor = ',
  #                              round(cor(match$pseudotime,
  #                                        as.numeric(as.character(match$truetime)),
  #                                        method = 'pearson'),2)) )) +

  if(is.null(ftheme)){
    p <- p + rj.ftheme
  }else {
    p <- p + ftheme
  }

  return(p)

}


#' #' GetPsePCASdf
#' #'
#' #' @description
#' #' get all param info from optim result
#' #'
#' #' @param object
#' #' @param param_meanvalue
#' #' @param t_pred_index
#' #' @param unified_direction
#' #' @param rev
#' #' @param labels
#' #'
#' GetPsePCASdf <- function(object,
#'                          param_meanvalue = NULL,
#'                          t_pred_index = 1,
#'                          unified_direction = FALSE,
#'                          rm_adjust_chain = FALSE,
#'                          rev = FALSE,
#'                          labels = NULL){
#'   # object = ct
#'   # t_pred_index=1
#'
#'   P <- getParams(object, "murp_number")
#'   L <- getParams(object, "trajectory_number")
#'   Q <- getParams(object, "murp_pc_number")
#'
#'   mamba = object@OptimResult$pse
#'   murp = object@MURP
#'   metadata = object@MetaData
#'   chains = mamba@chains
#'
#'   # 是否要和logpdf_ch结合过滤
#'   # logpdf_ch <- paste0("chain",1:10)[which(mamba@logpdf$cor$spearman[1,] < 0.7)]
#'   # filter_ch <- setdiff(mamba@t_pred[[1]]$filter_chain,logpdf_ch)
#'   if(mamba@t_pred[[t_pred_index]]$adjust){
#'     rev_ch <- mamba@t_pred[[t_pred_index]]$adjust_chain
#'   }
#'   if(rm_adjust_chain){
#'     filter_ch <- setdiff(mamba@t_pred[[t_pred_index]]$filter_chain,rev_ch)
#'   }else{
#'     filter_ch <- mamba@t_pred[[t_pred_index]]$filter_chain
#'   }
#'
#'   param_tb <- param_meanvalue$Tb[filter_ch,,drop=FALSE]
#'   param_T <- param_meanvalue$T[filter_ch,,drop=FALSE]
#'   param_X <- param_meanvalue$X[filter_ch,,drop=FALSE]
#'   param_lambda1 <- param_meanvalue$lambda1[filter_ch,, drop=FALSE]
#'   param_lambda2 <- param_meanvalue$lambda2[filter_ch,, drop=FALSE]
#'   param_lambda3 <- param_meanvalue$lambda3[filter_ch,, drop=FALSE]
#'   param_mt <- param_meanvalue$m_t[filter_ch,, drop=FALSE] %>% mean
#'   param_s2t <- param_meanvalue$s2_t[filter_ch,, drop=FALSE] %>% mean
#'   param_s2x <- param_meanvalue$s2_x[filter_ch,, drop=FALSE] %>% mean
#'   param_rho <- param_meanvalue$rho[filter_ch,, drop=FALSE] %>% mean
#'
#'   ### 用P确定C
#'   x <- param_meanvalue$C[filter_ch, , drop=FALSE]
#'   param_p <- param_meanvalue$p[filter_ch, , drop=FALSE]
#'   param_C <- ifelse(param_p>0.5, 1, 2)
#'   colnames(param_C) <- colnames(x)
#'
#'   ### 方向的统一性，需要同时翻转T和Tb
#'   if(unified_direction & !rm_adjust_chain){
#'     if(length(which(rev_ch %in% filter_ch))!=0 ){
#'       inter = intersect(rev_ch, filter_ch)
#'       param_tb[inter, ] = 1 - param_tb[inter, ]
#'       param_T[inter, ] = 1 - param_T[inter, ]
#'     }
#'   }
#'
#'   ### 方向翻转
#'   if(rev){
#'     param_T <- 1- param_T
#'     param_tb <- 1 - param_tb
#'   }
#'
#'   ### 取均值
#'   T <- apply(param_T, 2, mean)
#'   X <- apply(param_X, 2, mean)
#'
#'   ## 判断是不是单链
#'   if(nrow(param_tb)==1){
#'
#'     ###############################################
#'     ######           single chain             #####
#'     ###############################################
#'     tb_list = as.list( param_tb )
#'     l1_list = param_lambda1[1,]
#'     l2_list = param_lambda2[1,]
#'     l3_list = param_lambda3[1,]
#'     tmp = matrix(param_C[1,], nr=3)
#'
#'     c_combine = tmp[,]
#'     rownames(c_combine) = paste0("L_", 1:L)
#'
#'     c0_combine <- c_combine
#'     for(i in 1:L){
#'       ind = which(T < mean(tb_list[[i]]))
#'       if(length(ind)!=0)
#'         c0_combine[i,ind] = 0
#'     }
#'
#'   }else{
#'
#'     ###############################################
#'     ######             COMBINE                #####
#'     ###############################################
#'     ### 1. 对10个chain的Tb向量化并聚类，得到均值Tb
#'     param_tb_2 <- param_tb
#'     colnames(param_tb_2) <- paste0("L",1:L)
#'     tb_chain <- as.vector(param_tb_2)
#'     names(tb_chain) <- unlist(lapply(colnames(param_tb_2),
#'                                      function(x) paste0(rownames(param_tb_2),"_",x)))
#'     tb_pam <- cluster::pam(tb_chain, L, metric = "manhattan")
#'     pdf("2_pseudotime/2.2_binarytree/tb_cluster.pdf")
#'     plot(x = tb_chain, y = tb_chain, col = tb_pam$clustering)
#'     dev.off()
#'     # Tb <- tapply(tb_chain, tb_pam$clustering, mean)
#'     # names(Tb) <- paste0("Tb", 1:L)
#'
#'     ### 2. 准备一个漂亮的C
#'     C_list <- lapply(filter_ch, function(ch){
#'       tmp = matrix(param_C[ch, ], nr=L)
#'       rownames(tmp) = c(paste0(ch,"_L",1:L))
#'       tmp
#'     })
#'     names(C_list) <- c(filter_ch)
#'     c <- do.call(rbind, C_list)
#'
#'     ### 3. 根据聚类对C | Tb | lambda1 | lambda2 | lambda3分组
#'     cluster <- tb_pam$clustering
#'     split_c <- split(names(cluster), as.factor(cluster))
#'     group_list <- lapply(split_c, function(x){
#'       if(length(x)==1){
#'         tmp = matrix(c[x,],nr=1)
#'         dimnames(tmp) = list(x,colnames(c))
#'         tmp
#'       }else{
#'         c[x,]
#'       }
#'     })
#'
#'     param_tb_2 <- param_tb
#'     colnames(param_tb_2) <- paste0("L",1:L)
#'     tb <- melt(param_tb_2)
#'     tb <- data.frame(row.names = paste0(tb[,1],"_",tb[,2]), tb = tb[,3])
#'     tb_list <- split(tb[names(cluster),], as.factor(cluster))
#'
#'     param_l1_2 <- param_lambda1
#'     colnames(param_l1_2) <- paste0("L",1:L)
#'     l1 <- melt(param_l1_2)
#'     l1 <- data.frame(row.names = paste0(l1[,1],"_",l1[,2]), l1 = l1[,3])
#'     l1_list <- split(l1[names(cluster),], as.factor(cluster))
#'
#'     param_l2_2 <- param_lambda2
#'     colnames(param_l2_2) <- paste0("L",1:L)
#'     l2 <- melt(param_l2_2)
#'     l2 <- data.frame(row.names = paste0(l2[,1],"_",l2[,2]), l2 = l2[,3])
#'     l2_list <- split(l2[names(cluster),], as.factor(cluster))
#'
#'     param_l3_2 <- param_lambda3
#'     colnames(param_l3_2) <- paste0("L",1:L)
#'     l3 <- melt(param_l3_2)
#'     l3 <- data.frame(row.names = paste0(l3[,1],"_",l3[,2]), l3 = l3[,3])
#'     l3_list <- split(l3[names(cluster),], as.factor(cluster))
#'
#'     ### 5. 翻转
#'     group_rev_list <- lapply(group_list, function(x){
#'       x = as.matrix(x)
#'       x_ind = which(apply(x, 1, sd) != 0)
#'       x = x[x_ind, ]
#'       if(is.null(dim(x))){
#'         x = as.matrix(x) %>% t
#'         max_i = which.max(apply(cor(t(x)),1,mean, na.rm = TRUE))
#'         x1 = x[max_i,]
#'       }else{
#'         max_i = which.max(apply(cor(t(x)),1,mean, na.rm = TRUE))
#'         x1 = x[max_i,]
#'       }
#'       x_rev = apply(x, 1, function(y){
#'         y1 = y;
#'         y2 = y; y2[which(y==1)]=2; y2[which(y==2)]=1
#'         c1 = cor(x1, y1); c2 = cor(x1, y2)
#'         if(c2>c1) y2 else y1
#'       })
#'       t(x_rev)
#'     })
#'
#'     ### 6. 合并
#'     c_combine = lapply(group_rev_list, function(tb_C){
#'       apply(tb_C, 2, function(x){
#'         t = table(x)
#'         ind = which.max(t)
#'         as.numeric(names(t)[ind])
#'       })
#'     })
#'     c_combine = do.call(rbind, c_combine)
#'     rownames(c_combine) = paste0("L_", 1:L)
#'
#'     ### 7. 生成合并的c0_combine
#'     c0_combine <- c_combine
#'     for(i in 1:L){
#'       ind = which(T < mean(tb_list[[i]]))
#'       if(length(ind)!=0)
#'         c0_combine[i,ind] = 0
#'     }
#'
#'   }
#'
#'
#'   ###############################################
#'   #####             create sdf              #####
#'   ###############################################
#'   ### 1. 准备sdf
#'   ord <- 1:L
#'   Tb <- sapply(tb_list, mean)
#'   lambda1 <- sapply(l1_list, mean)
#'   lambda2 <- sapply(l2_list, mean)
#'   lambda3 <- sapply(l3_list, mean)
#'   mt <- param_mt
#'   s2t <- param_s2t
#'   s2x <- param_s2x
#'   rho <- param_rho
#'
#'   C <- t(c_combine)
#'   C0 <- t(c0_combine)
#'
#'   ord <- order(Tb, decreasing = FALSE)
#'   Tb <- Tb[ord]
#'   C <- C[,ord]
#'   C0 <- C0[,ord]
#'   lambda1 <- lambda1[ord]
#'   lambda2 <- lambda2[ord]
#'   lambda3 <- lambda3[ord]
#'
#'   sdf <- data.frame(T = T,
#'                     C0,
#'                     C,
#'                     matrix(rep(Tb, each=P),nr=P),
#'                     object@MURP$Recommended_K_cl$centers[,1:getParams(ct, "murp_pc_number")],
#'                     X = X,
#'                     matrix(rep(lambda1, each=P),nr=P),
#'                     matrix(rep(lambda2, each=P),nr=P),
#'                     matrix(rep(lambda3, each=P),nr=P),
#'                     mt = rep(mt, P),
#'                     s2t = rep(s2t, P),
#'                     s2x = rep(s2x, P),
#'                     rho = rep(rho, P))
#'   colnames(sdf) = c("T",
#'                     paste0("C0_", 1:L),
#'                     paste0("C_", 1:L),
#'                     paste0("Tb_", 1:L),
#'                     paste0("PC_", 1:Q),
#'                     "X",
#'                     paste0("lambda1_", 1:L),
#'                     paste0("lambda2_", 1:L),
#'                     paste0("lambda3_", 1:L),
#'                     "mt","s2t","s2x","rho")
#'   if(!is.null(labels)){
#'     tmp = data.frame(murp$murp_cellinfo[,labels])
#'     colnames(tmp) = labels
#'     sdf = cbind(sdf, tmp)
#'   }
#'   sdf = sdf[1:nrow(sdf),1:ncol(sdf)]
#'   sdf$unified_direction = unified_direction
#'
#'   ### 2. 准备sdf_orig
#'   orig_sdf <- data.frame(row.names = rownames(object@assay$data_matrix),
#'                          names = rownames(object@assay$data_matrix),
#'                          T = GetPseudoTt(CellNum = nrow(object@assay$data_matrix),
#'                                          pred_t = sdf$T,
#'                                          Cluster_Label = murp$Recommended_K_cl$cluster))
#'   for(l in 1:L){
#'     c0 = paste0("C0_",l)
#'     c = paste0("C_",l)
#'     orig_sdf$a = GetPseudoTt(CellNum = nrow(object@assay$data_matrix), pred_t = sdf[,c0], Cluster_Label = murp$Recommended_K_cl$cluster)
#'     orig_sdf$b = GetPseudoTt(CellNum = nrow(object@assay$data_matrix), pred_t = sdf[,c], Cluster_Label = murp$Recommended_K_cl$cluster)
#'     colnames(orig_sdf)[(ncol(orig_sdf)-1):ncol(orig_sdf)] = c(c0,c)
#'   }
#'
#'   if(!is.null(labels)){
#'     tmp = data.frame( object@MetaData[,labels])
#'     colnames(tmp) = labels
#'     orig_sdf = cbind(orig_sdf, tmp)
#'   }
#'
#'   ### 3. save
#'   save(sdf, file = "2_pseudotime/sdf.rda")
#'   save(orig_sdf, file = "2_pseudotime/orig_sdf.rda")
#'
#'   return(sdf)
#' }

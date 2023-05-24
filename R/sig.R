
#' ComputeInteractionEffectGene
#'
#' @param object
#'
#' @export
#'
ComputeInteractionEffectGene <- function(object,
                                         adjust_method = "bonferroni",
                                         save_mod = FALSE){


  P = getParams(object, "murp_number")
  L = getParams(object, "trajectory_number")
  Q = getParams(object, "murp_pc_number")

  ## 1. get meanvalue and track_sdf
  sdf = GetMURPInfo(object)

  ## 2. calculating
  cat(":: lm \n")
  gene <- colnames(object@MURP$Recommended_K_cl$centers)
  base_c <- "C0_"
  if(base_c=="C0_"){ enu = 3; lev = 0:2 }
  if(base_c=="C_"){ enu = 2; lev = 1:2 }

  mod_all <- lapply(1:L, function(l){

    cat("\n")
    cat("--", l, "\n")
    c = paste0(base_c, l)

    ## 1. 计算LM
    cat("////// Compute \n")
    modlist <- lapply(gene, function(g){
      tb = sdf[1,paste0("Tb_",l)]
      df = data.frame(T = sdf$T,
                      C = factor(sdf[,c],levels = lev),
                      g = object@MURP$Recommended_K_cl$centers[,g])
      tryCatch(
        {
          lm(g ~ C + T + C * T, data = df)
        },
        error = function(e){
          NULL
        }
      )
    })
    names(modlist) <- gene

    ## 2. 提取pvalue
    cat("////// Adjust Pvalue \n")
    pvalue_list <- lapply(gene, function(g){
      # cat(g, "\n")
      mod = modlist[[g]]
      mod_df = summary(mod)
      tryCatch(
        {
          pv_name = c("(Intercept)", paste0("C",1:(enu-1)),
                      "T", paste0("C",1:(enu-1), ":T"))
          pv = rep(NA, length(pv_name))
          names(pv) = pv_name
          tmp = mod_df$coefficients[,4]
          pv[names(tmp)] = tmp
          pv
        },
        error = function(e){
          NULL
        }
      )
    })

    ## 3. 整理 & p.adjust
    sig_df <- tryCatch( {
      df = do.call(rbind, pvalue_list) %>% data.frame
      colnames(df) = c("pvalue_intercept", paste0("pvalue_c",1:(enu-1)),
                       "pvalue_t", paste0("pvalue_c",1:(enu-1), ".t"))
      rownames(df) = gene

      df_adjust = df
      for(j in 1:ncol(df_adjust)){
        df_adjust[,j] = p.adjust(df_adjust[,j], method = adjust_method)
      }
      colnames(df_adjust) = c(paste0(adjust_method, "_intercept"),
                              paste0(adjust_method, "_c",1:(enu-1)),
                              paste0(adjust_method, "_t"),
                              paste0(adjust_method, "_c",1:(enu-1), ".t"))
      cbind(df_adjust, df)
    }, error = function(e){
      NULL
    })

    ## 4. 结果
    if(save_mod){
      list(mod = modlist, sig_df = sig_df)
    }else{
      list(mod = NULL, sig_df = sig_df)
    }
  })
  object@BeautGene$lm <- list(mod_all = mod_all)
  return(object)
}

#' ComputeInteractionEffectGene
#'
#' @param object
#'
#' @export
#'
ComputeAOVEffectGene <- function(object,
                                 adjust_method = "bonferroni",
                                 save_mod = FALSE){


  P = getParams(object, "murp_number")
  L = getParams(object, "trajectory_number")
  Q = getParams(object, "murp_pc_number")

  ## 1. get meanvalue and track_sdf
  sdf = GetMURPInfo(object)

  ## 2. calculating
  cat(":: anova c \n")
  gene <- colnames(object@MURP$Recommended_K_cl$centers)
  mod_all <- lapply(1:L, function(l){

    cat("\n")
    cat("--", l, "\n")

    base_c <- "C0_"
    if(base_c=="C0_"){ enu = 3; lev = 0:2 }
    if(base_c=="C_"){ enu = 2; lev = 1:2 }
    c = paste0(base_c, l)

    ## 1. 计算aov
    cat("////// Compute \n")
    modlist <- lapply(gene, function(g){
      # cat(g, "\n")
      tb = sdf[1,paste0("Tb_",l)]
      df = data.frame(T = factor(ifelse(sdf$T>=tb,2,1), levels = 1:2),
                      C = factor(sdf[,c],levels = lev),
                      g = object@MURP$Recommended_K_cl$centers[,g])
      tryCatch(
        {
          aov(g ~ C, data = df)
        },
        error = function(e){
          NULL
        }
      )
    })
    names(modlist) <- gene

    ## 2. 提取pvalue
    cat("////// Adjust Pvalue \n")
    pvalue_list <- lapply(gene, function(g){
      # cat(g, "\n")
      mod = modlist[[g]]
      mod_df = summary(mod)
      tryCatch(
        {

          pv_name = c("C", "Residuals")
          pv = rep(NA, length(pv_name))
          names(pv) = pv_name
          tmp = mod_df[[1]][,5]
          names(tmp) = gsub(" ","",rownames(mod_df[[1]]))
          pv[names(tmp)] = tmp
          pv[c("C")]

        },
        error = function(e){
          NULL
        }
      )
    })

    ## 3. 整理 & fdr
    sig_df <- tryCatch( {
      df = do.call(rbind, pvalue_list) %>% data.frame
      colnames(df) = c("pvalue_c")
      rownames(df) = gene

      df_adjust = df
      for(j in 1:ncol(df_adjust)){
        df_adjust[,j] = p.adjust(df_adjust[,j], method = adjust_method)
      }
      colnames(df_adjust) = paste0(adjust_method,"_c")
      cbind(df_adjust, df)
    }, error = function(e){
      NULL
    })

    ## 4. 结果
    if(save_mod){
      list(mod = modlist, sig_df = sig_df)
    }else{
      list(mod = NULL, sig_df = sig_df)
    }

  })
  object@BeautGene$aov_c <- list(mod_all = mod_all)
  return(object)
}

#' ComputeInteractionEffectGene
#'
#' @param object
#'
#' @export
#'
ComputeAOV2EffectGene <- function(object,
                                  adjust_method = "bonferroni",
                                  save_mod = FALSE){


  P = getParams(object, "murp_number")
  L = getParams(object, "trajectory_number")
  Q = getParams(object, "murp_pc_number")

  ## 1. get meanvalue and track_sdf
  sdf = GetMURPInfo(object)

  ## 2. calculating
  cat(":: anova c \n")
  gene <- colnames(object@MURP$Recommended_K_cl$centers)
  mod_all <- lapply(1:L, function(l){

    cat("\n")
    cat("--", l, "\n")

    base_c <- "C0_"
    if(base_c=="C0_"){ enu = 3; lev = 0:2 }
    if(base_c=="C_"){ enu = 2; lev = 1:2 }
    c = paste0(base_c, l)

    ## 1. 计算aov
    cat("////// Compute \n")
    modlist <- lapply(gene, function(g){
      # cat(g, "\n")
      tb = sdf[1,paste0("Tb_",l)]
      df = data.frame(T = factor(ifelse(sdf$T>=tb,2,1), levels = 1:2),
                      C = factor(sdf[,c],levels = lev),
                      g = object@MURP$Recommended_K_cl$centers[,g])
      df = df[which(df$C!=0),] ## 只保留后C1，C2
      df$C = factor(df$C, levels = c(1,2))
      tryCatch(
        {
          aov(g ~ C, data = df)
        },
        error = function(e){
          NULL
        }
      )
    })
    names(modlist) <- gene

    ## 2. 提取pvalue
    cat("////// Adjust Pvalue \n")
    pvalue_list <- lapply(gene, function(g){
      # cat(g, "\n")
      mod = modlist[[g]]
      mod_df = summary(mod)
      tryCatch(
        {

          pv_name = c("C", "Residuals")
          pv = rep(NA, length(pv_name))
          names(pv) = pv_name
          tmp = mod_df[[1]][,5]
          names(tmp) = gsub(" ","",rownames(mod_df[[1]]))
          pv[names(tmp)] = tmp
          pv[c("C")]

        },
        error = function(e){
          NULL
        }
      )
    })

    ## 3. 整理 & fdr
    sig_df <- tryCatch( {
      df = do.call(rbind, pvalue_list) %>% data.frame
      colnames(df) = c("pvalue_c")
      rownames(df) = gene

      df_adjust = df
      for(j in 1:ncol(df_adjust)){
        df_adjust[,j] = p.adjust(df_adjust[,j], method = adjust_method)
      }
      colnames(df_adjust) = paste0(adjust_method,"_c")
      cbind(df_adjust, df)
    }, error = function(e){
      NULL
    })

    ## 4. 结果
    if(save_mod){
      list(mod = modlist, sig_df = sig_df)
    }else{
      list(mod = NULL, sig_df = sig_df)
    }
  })
  object@BeautGene$aov_c2 <- list(mod_all = mod_all)
  return(object)
}

#' ComputeMeanDiff
#'
#' @param object
#'
#' @export
#'
ComputeMeanDiff <- function(object,
                            method = "wilcox.test",
                            adjust_method = "bonferroni"){

  P = getParams(object, "murp_number")
  L = getParams(object, "trajectory_number")
  Q = getParams(object, "murp_pc_number")
  # expr_dat = object@assay[[assay]]
  expr_dat = object@MURP$Recommended_K_cl$centers
  genes = colnames(expr_dat)

  ## 1. get meanvalue and track_sdf
  sdf = GetMURPInfo(object)

  ## 2. calculating
  # cat(":: compute fold change \n")
  result <- lapply(1:L, function(l){

    cat("\n")
    cat("---- Trajectory ", l, "\n\n")

    ## set C
    base_c <- "C0_"
    if(base_c=="C0_"){ enu = 3; lev = 0:2 }
    if(base_c=="C_"){ enu = 2; lev = 1:2 }
    c = paste0(base_c, l)

    ## mean difference
    cat("////// mean difference\n")
    tmp_list <- lapply(genes, function(g){
      tb = sdf[1,paste0("Tb_",l)]
      df = data.frame(C = factor(sdf[,c],levels = lev),
                      g = expr_dat[,g])
      df = df[which(df$C!=0),] ## 只保留后C1，C2
      df$C = factor(df$C, levels = c(1,2))

      ind1 = which(df$C==1)
      ind2 = which(df$C==2)
      diff = mean(df$g[ind1]) - mean(df$g[ind2])
      c(gene = g, diff = diff,
        status = ifelse(diff>0, "Up in C1", ifelse(diff<0, "Up in C2", "Non")) )
    })

    diff_df = do.call(rbind, tmp_list) %>% data.frame
    rownames(diff_df) = diff_df$gene
    diff_df[,2] = as.numeric(diff_df[,2])

    ## add diff test
    cat("//////", method, "\n")
    tmp_list <- lapply(genes, function(g){
      # cat(g, "\n")
      tb = sdf[1,paste0("Tb_",l)]
      df = data.frame(C = factor(sdf[,c],levels = lev),
                      g = expr_dat[,g])
      df = df[which(df$C!=0),] ## 只保留后C1，C2
      df$C = factor(df$C, levels = c(1,2))
      ind1 = which(df$C==1)
      ind2 = which(df$C==2)
      pval = tryCatch(
        {
          if(method=="wilcox.test"){
            mod = wilcox.test(df$g[ind1], df$g[ind2])
            p = mod$p.value
          }
          if(method=="t.test"){
            mod = t.test(df$g[ind1], df$g[ind2])
            p = mod$p.value
          }
          p
        },
        error = function(e){
          NA
        }
      )
      c(gene = g, pval = pval)
    })
    test_df = do.call(rbind, tmp_list) %>% data.frame
    rownames(test_df) = test_df$gene
    test_df[,2] = as.numeric(test_df[,2])

    ## p.adjust
    cat("////// adjust pvalue\n")
    test_df$p.adj = p.adjust(test_df[,2], method = adjust_method)

    ## combine
    df = cbind(diff_df, test_df[,-1])
    colnames(df) = c("gene", "diff", "status", paste0(method,".pvalue"), paste0(method,".",adjust_method))
    return(list(sig_df = df))

    # if(all(expr_dat>=0)){
    #   fclist <- lapply(genes, function(g){
    #     # cat(g, "\n")
    #     tb = sdf[1,paste0("Tb_",l)]
    #     df = data.frame(C = factor(sdf[,c],levels = lev),
    #                     g = expr_dat[,g])
    #     df = df[which(df$C!=0),] ## 只保留后C1，C2
    #     df$C = factor(df$C, levels = c(1,2))
    #
    #     ind1 = which(df$C==1)
    #     ind2 = which(df$C==2)
    #     fc = mean(df$g[ind1]+1)/mean(df$g[ind2]+1)
    #     lfc = log2(fc)
    #     c(gene = g, foldchange = fc, logfoldchange=lfc)
    #   })
    #   fcdf = do.call(rbind, fclist) %>% data.frame
    #   rownames(fcdf) = fcdf$gene
    #   fcdf[,2] = as.numeric(fcdf[,2])
    #   fcdf[,3] = as.numeric(fcdf[,3])
    # }else{
    #   fclist <- lapply(genes, function(g){
    #     # cat(g, "\n")
    #     tb = sdf[1,paste0("Tb_",l)]
    #     df = data.frame(C = factor(sdf[,c],levels = lev),
    #                     g = expr_dat[,g])
    #     df = df[which(df$C!=0),] ## 只保留后C1，C2
    #     df$C = factor(df$C, levels = c(1,2))
    #
    #     ind1 = which(df$C==1)
    #     ind2 = which(df$C==2)
    #     lfc = mean(df$g[ind1]) - mean(df$g[ind2])
    #     c(gene = g, logfoldchange=lfc)
    #   })
    #   fcdf = do.call(rbind, fclist) %>% data.frame
    #   rownames(fcdf) = fcdf$gene
    #   fcdf[,2] = as.numeric(fcdf[,2])
    # }
  })
  object@BeautGene$diff$mod_all <- result
  return(object)
}

#' ComputeMeanDiff
#'
#' @param object
#'
#' @export
#'
ComputeCorGeneTraj <- function(object,
                               method = "spearman",
                               adjust_method = "bonferroni"){

  P = getParams(object, "murp_number")
  L = getParams(object, "trajectory_number")
  Q = getParams(object, "murp_pc_number")
  # expr_dat = object@assay[[assay]]
  expr_dat = ct@MURP$Recommended_K_cl$centers
  genes = colnames(expr_dat)

  ## 1. get meanvalue and track_sdf
  sdf = GetMURPInfo(object)

  ## 2. calculating
  result <- lapply(1:L, function(l){

    cat("\n")
    cat("---- Trajectory ", l, "\n\n")
    w = sdf[, paste0("W_",l)]
    ## set C
    # base_c <- "C0_"
    # if(base_c=="C0_"){ enu = 3; lev = 0:2 }
    # if(base_c=="C_"){ enu = 2; lev = 1:2 }
    # c = paste0(base_c, l)

    ## cor test
    cat("////// cor test\n")
    tmp_list <- lapply(genes, function(g){
      r = cor.test(expr_dat[,g], w, method = method)
      c(gene = g, cor = r$estimate, pvalue = r$p.value)
    })
    cor_df = do.call(rbind, tmp_list) %>% data.frame
    dimnames(cor_df) = list(cor_df$gene,c("gene","cor","pvalue"))
    cor_df[,2] = as.numeric(cor_df[,2])

    ## p.adjust
    cat("////// adjust pvalue\n")
    cor_df$p.adj = p.adjust(cor_df$pvalue, method = adjust_method)

    ## combine
    colnames(cor_df) = c("gene", paste0(method,"_cor"), paste0(method,"_pvalue"), paste0(method,"_",adjust_method))
    return(list(sig_df = cor_df))

  })
  object@BeautGene$cor_gene_traj$mod_all <- result
  return(object)
}

#' ComputeMeanDiff
#'
#' @param object
#'
#' @export
#'
ComputeCorGeneT <- function(object,
                               method = "spearman",
                               adjust_method = "bonferroni"){

  sdf = GetMURPInfo(object)

  # gene weight
  gene_weight = GetGeneWeight(object)
  # gene_sd = apply(object@MURP$Recommended_K_cl$centers, 2, sd)
  # gene_weight = sapply(1:3, function(i)gene_weight[,i]*gene_sd)
  colnames(gene_weight) = paste0("L",1:3)
  # gene_num = (nrow(gene_weight) * weight_ratio) %>% round

  expr_mat = object@MURP$Recommended_K_cl$centers
  t = sdf$T
  lapply(colnames(expr_mat), function(g){
    r = cor.test(t, expr_mat[,g], method = "spearman")
    c(gene = g, cor = r$estimate, pvalue = r$p.value)
  }) -> tmp_list

  cor_df = do.call(rbind, tmp_list) %>% data.frame
  dimnames(cor_df) = list(cor_df$gene,c("gene","cor","pvalue"))
  cor_df[,2] = as.numeric(cor_df[,2])

  ## p.adjust
  cat("////// adjust pvalue\n")
  cor_df$p.adj = p.adjust(cor_df$pvalue, method = adjust_method)
  colnames(cor_df) = c("gene", paste0(method,"_cor"), paste0(method,"_pvalue"), paste0(method,"_",adjust_method))

  object@BeautGene$cor_gene_pseudot <- cor_df
  return(object)

}

#' ComputeSigGene
#' @description
#' compute significance genes with lm or aov
#'
#' @param cs celltrek-gene object
#' @param ct_method dge method
#' @param sig_method statistics label
#' @param sig_th filter threshold
#'
#' sig = TRUE
#' ct_method = "lm"
#' sig_method = "pvalue"
#' sig_th = 0.05
#'
#' @export
#'
ComputeSigGene <- function(cs,
                           ct_method,
                           sig_method,
                           sig_th){

  obj = cs[[ct_method]]
  L = length(obj$mod_all)

  ## 按照列取满足阈值的gene
  sig_gene <- lapply(1:L, function(l){

    # if(ct_method == "diff"){
    #   sig_df = obj$mod_all[[l]]
    # }else{
    #   sig_df = obj$mod_all[[l]]$sig_df
    # }
    sig_df = obj$mod_all[[l]]$sig_df
    tryCatch(
      {
        tmp = sig_df[,grep(sig_method,colnames(sig_df)),drop=FALSE]
        tmp_list = lapply(1:ncol(tmp), function(i){
          x = tmp[,i]
          ind = order(x)
          x2 = x[ind]
          tmp2 = tmp[ind,,drop=FALSE]
          rownames(tmp2)[which(x2 < sig_th)]
        })
        names(tmp_list) = stringr::str_split_fixed(colnames(tmp),"_",2)[,2]
        if(ct_method=="lm"){
          tmp_list$ct = tmp_list[grep("\\.t",names(tmp_list))] %>% unlist %>% unique
          tmp_list$c = tmp_list[grep("^c[0-9]$",names(tmp_list))] %>% unlist %>% unique
        }
        tmp_list
      },
      error = function(e){
        NULL
      })
  })


  ## 统计个数
  sig_gene_num <- do.call(rbind,
                          lapply(1:L, function(l){
                            x = sig_gene[[l]]
                            x_len = sapply(x, length)
                            if(length(x_len)==0) NA else x_len
                          }))
  rownames(sig_gene_num) = paste0("L",1:L)

  ## 保存结果
  result = list(sig_gene = sig_gene,
                sig_gene_num = sig_gene_num)

  obj = Save2CS(obj,
                "sig_result",
                paste0(sig_method,"_", sig_th),
                result)

  # if("sig_result" %nin% names(obj)){
  #   obj$sig_result = list()
  #   obj$sig_result[[1]] = result
  #   names(obj$sig_result) = paste0(sig_method,"_", sig_th)
  # }else{
  #   tmp_n = paste0(sig_method,"_", sig_th)
  #   tmp_all = c(names(obj$sig_result), tmp_n)
  #   obj$sig_result[[length(obj$sig_result)+1]] = result
  #   names(obj$sig_result) = make.unique(tmp_all)
  # }
  cs[[ct_method]] = obj

  return(cs)
}

#' AdjustPvalue
#' @description
#' adjust pvalue using different method
#'
#' @param cs celltrek-gene object
#' @param ct_method dge method
#' @param sig_method statistics label
#' @param adjust_method filter threshold
#'
#' ct_method = "lm"
#' sig_method = "pvalue"
#'
#' @export
AdjustPvalue <- function(cs,
                         ct_method,
                         adjust_method){

  obj = cs[[ct_method]]
  L = length(obj$mod_all)

  # adjust
  for(i in 1:L){
    sig_df = obj$mod_all[[i]]$sig_df

    # sig_df = obj$mod_all[[i]]$sig_df
    ind = grep("pvalue",colnames(sig_df))
    for(j in ind){
      sig_df[,j] = p.adjust(sig_df[,j], method = adjust_method)
    }
    sig_df_adjust = sig_df[,ind,drop=FALSE]
    colnames(sig_df_adjust) = gsub("pvalue",adjust_method,colnames(sig_df_adjust))
    obj$mod_all[[i]]$sig_df = cbind(sig_df, sig_df_adjust)
  }

  cs[[ct_method]] = obj
  return(cs)

}

#' ComputeSigData
#' @description
#' get new sdata and cl according to differential gene
#'
#' @param ct celltrek object
#' @param cs celltrek-gene object
#' @param ct_method dge method
#' @param sig_label dge factor
#' @param sig_method statistics label
#' @param sig_th filter threshold
#'
#' sig = TRUE
#' ct_method = "lm"
#' sig_method = "pvalue"
#' sig_th = 0.05
#'
#' @export
#'
ComputeSigData <- function(ct,
                           cs,
                           sig_label,
                           ct_method,
                           sig_method,
                           sig_th){


  obj = cs[[ct_method]]
  if( ("sig_result" %nin% names(obj)) | (paste0(sig_method,"_",sig_th) %nin% names(obj$sig_result) ) ){
    cs = ComputeSigGene(cs, ct_method, sig_method, sig_th)
  }

  ## get dge gene
  L = length(obj$mod_all)
  cat("L", l, "-", ct_method, "-", sig_label, "-", sig_method, "-", sig_th, "\n")
  ct_gene <- lapply(1:L, function(l){
    cat(l, "\n")
    GetSigGene(cs,
               traj = l,
               ct_method = ct_method,
               sig_method= sig_method,
               sig_label = sig_label,
               sig_th = sig_th)

  }) %>% unlist %>% unique

  ## get new cl
  # Y1 = ct@MURP$Recommended_K_cl$centers
  # Yx = ct@MURP$centersPCA$x
  # U = ct@MURP$centersPCA$rotation
  # all.equal(Y1 %*% U[,1:3], Yx[,1:3])

  # new_centers <- ct@MURP$Recommended_K_cl$centers[,ct_lm_gene]
  # new_U <- ct@MURP$centersPCA$rotation[ct_lm_gene,]
  # cl <- new_centers %*% new_U

  sd_cutoff <- 1
  center <- FALSE
  scale <- FALSE
  new_centers <- ct@MURP$Recommended_K_cl$centers[,ct_gene]
  new_pca <- prcomp(new_centers, center = center, scale. = scale)
  npc <- GetNPC(pca_result = new_pca, sd_cutoff = sd_cutoff)
  cl <- new_pca$x[,1:npc]
  sdata <- ct@assay$sdata[,ct_gene]

  result <- list(gene = ct_gene,
                 cl = cl,
                 sdata = sdata)

  obj = Save2CS(obj,
                "new_data",
                paste0(sig_label, "_", sig_method,"_", sig_th),
                result)
  cs[[ct_method]] = obj

  return(cs)

  # if(!("new_data" %in% names(obj))){
  #   obj$new_data = list()
  #   obj$new_data[[1]] = result
  #   names(obj$new_data) = paste0(sig_label, "_", sig_method,"_", sig_th)
  # }else{
  #
  #   tmp_name = paste0(sig_label, "_", sig_method,"_", sig_th)
  #   if(tmp_name %in% names(obj$new_data)){
  #     obj$new_data[[tmp_name]] = result
  #   }else{
  #     obj$new_data[[length(obj$new_data)+1]] = result
  #     names(obj$new_data)[length(obj$new_data)] = tmp_name
  #   }
  # }

}


#' GetSigGene
#'
#' @param cs celltrek-gene object
#' @param ct_method dge method
#' @param sig_label dge factor
#' @param sig_method statistics label
#' @param sig_th filter threshold
#'
#' sig = TRUE
#' ct_method = "lm"
#' sig_method = "pvalue"
#' sig_th = 0.05
#'
#' @export
#'
GetSigGene <- function(cs,
                       traj,
                       ct_method,
                       sig_method,
                       sig_label,
                       sig_th){
  l = traj
  obj = cs[[ct_method]]
  tag = paste0(sig_method,"_",sig_th)

  ## 如果没有现场算
  if( !(tag %in% names(obj$sig_result)) ){
    cs = ComputeSigGene(cs, ct_method, sig_method, sig_th)
  }

  fgene = obj$sig_result[[tag]]$sig_gene
  if(ct_method=="diff"){
    sig_gene = fgene[[l]][[1]]
  }else{
    sig_gene = fgene[[l]][[sig_label]]
  }

  return(sig_gene)
}


#' GetSigDF
#'
#' @param cs celltrek-gene object
#' @param ct_method dge method
#' @param traj trajectory number
#'
#' sig = TRUE
#' ct_method = "lm"
#' sig_method = "pvalue"
#' sig_th = 0.05
#'
#' @export
#'
GetSigDF <- function(cs,
                     traj,
                     ct_method){
  if(ct_method=="cor_gene_pseudot"){
    sig_df = cs[[ct_method]]
  }else{
    sig_df = cs[[ct_method]]$mod_all[[traj]]$sig_df
  }
  return(sig_df)
}

#' Save2CS
#' save result to cs
#'
#' @param obj
#' @param base
#' @param name
#' @param result
#'
#' @export
#'
Save2CS <- function(obj,
                    base,
                    name,
                    result){

  if(base %nin% names(obj)){
    obj[[base]] = list()
    obj[[base]][[1]] = result
    names(obj[[base]]) = name
  }else{

    if(name %in% names(obj[[base]])){
      obj[[base]][[name]] = result
    }else{
      obj[[base]][[length(obj[[base]])+1]] = result
      names(obj[[base]])[length(obj[[base]])] = name
    }
  }
  return(obj)
}

#' SaveSigSDF
#'
#' @param ct celltrek object
#' @param cs celltrek-gene object
#' @param ct_method dge method
#' @param sig_label dge factor
#' @param sig_method statistics label
#' @param sig_th filter threshold
#'
#' sig = TRUE
#' ct_method = "lm"
#' sig_method = "pvalue"
#' sig_th = 0.05
#'
#' @export
#'
SaveSigSDF <- function(ct,
                       cs,
                       sdf,
                       sig_label,
                       ct_method,
                       sig_method,
                       sig_th){

  obj = cs[[ct_method]]
  if( ("sig_result" %nin% names(obj)) | (paste0(sig_method,"_",sig_th) %nin% names(obj$sig_result)) ){
    cs = ComputeSigGene(cs, ct_method, sig_method, sig_th)
  }
  if( ("new_data" %nin% names(obj)) | (paste0(sig_method,"_",sig_th) %nin% names(obj$new_data)) ){
    cs = ComputeSigData(ct, cs, sig_label, ct_method, sig_method, sig_th)
  }

  new_data_index = paste0(sig_label, "_",sig_method,"_",sig_th)
  sdata = cs[[ct_method]]$new_data[[new_data_index]]$sdata
  cl = cs[[ct_method]]$new_data[[new_data_index]]$cl
  colnames(cl) = paste0("PC_",1:ncol(cl))
  sdf[, grep("PC", colnames(sdf))] = NULL
  sdf = cbind(sdf, cl)

  save(sdf, file = "sdf.rda")
  save(cl, file = paste0('cl_', nrow(cl) ,'cell.rda'))
  save(sdata, file = "sdata.rda")

  return(cs)
}

#' WriteGene2TXT
#'
#' @param ct celltrek object
#' @param cs celltrek-gene object
#' @param ct_method dge method
#' @param sig_label dge factor
#' @param sig_method statistics label
#' @param sig_th filter threshold
#'
#' sig = TRUE
#' ct_method = "lm"
#' sig_method = "pvalue"
#' sig_th = 0.05
#'
#' @export
#'
WriteGene2TXT <- function(cs,
                          ct_method,
                          sig_method,
                          sig_label,
                          sig_th,
                          ENTREZID = FALSE){

  base_n = paste0(ct_method,"_",sig_method,"_", sig_label,"_", sig_th)
  obj = cs[[ct_method]]
  fgene = obj$sig_result[[paste0(sig_method,"_",sig_th)]]$sig_gene
  L = length(fgene)

  for(l in 1:L){
    tryCatch({
      if(ct_method=="diff"){
        gene = fgene[[l]][[1]]
      }else{
        gene = fgene[[l]][[sig_label]]
      }
      write.table(gene,
                  file = paste0("4_differential_genes/",base_n,"_SYMBOL_L",l,".txt"),
                  row.names = FALSE, col.names = FALSE, quote = FALSE)
      if(ENTREZID){
        marker = bitr(gene, fromType = "SYMBOL", toType = c("ENSEMBL", "ENTREZID"), OrgDb = org.Hs.eg.db)
        write.table(marker$ENTREZID,
                    file = paste0("4_differential_genes/", base_n,"_ENTREZID_L",l,".txt"),
                    row.names = FALSE, col.names = FALSE, quote = FALSE)
      }

    },
    error = function(e){
      NULL
    })
  }
}

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' ComputeWeightGene
#'
#' @description
#' 原来的0_2_track_main_Q.jl以及老版的kalman.jl得到的结果计算权重基因的方法
#'
#' @param object celltrak object
#' @param iter_range iterations range
#'
# ComputeWeightGene <- function(object = NULL,
#                               iter_range = NULL){
#
#   P = getParams(object, "murp_number")
#   L = getParams(object, "trajectory_number")
#   Q = getParams(object, "murp_pc_number")
#
#   ## 1. get meanvalue and track_sdf
#   load(paste0("3_tracking/sdf_",iter_range[1], "_",iter_range[2],".rda"))
#   load(paste0("3_tracking/param_meanvalue_",iter_range[1], "_", iter_range[2],".rda"))
#
#   ## 2. get R
#   if("R" %in% names(param_meanvalue)){
#     W <- matrix(apply(param_meanvalue$W, 2, mean), nc = L) %>% as.data.frame
#     colnames(W) <- paste0("W_", 1:L)
#     R <- matrix(apply(param_meanvalue$R, 2, mean), nc = Q) %>% as.data.frame
#     rownames(R) <- paste0("R_", 1:Q)
#     R <- as.matrix(R[,1:L])
#
#   }else{
#     W = sdf[,grep("W_[1-9]$",colnames(sdf))] %>% t
#     Yx = sdf[,grep("PC_",colnames(sdf))] %>% t
#     R = do.call(cbind,
#                 lapply(1:L, function(l){
#                   ginv(Yx %*% t(Yx)) %*% Yx %*% W[l,]
#                 }) )
#     R <- as.matrix(R[,1:L])
#   }
#
#   ## 3. gene weight
#   U <- object@MURP$centersPCA$rotation[,1:Q]
#   weight <- U %*% R
#
#   ## 4. gene weight order
#   gene_weight_order <- apply(abs(weight),2,function(x){
#     rownames(weight)[order(x,decreasing = TRUE)] })
#
#   ## 5. save to object
#   gene_weight <- object@OptimResult$track@gene_weight
#   len = length(gene_weight)
#   gene_weight[[len+1]] <- weight
#   names(gene_weight)[len+1] <- paste0("iter_", iter_range[1], "_", iter_range[2])
#   object@OptimResult$track@gene_weight = gene_weight
#
#   ## 6. save
#   save(weight, file = paste0("3_tracking/3.3_gene_weight/gene_weight_",iter_range[1],"_",iter_range[2],".rda"))
#   write.csv(gene_weight_order, file = paste0("3_tracking/3.3_gene_weight/gene_weight_order_", iter_range[1],"_",iter_range[2],".csv") )
#
#   return(object)
# }

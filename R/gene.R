
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#                                    gpr
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' gaussian_process_regression
#'
#' @description
#' Use Gaussian process regression to calculate the smooth expression of any gene/trajectory score along a trajectory.
#' @param s covariance matrix
#' @param P murp number
#' @param gene_expr expression vector of any gene in the 'murp' dataset.
#' @param n.samples The number of points included in the regression results
#' @export
gaussian_process_regression <- function(s, P,
                                        gene_expr,
                                        n.samples = 100){

  P_ = nrow(s)-P

  # compute k
  k.xx <- s[1:P, 1:P]
  k.xxs <- s[1:P, (P+1):ncol(s)]
  k.xsx <- s[(P+1):ncol(s), 1:P]
  k.xsxs <- s[(P+1):ncol(s), (P+1):ncol(s)]

  # update mu and sigma
  cov.f.star <- k.xsxs - k.xsx %*% solve(k.xx) %*% k.xxs
  f.star.bar <- k.xsx %*% solve(k.xx) %*% gene_expr

  # mvrnorm
  values <- matrix(rep(0, P_*n.samples), ncol=n.samples)
  for (n in 1:n.samples) {
    cat("sample: ", n, "\n")
    if(!is.positive.definite(cov.f.star)){
      cov.f.star = make.positive.definite(cov.f.star, 1e-6)
    }
    values[,n] <- mvrnorm(1, f.star.bar, cov.f.star)
  }
  colnames(values) = paste0("sample_",1:n.samples)
  return(values)

}

#' GetMURPGene
#'
#' @description
#' Map expression vectors of all genes onto a subset of MURP
#' @param object MGPfact object
#' @export
GetMURPGene <- function(object){
  if("data_matrix_all_gene" %nin% names(object@assay)){
    object@MURP$data_matrix = object@MURP$Recommended_K_cl$center
  }else{
    dat = object@assay$data_matrix_all_gene
    c = object@MURP$Recommended_K_cl$cluster
    c = factor(c, levels = 1:length(unique(c)) )
    tapply(1:nrow(dat), as.factor(c),
           function(x, y) apply(dat[x,,drop=F],2,mean) ) -> tmp
    tmp = do.call(rbind, tmp)
    object@MURP$data_matrix = tmp
  }
  return(object)
}

#' GetTrajectoryCurve
#'
#' @description
#' Smooth any trajectory score using Gaussian process regression
#' @param object MGPfact object
#' @param n.samples The number of points included in the regression results
#' @export
GetTrajectoryCurve <- function(object,
                               n.samples = 30){

  P = getParams(object, "murp_number")
  L = getParams(object, "trajectory_number")
  Q = getParams(object, "murp_pc_number")

  ## load cov matrix and ssdf
  S_list_ = object@GPR$reg_cov$cov_l
  sdf = GetMURPInfo(object)

  ## curve
  tmp_list <- lapply(1:L, function(l){
    x = gaussian_process_regression(s = S_list_[[l]], P = P, n.samples = n.samples,
                                    gene_expr = sdf[,paste0("W_",l)])
    apply(x, 1, mean)
  })
  trajectory_curve = do.call(rbind, tmp_list)
  rownames(trajectory_curve) = paste0("W_",1:L)

  ## save to ssdf
  object@GPR$trajectory_curve = trajectory_curve
  object@GPR$curve_sdf = cbind(object@GPR$curve_sdf,t(trajectory_curve))

  save(trajectory_curve, file = "2_pseudotime/2.4_cov_matrix/trajectory_curve.rda")

  return(object)
}

#' gaussian_process_regressionDF
#'
#' @description
#' Use Gaussian process regression to calculate the smooth
#' expression of any gene/trajectory score along a trajectory.
#' @param object MGPfact object
#' @param n.samples The number of points included in the regression results
#' @param genes a character contain gene names
#' @export
GetGeneCurve <- function(object,
                         genes = NULL,
                         n.samples = 30){

  ## parameters
  P = getParams(object, "murp_number")
  L = getParams(object, "trajectory_number")
  Q = getParams(object, "murp_pc_number")

  ## curve
  S_list_ = object@GPR$reg_cov$cov_l
  curve_gene_expr_list <- lapply(1:length(S_list_), function(l){
    tmp_list <- lapply(genes, function(g){
      x = gaussian_process_regression(s = S_list_[[l]],
                                      P = P,
                                      n.samples = n.samples,
                                      gene_expr = object@MURP$data_matrix[,g])
      apply(x, 1, mean)
    })
    g_expr = do.call(rbind, tmp_list )
    rownames(g_expr) = genes
    g_expr
  })

  ## save plot gene df
  ssdf = object@GPR$curve_sdf
  if(is.null(ssdf)){
    object=GetCurveSdf(object)
    ssdf = object@GPR$curve_sdf
  }
  plot_gene_df <- data.frame(ssdf,
                             curve_gene_expr_list[[1]] %>% t,
                             curve_gene_expr_list[[2]] %>% t,
                             curve_gene_expr_list[[3]] %>% t)
  colnames(plot_gene_df) <- c(colnames(ssdf),
                              paste0("L1_", rownames(curve_gene_expr_list[[1]])),
                              paste0("L2_", rownames(curve_gene_expr_list[[2]])),
                              paste0("L3_", rownames(curve_gene_expr_list[[3]])))
  object@GPR$curve_sdf = plot_gene_df
  save(plot_gene_df, file = paste0("2_pseudotime/2.4_cov_matrix/plot_gene_df.rda"))

  ## save gpr
  object@GPR$gene_expr = curve_gene_expr_list
  save(curve_gene_expr_list, file = paste0("2_pseudotime/2.4_cov_matrix/curve_gene_expr_list.rda"))

  return(object)
}

#' GetCurveSdf
#'
#' @description
#' Integrate the smooth results from Gaussian process regression
#' @param object MGPfact object
#' @export
GetCurveSdf <- function(object){

  L = getParams(object, "trajectory_number")
  s_ = object@GPR$reg_cov

  C0 <- t(s_$C)
  for(i in 1:L){
    ind = which(s_$T_ < s_$Tb[i])
    if(length(ind)!=0)
      C0[ind,i] = 0
  }

  ssdf <- data.frame(T = s_$T,
                     C0,
                     t(s_$C),
                     matrix(rep(s_$Tb, each=s_$P),nrow=s_$P),
                     X = s_$X,
                     matrix(rep(s_$lambda, each=s_$P),nrow=s_$P),
                     mt = rep(s_$m_t, s_$P),
                     s2t = rep(s_$s2_t, s_$P),
                     s2x = rep(s_$s2_x, s_$P),
                     rho = rep(s_$rho, s_$P))

  colnames(ssdf) = c("T",
                     paste0("C0_", 1:L),
                     paste0("C_", 1:L),
                     paste0("Tb_", 1:L),
                     "X",
                     paste0("lambda1_", 1:L),
                     paste0("lambda2_", 1:L),
                     paste0("lambda3_", 1:L),
                     "mt","s2t","s2x","rho")

  save(ssdf, file = "2_pseudotime/2.4_cov_matrix/ssdf.rda")
  object@GPR$curve_sdf = ssdf
  return(object)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#                              expression plot
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' GeneExprMURP
#'
#' @description
#' Obtain the expression vector of a gene in MURP
#' @param gene any gene
#' @param murp murp result
#' @param expr expression matrix
#' @param methods mean/sum expression
#' @export
#'
GeneExprMURP <- function(gene = NULL,
                         expr = NULL,
                         murp = NULL,
                         methods = NULL){

  if(!(gene %in% rownames(expr))){
    stop(gene, " not found in expr")
  }
  if(!(methods %in% c("mean","sum"))){
    stop("wrong methods")
  }
  if(is.null(methods)){
    methods = "sum"
  }

  if(length(gene)==1){
    value = tapply(1:ncol(expr),
                   as.factor(murp$Recommended_K_cl$cluster),
                   function(x, y){
                     if(methods=="sum"){
                       sum(expr[gene,x])
                     }else{
                       mean(expr[gene,x])
                     }
                   })
  }else{
    tmp = tapply(1:ncol(expr),
                 as.factor(murp$Recommended_K_cl$cluster),
                 function(x, y){
                   if(methods=="sum"){
                     # sum(count[gene,x])
                     apply(expr[gene,x], 1, sum)
                   }else{
                     apply(expr[gene,x], 1, mean)
                   }
                 })
    value = do.call(rbind, tmp)
  }

  return(value)
}

#' GeneBar
#'
#' @description
#' Bar plot of gene expression in different branches.
#' @param object MGPfact object
#' @param gene any gene
#' @param murp murp result
#' @param group group name
#' @param exprs expression matrix
#' @param expr_method mean/sum expression
#' @param barcolor color of bar
#' @param errorbar logical value, whether to add errorbar
#' @export
GeneBar <- function(object,
                    gene = NULL,
                    group = NULL,
                    exprs = NULL,
                    murp = NULL,
                    expr_method = "mean",
                    barcolor = colorRampPalette(pal_rickandmorty()(12))(12),
                    errorbar = FALSE){

  expr = GeneExprMURP(gene = gene,
                      expr = exprs,
                      murp = object@MURP,
                      methods = expr_method)

  cdf <- data.frame(group, expr)
  gdf1 <- cdf %>% select("expr", "group") %>% group_by(group ) %>%
    summarise( mean(expr), sum(expr)) %>%  `colnames<-`(c("group", "mean", "sum")) %>% data.frame
  gdf2 <- cdf %>% select("expr", "group") %>% group_by(group ) %>%
    summarySE(measurevar = "expr", groupvars = "group") %>%  data.frame

  gdf <- merge(gdf1, gdf2, by = "group")
  gdf[,"group"] <- factor(gdf[,"group"], levels = c(0,1,2))

  p = ggplot(gdf, aes_string(x = "group", y = expr_method, fill = "group")) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = barcolor) +
    guides(fill = "none") +
    rj.ftheme

  if(expr_method == "mean"){
    p = p +
      labs(x = "", y = "Mean Expression",
           fill = "Group", title = gene)
  }

  if(expr_method == "sum"){
    p = p +
      labs(x = "", y = "Total Expression",
           fill = "Group", title = gene)
  }

  if(errorbar){
    p = p +
      geom_errorbar(aes(ymin = get(expr_method) - .data$ci, ymax = get(expr_method) + .data$ci),
                    width = .1, position = position_dodge(.5))
  }

  return(p)
}

#' GeneVln
#'
#' @description
#' violin plot of gene expression in different branches.
#' @param expr expression matrix including several genes
#' @param group group name
#' @param title plot title
#' @param ylab Y-axis title
#' @param vlncolor color of violin
#' @param errorbar logical value, whether to add errorbar
#' @param max_y  Y-axis maximum value.
#' @param pvalue logical value, whether to add pvalue label
#' @export
GeneVln <- function(expr = NULL,
                    group = NULL,
                    title = "x",
                    ylab = NULL,
                    vlncolor = colorRampPalette(pal_rickandmorty()(12))(12),
                    errorbar = FALSE,
                    max_y = 3,
                    pvalue = TRUE){

  if(is.null(ylab)){ ylab = "Expr"}
  gdf = data.frame(group = group, expr = expr)
  gdf$group = factor(gdf$group, levels = c(0,1,2))

  p = ggplot(gdf,aes(x = factor(group,levels = 0:2), y = expr, group = group)) +
    geom_jitter(aes(colour = group), alpha = 0.05) +
    geom_violin(aes(fill = group), alpha = 0.7) +
    geom_boxplot(width = 0.08) +
    scale_colour_manual(values = vlncolor, drop = F) +
    scale_fill_manual(values = vlncolor, drop = F) +
    scale_x_discrete(breaks = c(0:2), labels = paste0("Phase ",0:2)) +
    scale_y_continuous(limits = c(min(expr), max(expr) + max_y)) +
    labs(x = NULL, y = ylab, fill = "Group", title = title) +
    guides(colour = "none", fill = "none") +
    rj.ftheme

  if(pvalue){
    if(length(unique(group))==2){
      p = p +
        stat_compare_means(comparisons = list(c("1","2")),
                           method = "wilcox.test", label = "p.signif",hide.ns = T)
    }else{
      p = p +
        stat_compare_means(comparisons = list(c("1","2"),c("0","1"),c("0","2")),
                           method = "wilcox.test", label = "p.signif",hide.ns = T)
    }
  }

  if(errorbar){
    p = p +
      stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), geom="pointrange")
  }

  return(p)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#                             beaut area plot
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' BeautAreaMulti
#'
#' @description
#' Display the expression density variation of any one gene across the time domain,
#' combined with the density variation of cell types
#' @param df Data frame containing all genes and trajectory information,
#' transformed into a long matrix with separate columns for pseudotime and branch label.
#' @param mdf dataframe for any one gene on any one trajectory,
#' including T, branch label, celltype, gene expression
#' @param g any gene name
#' @param cut_tb logical value, Whether to remove the section before the bifurcation
#' @param cols color of gene expression density in different branches
#' @param y_title_vjust Distance of y-axis labels from the panel.
#' @param max_str Maximum width for cell type strings
#' @export
BeautAreaDrq <- function(df, mdf, g, cut_tb = FALSE,
                         cols = c("#374e55ff","#df8f44ff"),
                         y_title_vjust = -55,
                         max_str = 50){

  # g = gsub("-","\\.",g)
  df = df[which(df$variable==g),]
  mdf = mdf[which(mdf$variable==g),]
  tb = df$Tb[1]

  #### cut_tb
  if(cut_tb){
    df = df[which(df$T>tb),]
    mdf = mdf[which(mdf$T>tb),]
  }

  #### scale_y
  # y_ma = max(abs(df$value)) %>% round(1)
  y_ma = max(abs(df$value)) %>% ceiling
  y_mi = -y_ma
  mseq = seq(y_mi, y_ma, length = 5) %>% round(2)
  mseq = mseq[2:4]
  ####  scale_x
  x_mi = min(df$T) %>% round(2)
  x_ma = 1

  ####  pmain_1
  celltype_max_length = max(nchar(as.character(unique(mdf$celltype))))
  pmain_1 <- ggplot()+

    # 1. plot
    geom_ribbon(data = df, size = 2, alpha = 0.7, ymin = y_mi,
                aes(x = .data$T, ymax = .data$value, fill = .data$C)) +

    # 2. scale
    scale_fill_manual(values = cols, limits=levels(df$C), drop=TRUE) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0),
                       limits = c(y_mi, y_ma), breaks = mseq, labels = mseq, position = "left") +

    # 3. signig
    annotate(geom="text", x = quantile(df$T,probs=seq(0,1,0.05))[8], y = y_ma-2*y_ma/10, size = 5,
             label = paste0("Anova, ",
                            # ifelse(df$sig[1]<2.22e-16, "p-value < 2.22e-16", paste0("p-value = ", round(df$sig[1],7))) )) +
                            ifelse(df$sig[1]<2.22e-16, "p < 2.22e-16", paste0("p = ", format(df$sig[1],digits=2))) )) +

    # 4. theme
    guides(color = "none",
           fill = guide_legend(override.aes=list(alpha = 0.5, size = 3), nrow = 2)) +
    labs(x = "", y = "Scaled Expression", fill = paste0(g)) +
    theme(plot.margin = unit(c(0,0.5,0,0.5), "cm"),
          panel.background = element_rect(fill='transparent', color="black"),
          panel.border = element_rect(fill='transparent', color='black'),
          panel.grid.minor=element_blank(),
          panel.grid.major=element_blank(),
          axis.title.y = element_text(size = 13, color='black', angle = 90, hjust = 0.5,
                                      vjust = y_title_vjust),
          axis.ticks = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 13, color='black', angle = 0, hjust = 1),
          legend.background = element_rect(fill='transparent', color= NULL),
          legend.position=c(0.85, 0.8),
          legend.key = element_rect(fill='transparent', color = NULL),
          legend.text = element_text(vjust = 0.4, size = 13, colour = 'black'),
          legend.title = element_text(vjust = 0.4, size = 13, colour = 'black') )
  if(!cut_tb){
    pmain_1 = pmain_1 +
      geom_vline(xintercept = tb, color = "darkgreen", linetype = "dashed")
  }
  ####  xdens
  xdens_1 <- ggplot() +
    geom_density_ridges(data = mdf,
                        aes(x = .data$T, y  = .data$celltype, fill = .data$C, height = .data$..density..),
                        # stat = "density",
                        panel_scaling = FALSE, adjust = 2,
                        alpha = 0.7,  size = 0.01) +
    labs(y="", fill = "") +
    guides(fill = "none") +
    scale_x_continuous(limits = c(x_mi,1), expand = c(0,0)) +
    scale_y_discrete(expand = c(0.05, 0), position = "left",
                     labels=ggplot2:::parse_safe) +
                     # labels=function(x) stringr::str_wrap(x, width=max_str)) +
    scale_fill_manual(values = cols,
                      limits=levels(mdf$C), drop=TRUE) +
    theme(plot.margin = unit(c(0.5,0,0,0.2), "cm"),
          plot.title = element_text(size = 13, hjust = 0),
          strip.background = element_blank(),
          strip.text = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(color='black', size = 0.3),
          axis.title = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y.left = element_text(size = 13, color='black', angle = 0, hjust = 1))
  # p <- xdens_1 / pmain_1 + plot_layout(heights = c(.3, .5))
  p <- xdens_1 / pmain_1 + plot_layout(heights = c(.3, .4))
  return(p)
}
# -------------
#' #' GeneCurveHeatmap
#' #'
#' #' @description
#' #'
#' #' @param set_list
#' #' @export
#' GeneCurveHeatmap <- function(x,
#'                              expr_scale = FALSE,
#'                              cluster_rows = TRUE,
#'                              hclust_method = "ward.D2",
#'                              num_clusters = 6,
#'                              hmcols = NULL,
#'                              add_annotation_row = NULL,
#'                              add_annotation_col = NULL,
#'                              row_dist = "euclidean",
#'                              show_rownames = FALSE,
#'                              scale_max = 3,
#'                              scale_min = -3,
#'                              return_heatmap = FALSE){
#'
#'   # scale_max = 3
#'   # scale_min = -3
#'   # cluster_rows = TRUE
#'   # hclust_method = "ward.D2"
#'   # num_clusters = 3
#'   # num_clusters = min(num_clusters, nrow(heatmap_matrix))
#'   # add_annotation_row = NULL
#'   # add_annotation_col = NULL
#'   # show_rownames = F
#'
#'   # m = log10(x + 1)
#'   if(expr_scale){
#'     m = Matrix::t(scale(Matrix::t(x), center = TRUE))
#'   }else{m = x}
#'   m[m > scale_max] = scale_max
#'   m[m < scale_min] = scale_min
#'   heatmap_matrix <- m
#'
#'   if(is.null(row_dist)){
#'     row_dist <- as.dist((1 - cor(Matrix::t(heatmap_matrix)))/2)
#'     row_dist[is.na(row_dist)] <- 1
#'   }
#'
#'   if (is.null(hmcols)) {
#'     # bks <- seq(-3.1, 3.1, by = 0.1)
#'     bks <- seq(scale_min, scale_max, by = 0.1)
#'     hmcols <- colorRamps::blue2green2red(length(bks) - 1)
#'   } else {
#'     bks <- seq(scale_min, scale_max, length.out = length(hmcols))
#'   }
#'
#'   ph <- pheatmap(heatmap_matrix,
#'                  useRaster = T,
#'                  cluster_cols = FALSE,
#'                  cluster_rows = cluster_rows,
#'                  show_rownames = F,
#'                  show_colnames = F,
#'                  clustering_distance_rows = row_dist,
#'                  clustering_method = hclust_method,
#'                  cutree_rows = num_clusters,
#'                  silent = TRUE,
#'                  filename = NA,
#'                  breaks = bks,
#'                  border = FALSE,
#'                  border_color = NA,
#'                  color = hmcols)
#'   if (cluster_rows & (!is.na(num_clusters))) {
#'     annotation_row <- data.frame(Cluster = factor(cutree(ph$tree_row,
#'                                                          num_clusters)))
#'   } else {
#'     annotation_row <- NULL
#'   }
#'
#'   if (!is.null(add_annotation_row)) {
#'     old_colnames_length <- ncol(annotation_row)
#'     annotation_row <- cbind(annotation_row, add_annotation_row[row.names(annotation_row),
#'     ])
#'     colnames(annotation_row)[(old_colnames_length + 1):ncol(annotation_row)] <- colnames(add_annotation_row)
#'   }
#'   if (!is.null(add_annotation_col)) {
#'     annotation_col <- add_annotation_col
#'   } else {
#'     annotation_col <- NA
#'   }
#'
#'   feature_label <- row.names(heatmap_matrix)
#'   row_ann_labels <- row.names(heatmap_matrix)
#'   if (!is.null(annotation_row))
#'     row_ann_labels <- row.names(annotation_row)
#'   row.names(heatmap_matrix) <- feature_label
#'
#'   colnames(heatmap_matrix) <- c(1:ncol(heatmap_matrix))
#'   ph_res <- pheatmap(heatmap_matrix[, ],
#'                      useRaster = T,
#'                      cluster_cols = FALSE,
#'                      cluster_rows = cluster_rows,
#'                      show_rownames = show_rownames,
#'                      show_colnames = F,
#'                      clustering_distance_rows = row_dist,
#'                      clustering_method = hclust_method,
#'                      cutree_rows = num_clusters,
#'                      annotation_row = annotation_row,
#'                      annotation_col = annotation_col,
#'                      treeheight_row = 20,
#'                      breaks = bks,
#'                      fontsize = 6,
#'                      color = hmcols,
#'                      border = FALSE,
#'                      border_color = NA,
#'                      silent = TRUE,
#'                      filename = NA)
#'
#'   grid::grid.rect(gp = grid::gpar("fill", col = NA))
#'   grid::grid.draw(ph_res$gtable)
#'   if (return_heatmap) {
#'     return(ph_res)
#'   }
#'
#' }
#' #' GeneCurveHeatmap2
#' #'
#' #' @description
#' #' @export
#' GeneCurveHeatmap2 <- function(x,
#'                               expr_scale = FALSE,
#'                               row_title = NULL,
#'                               cluster_rows = TRUE,
#'                               cluster_columns = FALSE,
#'                               column_title = paste0("L",l),
#'                               hclust_method = "ward.D2",
#'                               hmcols = NULL,
#'                               top_annotation = NULL,
#'                               right_annotation = NULL,
#'                               left_annotation = NULL,
#'                               row_dist = "euclidean",
#'                               col_split = NULL,
#'                               show_rownames = FALSE,
#'                               show_colnames = FALSE,
#'                               fontsize = 6,
#'                               scale_max = 3,
#'                               scale_min = -3,
#'                               pic_height = unit(6, "cm"),
#'                               pic_width = unit(6, "cm")){
#'   require(ComplexHeatmap)
#'   require(RColorBrewer)
#'
#'   # scale data
#'   if(expr_scale){
#'     m = Matrix::t(scale(Matrix::t(x), center = TRUE))
#'   }else{m = x}
#'   m[m > scale_max] = scale_max
#'   m[m < scale_min] = scale_min
#'   heatmap_matrix <- m
#'
#'   # row dist
#'   if(is.null(row_dist)){
#'     row_dist <- as.dist((1 - cor(Matrix::t(heatmap_matrix)))/2)
#'     row_dist[is.na(row_dist)] <- 1
#'   }
#'
#'   # color
#'   if (is.null(hmcols)) {
#'     bks <- seq(scale_min, scale_max, by = 0.1)
#'     hmcols <- colorRamps::blue2green2red(length(bks) - 1)
#'   } else {
#'     bks <- seq(scale_min, scale_max, length.out = length(hmcols))
#'   }
#'
#'   # plot
#'   ph <- Heatmap(heatmap_matrix,
#'                 use_raster = TRUE,
#'                 heatmap_legend_param = list(
#'                   title = "Expression",
#'                   title_position = "topcenter",
#'                   legend_direction = "horizontal", # horizontal vertical
#'                   at = c(scale_min, 0, scale_max),
#'                   labels = c("low", "0", "high")
#'                 ),
#'                 top_annotation = top_annotation,
#'                 right_annotation = right_annotation,
#'                 left_annotation = left_annotation,
#'                 name = paste0("all"),
#'                 col = hmcols,
#'                 column_title = column_title,
#'                 height = pic_height,
#'                 width = pic_width,
#'                 border = FALSE,
#'                 show_row_dend = FALSE,
#'                 cluster_rows = cluster_rows,
#'                 cluster_columns = cluster_columns,
#'                 clustering_method_rows = hclust_method,
#'                 clustering_method_columns = hclust_method,
#'                 clustering_distance_rows = row_dist,
#'                 row_names_gp = gpar(fontsize = fontsize),
#'                 column_split = col_split,
#'                 column_gap = unit(1.5, "mm"),
#'                 show_row_names = show_rownames,
#'                 show_column_names = show_colnames,
#'                 row_title = row_title)
#'
#'   return(ph)
#'
#' }
#' BeautAreaMulti
#'
#' @description
#'
#' @param df
#' @param mdf
#' @param g
#' @param cut_tb
#' @export
# BeautAreaMulti <- function(df, mdf, cut_tb = TRUE,
#                            right_space = 2){
#
#   # df = df[which(df$variable==g),]
#   # mdf = mdf[which(mdf$variable==g),]
#
#   if(cut_tb){
#     df = df[which(df$T>tb),]
#     mdf = mdf[which(mdf$T>tb),]
#   }
#
#   #### split df
#   mdf1 <- mdf[which(mdf$C!=2),]; mdf1$celltype = as.factor(mdf1$celltype); mdf1$C = 1
#   mdf2 <- mdf[which(mdf$C!=1),]; mdf2$celltype = as.factor(mdf2$celltype); mdf2$C = 2
#   celltype_col <- pal_jco("default")(10)[1:enum(mdf$celltype)]
#
#   #### scale_y
#   # y_ma = max(abs(df$value)) %>% round(1)
#   y_ma = max(abs(df$value)) %>% ceiling
#   y_mi = -y_ma
#   mseq = seq(y_mi, y_ma, length = 5)
#
#   #### scale_x
#   x_mi = min(df$T) %>% round(2)
#   x_ma = 1
#
#   #### plot
#   pmain_1 <- ggplot() +
#     geom_ribbon(data = df, size = .9, alpha = 0.7, ymin = y_mi,
#                 aes(x = T, ymax = value)) +
#     scale_y_continuous(limits = c(y_mi, y_ma),position = "right") +
#     scale_x_continuous(limits = c(x_mi, 1),expand = c(0.02, 0)) +
#     facet_grid(variable ~ C, labeller = labeller(.cols = c("1" = "C = 1", "2" = "C = 2"))) +
#     theme(plot.title = element_text(size = 15, hjust = 0),
#           plot.margin = unit(c(0,right_space,0,0.2), "cm"),
#           strip.text.x = element_text(size = 12, color='black'),
#           strip.text.y.right = element_text(size = unit(10,"pt"), color='black', angle = 0, hjust = 0),
#           # margin = margin(r = main_l,l = 5)),
#           strip.background.y = element_blank(),
#           # strip.background.x = element_rect(fill='transparent', color='black'),
#           panel.background = element_rect(fill='transparent', color="black"),
#           panel.border = element_rect(fill='transparent', color='black'),
#           panel.grid.minor=element_blank(),
#           panel.grid.major=element_blank(),
#           axis.title = element_blank(),
#           axis.ticks = element_blank(),
#           axis.text = element_blank(),
#           legend.position = "none",
#           legend.key = element_rect(fill='transparent', color="white"),
#           legend.text = element_text(vjust = 0.4, size = 13, colour = 'black'),
#           legend.title = element_text(vjust = 0.4, size = 13, colour = 'black') )
#
#   xdens_1 <- ggplot() +
#     geom_density_ridges2(data = rbind(mdf1,mdf2),
#                          aes(x = T, y  = celltype, fill = celltype, height = ..density..),
#                          stat = "density", panel_scaling = FALSE, adjust = 1.5,
#                          # rel_min_height = 0.05,
#                          alpha = 0.7,  size = 0.2) +
#     facet_grid(~C) +
#     labs(y="", fill = "") +
#     guides(fill = "none") +
#     scale_x_continuous(limits = c(x_mi,1), expand = c(0.03, 0)) +
#     scale_y_discrete(position = "right") +
#     # scale_y_discrete(position = "right", labels=function(x) stringr::str_wrap(x, width=10)) +
#     scale_fill_jco() +
#     theme(plot.margin = unit(c(0.5,0,0,0.2), "cm"),
#           plot.title = element_text(size = 15, hjust = 0),
#           strip.background = element_blank(),
#           strip.text = element_blank(),
#           panel.background = element_blank(),
#           panel.border = element_blank(),
#           panel.grid.minor=element_blank(),
#           panel.grid.major=element_blank(),
#           axis.title = element_blank(),
#           axis.ticks = element_blank(),
#           axis.text.x = element_blank(),
#           axis.text.y.right = element_text(size = 10, color='black', angle = 0, hjust = 0,
#                                            margin = margin(r = 15, l = unit(5,"pt"))))
#
#   empty <- ggplot() +
#     theme(panel.background=element_blank(),
#           axis.title.x=element_blank(),
#           axis.title.y=element_blank(),
#           axis.text.x=element_blank(),
#           axis.text.y=element_blank(),
#           axis.ticks=element_blank())
#
#   main_h = enum(df$variable)*0.5
#   p <- grid.arrange(xdens_1, empty, pmain_1,
#                     ncol=2, nrow=2,
#                     widths=c(4,0), heights=c(2, main_h)) %>% as.ggplot
#
#   return(p)
# }
#' BeautArea
#'
#' @description
#'
#' @param df
#' @param mdf
#' @param g
#' @param cut_tb
#'
# BeautArea <- function(df, mdf, g, cut_tb = TRUE,
#                       cols = c("#374e55ff","#df8f44ff"),
#                       top_air = 3,
#                       y_title_hjust = 1.5,
#                       celltype_legend_ncol = 1,
#                       y_title_vjust = -5){
#
#   require(ggplotify)
#
#   df = df[which(df$variable==g),]
#   mdf = mdf[which(mdf$variable==g),]
#
#   #### cut_tb
#   if(cut_tb){
#     df = df[which(df$T>tb),]
#     mdf = mdf[which(mdf$T>tb),]
#   }
#
#   #### scale_y
#   # y_ma = max(abs(df$value)) %>% round(1)
#   y_ma = max(abs(df$value)) %>% ceiling
#   y_mi = -y_ma
#   mseq = seq(y_mi, y_ma, length = 5)
#
#   ####  scale_x
#   x_mi = min(df$T) %>% round(2)
#   x_ma = 1
#
#   ####  split
#   df1 <- df[which(df$C==1),]
#   df2 <- df[which(df$C==2),]
#   mdf1 <- mdf[which(mdf$C!=2),]; mdf1$celltype = as.factor(mdf1$celltype)
#   mdf2 <- mdf[which(mdf$C!=1),]; mdf2$celltype = as.factor(mdf2$celltype)
#   # celltype_col <- pal_jco("default")(10)[1:enum(mdf$celltype)]
#   # celltype_col <- pal_d3("category20")(20)[1:enum(mdf$celltype)]
#   # celltype_col <- rev(pal_d3("category20")(20))[1:enum(mdf$celltype)]
#   celltype_col <- pal_d3("category20")(20)[1:length(levels(mdf$celltype))]
#   names(celltype_col) <- levels(mdf$celltype)
#
#   ####  pmain_1
#   pmain_1 <- ggplot()+
#
#     # 1. plot
#     # geom_area(data = df1, size = .9, alpha = 0.7, position = "identity",
#     #           aes(x = T, y = value, fill = C)) +
#     geom_ribbon(data = df1, size = .9, alpha = 0.5, ymin = y_mi,
#                 aes(x = T, ymax = value, fill = C)) +
#     geom_point(data = mdf1, aes(x = T, y = value, color = celltype), alpha = 0) +
#
#     # 2. scale
#     scale_fill_manual(values = cols, #c("#374e55ff","#df8f44ff"),  c("#ff940a", "#8cc68c")
#                       limits=levels(df1$C), drop=TRUE) +
#     scale_y_continuous(limits = c(y_mi,y_ma), breaks = mseq[2:5], labels = mseq[2:5],expand = c(0, 0)) +
#     scale_x_continuous(limits = c(x_mi,x_ma)) +
#
#     # 3. theme
#     guides(color = "none",
#            fill = guide_legend(override.aes=list(alpha = 0.5, size = 3),
#                                # nrow = length(levels(mdf$celltype)),
#                                nrow = 2)) +
#     labs(x = "", y = "", fill = paste0(g)) +
#     theme(plot.margin = unit(c(top_air, 0.5,0,0.5), "cm"),
#           plot.title = element_text(size = 15, hjust = 0),
#           panel.background = element_rect(fill='transparent', color="black"),
#           panel.border = element_rect(fill='transparent', color='black'),
#           panel.grid.minor=element_blank(),
#           panel.grid.major=element_blank(),
#           axis.title = element_blank(),
#           axis.ticks = element_blank(),
#           axis.text.x = element_blank(),
#           axis.text.y = element_text(vjust = 0.5, size = 16, colour = 'black'),
#           legend.position = "top",
#           legend.key = element_rect(fill='transparent', color="white"),
#           legend.text = element_text(vjust = 0.4, size = 16, colour = 'black'),
#           legend.title = element_text(vjust = 0.4, size = 16, colour = 'black') )
#
#   #### pmain_2
#   pmain_2 <- ggplot()+
#
#     # 1. plo2
#     # geom_area(data = df2, size = .9, alpha = 0.7,
#     #           aes(x = T, y = value, fill = C)) +
#     geom_ribbon(data = df2, size = .9, alpha = 0.5, ymin = -y_mi,
#                 aes(x = T, ymax = value, fill = C)) +
#     geom_point(data = mdf2, aes(x = T, y = value, color = celltype), alpha = 0) +
#
#     # 2. scale
#     scale_color_manual(values = celltype_col, breaks = names(celltype_col), drop = FALSE) +
#     scale_fill_manual(values = cols,  limits=levels(df2$C), drop=TRUE) +
#     scale_y_reverse( limits = c(y_ma, y_mi), expand = c(0, 0),breaks = rev(mseq), labels = rev(mseq)) +
#     scale_x_continuous(limits = c(x_mi,x_ma)) +
#
#     # 3. theme
#     guides(fill = "none",
#            color = guide_legend("",override.aes=list(alpha = 1, size = 4,
#                                                      color = celltype_col),ncol = celltype_legend_ncol )) +
#     labs(x = "", y = "", color = "") +
#     theme(plot.margin = unit(c(-0.1,0.5,0,0.5), "cm"),
#           panel.background = element_rect(fill='transparent', color="black"),
#           panel.border = element_rect(fill='transparent', color='black'),
#           panel.grid.minor=element_blank(),
#           panel.grid.major=element_blank(),
#           axis.title = element_blank(),
#           axis.ticks = element_blank(),
#           axis.text.x = element_blank(),
#           axis.text.y = element_text(vjust = 0.5, size = 16, colour = 'black'),
#           legend.position = "bottom",
#           legend.key = element_rect(fill='transparent', color="white"),
#           legend.text = element_text(vjust = 0.4, size = 16, colour = 'black'),
#           legend.title = element_text(vjust = 0.4, size = 16, colour = 'black'))
#
#   ## xdens
#   xdens_1 <- axis_canvas(pmain_1, axis = "x") +
#     geom_density(data = mdf1, aes(x = T, y = ..scaled.., fill = celltype),
#                  alpha = 0.7, size = 0.1, adjust = 1.5) +
#     scale_x_continuous(limits = c(x_mi,x_ma)) +
#     scale_fill_manual(values = celltype_col,drop=F)
#   # scale_fill_d3("category20c")
#   xdens_2 <- axis_canvas(pmain_2, axis = "x") +
#     geom_density(data = mdf2, aes(x = T, y = ..scaled.., fill = celltype),
#                  alpha = 0.7, size = 0.1, adjust = 1.5) +
#     scale_y_reverse(expand = c(0, 0) ) +
#     scale_x_continuous(limits = c(x_mi,x_ma)) +
#     scale_fill_manual(values = celltype_col,drop=F)
#   # scale_fill_d3("category20")
#
#   p1 <- insert_xaxis_grob(pmain_1, xdens_1,
#                           grid::unit(.3, "null"), position = "top") %>% as.ggplot
#   p2 <- insert_xaxis_grob(pmain_2, xdens_2,
#                           grid::unit(.3, "null"), position = "bottom") %>% as.ggplot
#
#   p <- p1/p2 + labs(y = "Scaled Expression")
#   p$theme$axis.title = NULL
#   p$theme$axis.title.x = element_blank()
#   p$theme$axis.title.y.left = element_text(vjust = y_title_vjust,hjust = y_title_hjust, angle = 90, size = 16, colour = 'black')
#
#   return(p)
# }


#'PCANestedGridPlot
#'
#' @description
#' Partial zoom
#'
#' @param pca_result pca result
#' @param sd_cutoff choose the right number of pc according sd
#' @param max_pc maximum number of pc
#'
#' @export
#'
PCANestedGridPlot <- function(pca_result = NULL,
                              sd_cutoff = 1,
                              max_pc = NULL){

  #pca <- MGPfact_object@DimReduc$PCA$pca_result
  pca <- pca_result
  npc <- length(which(pca$sdev >= sd_cutoff))
  npc <- ifelse(npc == 0, npc+1, npc)
  npc <- ifelse(npc == 1, npc+1, npc)
  npc <- ifelse(npc > 10, 10, npc)


  df <- as.data.frame(t(summary(pca)$importance))
  df$index <- 1:nrow(df)
  if(!is.null(max_pc)){
    if(is.numeric(max_pc) & max_pc < nrow(df)){
      df <- df[1:max_pc,]
    }
  }
  colnames(df) <- c("Standard_Deviation","Proportion_of_Variance","Cumulative_Proportion","index")
  dfs <- df[1:20,]
  i_1 <- which(summary(pca)$importance[3,]>=0.8)[1]
  if(i_1==1){
    df_cp=df[i_1:(i_1+5), ]
  }else if(i_1<=5){
    df_cp=df[1:(i_1+10), ]
  }else{
    df_cp=df[(i_1-5):(i_1+10), ]
  }

  rj.ftheme <-   theme(panel.background = element_rect(fill='transparent', color='black'),
                       panel.grid.major=element_blank(),
                       panel.grid.minor=element_blank(),
                       #panel.border = element_rect(fill='transparent', color='transparent'),
                       plot.title = element_text(size = 12), # centers, hjust = 0.5
                       plot.subtitle= element_text(size = 10),
                       #plot.caption = element_text()，
                       legend.key = element_rect( fill = "white"),
                       #axis.title=element_text(size = 12),
                       #axis.ticks=element_blank(),
                       #axis.line = element_line(colour = 'black'),
                       axis.ticks = element_blank(),
                       #axis.text.y = element_blank(),
                       axis.text.x = element_text(vjust=0.5, size = 8) )
  rj.ftheme.small <-   theme(panel.background = element_rect(fill='transparent', color='black'),
                             panel.grid.major=element_blank(),
                             panel.grid.minor=element_blank(),
                             #panel.border = element_rect(fill='transparent', color='transparent'),
                             plot.title = element_text(size = 8), # centers, hjust = 0.5
                             plot.subtitle= element_text(size = 10),
                             #plot.caption = element_text()，
                             legend.key = element_rect( fill = "white"),
                             #axis.title=element_text(size = 12),
                             #axis.ticks=element_blank(),
                             #axis.line = element_line(colour = 'black'),
                             axis.ticks = element_blank(),
                             #axis.text.y = element_blank(),
                             axis.text.x = element_text(vjust=0.5, size = 8)   )

  p1_1 <- ggplot(df, aes(x = .data$index, y = .data$Standard_Deviation)) +
    geom_point(shape = 0, color = "#33681e", size = 2.5, alpha = 0.6) +
    geom_hline(aes(yintercept = sd_cutoff), colour = "#990000", linetype = "dashed") +
    ylim(c(0, max(df$Standard_Deviation))) +
    labs(title = "PC ~ Standard Deviation ",
         subtitle = paste0("Standard Deviation >= ", sd_cutoff, ": PC1 ~ PC", tail(which(df$Standard_Deviation>=1),1) ),
         x = "pc index",
         y = "Standard Deviation") +
    rj.ftheme

  p2_1 <- ggplot(df, aes(x = .data$index, y = .data$Proportion_of_Variance)) +
    geom_point(shape = 1, color = "#6c4c40", size = 2.5, alpha = 0.6) +
    #geom_hline(aes(yintercept = 0.8), colour = "#990000", linetype = "dashed") +
    ylim(c(0, max(df$Proportion_of_Variance))) +
    labs(title = "PC ~ Proportion of Variance ",
         #subtitle = "cutoff: 1",
         x = "pc index",
         y = "Proportion of Variance") +
    rj.ftheme

  p3_1 <- ggplot(df, aes(x = .data$index, y = .data$Cumulative_Proportion)) +
    geom_point(shape = 2, color = "#00579a", size = 2.5, alpha = 0.6) +
    geom_hline(aes(yintercept = 0.8), colour = "#990000", linetype = "dashed") +
    ylim(c(0, max(df$Cumulative_Proportion))) +
    labs(title = "PC Number ~ Cumulative Proportion ",
         subtitle = paste0("Cumulative Proportion >= 0.8: ",
                           rownames(df)[which(df$Cumulative_Proportion>=0.8)[1]]),
         x = "pc number",
         y = "Cumulative Proportion") +
    rj.ftheme

  #
  p1_2 <- ggplot(dfs, aes(x = .data$index, y = .data$Standard_Deviation)) +
    geom_point(shape = 3, color = "sienna1", size = 1.5, alpha = 0.7) +
    geom_hline(aes(yintercept = sd_cutoff), colour = "#990000", linetype = "dashed") +
    labs(title = paste0(rownames(dfs)[1]," - ",rownames(dfs)[nrow(dfs)]), x = "", y = "") +
    rj.ftheme.small

  p2_2 <- ggplot(dfs, aes(x = .data$index, y = .data$Proportion_of_Variance)) +
    geom_point(shape = 3, color = "sienna1", size = 1.5, alpha = 0.7) +
    #geom_hline(aes(yintercept = 1), colour = "#990000", linetype = "dashed") +
    labs(title = paste0(rownames(dfs)[1]," - ",rownames(dfs)[nrow(dfs)]), x = "", y = "") +
    rj.ftheme.small

  p3_2 <- ggplot(df_cp, aes(x = .data$index, y = .data$Cumulative_Proportion)) +
    geom_point(shape = 3, color = "sienna1", size = 1.5, alpha = 0.7) +
    geom_hline(aes(yintercept = 0.8), colour = "#990000", linetype = "dashed") +
    labs(title = paste0(rownames(df_cp)[1]," - ",rownames(df_cp)[nrow(df_cp)]),
         x = "", y = "") +
    rj.ftheme.small

  g1 <- ggplotGrob(p1_2)

  if(max(df$Standard_Deviation)/4*2.4 < sd_cutoff & max(df$Standard_Deviation) > sd_cutoff){
    p1 <- p1_1 + annotation_custom(g1,
                                   xmin = nrow(df)/4*2, xmax = nrow(df),
                                   ymin = max(df$Standard_Deviation)/4*0.3, ymax = max(df$Standard_Deviation)/4*2.4)
  }else {
    p1 <- p1_1 + annotation_custom(g1,
                                   xmin = nrow(df)/4*2, xmax = nrow(df),
                                   ymin = max(df$Standard_Deviation)/4*1.7, ymax = max(df$Standard_Deviation))
  }


  g2 <- ggplotGrob(p2_2)
  p2 <- p2_1 + annotation_custom(g2,
                                 xmin = nrow(df)/4*1.9, xmax = nrow(df),
                                 ymin = max(df$Proportion_of_Variance)/4*1.8, ymax = max(df$Proportion_of_Variance))

  g3 <- ggplotGrob(p3_2)
  p3 <- p3_1 + annotation_custom(g3,
                                 xmin = nrow(df)/4*1.9, xmax = nrow(df),
                                 ymin = 0, ymax = max(df$Cumulative_Proportion)/4*2.2)

  p <- list(p1,p2,p3)

  return(p)
}

#' PlotLabelSplit
#'
#' @description
#' plot label in tSNE/UMAP/DM/PCA split
#'
#' @param object MGPfact object
#' @param labels attributes in metdata
#' @param width canvas width
#'
#' @export
#'
PlotLabelSplit <- function(object = NULL,
                           labels = NULL,
                           width = 9){

  pca = is.null(object@DimReducs@PCA)
  tsne = is.null(object@DimReducs@tSNE)
  umap = is.null(object@DimReducs@UMAP)
  dm = is.null(object@DimReducs@DM)
  TTAG = c("PCA" = pca, "tSNE" = tsne, "UMAP" = umap, "DM" = dm)
  stag = names(which(!TTAG))

  ### 分别在四个维度上可视化不同的label
  len = length(labels)
  tmp = data.frame(object@MetaData[,labels])
  colnames(tmp) = labels
  for(i in 1:ncol(tmp)) tmp[,i] = as.factor(tmp[,i])
  df = tmp
  for(c in stag){
    if(c=="PCA"){
      df <- data.frame(PCA1 = object@DimReducs@PCA$x[,1],
                       PCA2 = object@DimReducs@PCA$x[,2],
                       df, stringsAsFactors = FALSE)
    }
    if(c=="tSNE"){
      df <- data.frame(tSNE1 = object@DimReducs@tSNE$Y[,1],
                       tSNE2 = object@DimReducs@tSNE$Y[,2],
                       df, stringsAsFactors = FALSE)
    }
    if(c=="UMAP"){
      df <- data.frame(UMAP1 = object@DimReducs@UMAP$layout[,1],
                       UMAP2 = object@DimReducs@UMAP$layout[,2],
                       df, stringsAsFactors = FALSE)
    }
    if(c=="DM"){
      df <- data.frame(DM1 = object@DimReducs@DM$X[,1],
                       DM2 = object@DimReducs@DM$X[,2],
                       df, stringsAsFactors = FALSE)
    }
  }
  # df <- data.frame(df, tmp, stringsAsFactors = FALSE)

  dir1 = "1_murp/reduction_split/"
  dir.create(dir1)

  for(c in stag){

    cat("Plot: ", c, "\n")
    a = paste0(c,1)
    b = paste0(c,2)

    fat = ifelse(len==1, 5, 8)
    ### 1. all point
    plist <- lapply(labels, function(lab){
      ggplot(df, aes_string(x = a, y = b, colour = lab) ) +
        geom_point(alpha = 0.6, size = 1.5) +
        scale_colour_d3("category20") +
        guides(colour = guide_legend(override.aes=list(size=4))) +
        rj.ftheme
    })
    p <- wrap_plots(plist, ncol = 1)
    ggplot2::ggsave(paste0(dir1, c,"_celltype.pdf"), p, width = width, height = ceiling(len/2)*fat, limitsize = FALSE)
    # ggplot2::ggsave(paste0(dir1, c,"_celltype.png"), p, width = width, height = ceiling(len/2)*8, limitsize = FALSE)

    ### 2. murp
    df1 <- df
    tmp <- tapply(1:nrow(df1),
                  object@MURP$Recommended_K_cl$cluster,
                  function(x, y){
                    apply(df1[x, c(a, b)],2,mean)
                  })
    df2 <- do.call(rbind, tmp) %>% as.data.frame
    tmp <- data.frame(object@MURP$murp_cellinfo[,labels])
    colnames(tmp) <- labels
    df2 <- cbind(df2, tmp)
    df2$label <- 1:nrow(df2)
    plist <- lapply(labels, function(lab){
      ggplot() +
        geom_point(data = df1, aes_string(x = a, y = b),
                   colour = "grey80", alpha = 0.2, size = 2.3) +
        geom_point(data = df2, aes_string(x = a, y = b, colour = lab),
                   alpha = 0.6, size = 2.3) +
        scale_colour_d3("category20") +
        guides(colour = guide_legend(override.aes=list(size=4))) +
        rj.ftheme
    })
    p <- wrap_plots(plist, ncol = 1)
    ggplot2::ggsave(paste0(dir1, c, "_murp_celltype.pdf"), p,
                    width = width, height = ceiling(len/2)*fat, limitsize = FALSE)
    # ggplot2::ggsave(paste0(dir1, c, "_murp_celltype.png"), p, width = width, height = ceiling(len/2)*8, limitsize = FALSE)

    ### 3. label
    plist <- lapply(labels, function(lab){
      ggplot() +
        geom_point(data = df1, aes_string(x = a, y = b),
                   colour = "grey80", alpha = 0.2, size = 2.3) +
        geom_point(data = df2, aes_string(x = a, y = b, colour = lab),
                   alpha = 0.6, size = 2.3) +
        geom_text(data = df2, aes_string(x = a, y = b, label = "label"),
                  size = 2, check_overlap = FALSE) +
        scale_colour_d3("category20") +
        guides(colour = guide_legend(override.aes=list(size=4))) +
        rj.ftheme
    })
    p <- wrap_plots(plist, ncol = 1)
    ggplot2::ggsave(paste0(dir1, c, "_murp_label.pdf"), p,
                    device = cairo_pdf,
                    width = width, height = ceiling(len/2)*fat, limitsize = FALSE)
    # ggplot2::ggsave(paste0(dir1, c, "_murp_label.png"), p, width = width, height = ceiling(len/2)*8, limitsize = FALSE)

  }

}


#' PlotLabelMerge
#'
#' @description
#' plot label in tSNE/UMAP/DM/PCA merge
#'
#' @param object MGPfact object
#' @param labels attributes in metdata
#' @param width canvas width
#'
#' @export
#'
PlotLabelMerge <- function(object = NULL,
                           labels = NULL,
                           width = 17){

  pca = is.null(object@DimReducs@PCA)
  tsne = is.null(object@DimReducs@tSNE)
  umap = is.null(object@DimReducs@UMAP)
  dm = is.null(object@DimReducs@DM)
  TTAG = c("PCA" = pca, "tSNE" = tsne, "UMAP" = umap, "DM" = dm)
  stag = names(which(!TTAG))

  ### 分别在四个维度上可视化不同的label
  len = length(labels)
  tmp = data.frame(object@MetaData[,labels])
  colnames(tmp) = labels
  for(i in 1:ncol(tmp)) tmp[,i] = as.factor(tmp[,i])
  df = tmp
  for(c in stag){
    if(c=="PCA"){
      df <- data.frame(PCA1 = object@DimReducs@PCA$x[,1],
                       PCA2 = object@DimReducs@PCA$x[,2],
                       df, stringsAsFactors = FALSE)
    }
    if(c=="tSNE"){
      df <- data.frame(tSNE1 = object@DimReducs@tSNE$Y[,1],
                       tSNE2 = object@DimReducs@tSNE$Y[,2],
                       df, stringsAsFactors = FALSE)
    }
    if(c=="UMAP"){
      df <- data.frame(UMAP1 = object@DimReducs@UMAP$layout[,1],
                       UMAP2 = object@DimReducs@UMAP$layout[,2],
                       df, stringsAsFactors = FALSE)
    }
    if(c=="DM"){
      df <- data.frame(DM1 = object@DimReducs@DM$X[,1],
                       DM2 = object@DimReducs@DM$X[,2],
                       df, stringsAsFactors = FALSE)
    }
  }

  dir2 = "1_murp/reduction_merge/"
  dir.create(dir2)

  ht = ifelse(length(stag)<3, 4, 8)
  ## 四个维度合并
  for(lab in labels){

    cat("Plot: ", lab, "\n")
    if(is.null(width)){
      d2_levels = unique(tmp[,lab])
      width = 9.5 + max(strwidth(d2_levels, "inches"))
    }
    ### 1. all point
    plist <- lapply(stag, function(c){
      cat(c, "\n")
      a = paste0(c,1)
      b = paste0(c,2)
      p1 <- ggplot(df, aes_string(x = a, y = b, colour = lab)) +
        geom_point(alpha = 0.6, size = 1.5) +
        scale_colour_d3("category20") +
        guides(colour = guide_legend(override.aes=list(size=4))) +
        rj.ftheme
    })
    p <- wrap_plots(plist, ncol = 2, guides = "collect")
    ggplot2::ggsave(paste0(dir2, "1_", lab,".pdf"), p, device = cairo_pdf, width = width, height = ht)
    # ggplot2::ggsave(paste0(dir2, "1_", lab,".png"), p, width = width, height = 8)

    ### 2. murp
    plist <- lapply(stag, function(c){
      cat(c, "\n")
      a = paste0(c,1)
      b = paste0(c,2)

      df1 <- df
      tmp <- tapply(1:nrow(df1),
                    object@MURP$Recommended_K_cl$cluster,
                    function(x, y){
                      apply(df1[x, c(a, b)],2,mean)
                    })
      df2 <- do.call(rbind, tmp) %>% as.data.frame
      tmp <- data.frame(object@MURP$murp_cellinfo[,labels])
      colnames(tmp) <- labels
      df2 <- cbind(df2, tmp)
      df2$label <- 1:nrow(df2)
      # df2$seurat_clusters <- factor(df2$seurat_clusters,levels = c(0,1,4,6,7,8,10))
      ggplot() +
        geom_point(data = df1, aes_string(x = a, y = b),
                   colour = "grey80", alpha = 0.2, size = 2.3) +
        geom_point(data = df2, aes_string(x = a, y = b, colour = lab),
                   alpha = 0.6, size = 2.3) +
        scale_colour_d3("category20") +
        guides(colour = guide_legend(override.aes=list(size=4))) +
        rj.ftheme
    })
    p <- wrap_plots(plist, ncol = 2, guides = "collect")
    ggplot2::ggsave(paste0(dir2, "2_murp_", lab,".pdf"), p, device = cairo_pdf, width = width, height = ht)
    # ggplot2::ggsave(paste0(dir2, "2_murp_", lab,".png"), p, width = width, height = 8)

    ### 3. murp_label
    plist <- lapply(stag, function(c){
      cat(c, "\n")
      a = paste0(c,1)
      b = paste0(c,2)

      df1 <- df
      tmp <- tapply(1:nrow(df1),
                    object@MURP$Recommended_K_cl$cluster,
                    function(x, y){
                      apply(df1[x, c(a, b)],2,mean)
                    })
      df2 <- do.call(rbind, tmp) %>% as.data.frame
      tmp <- data.frame(object@MURP$murp_cellinfo[,labels])
      colnames(tmp) <- labels
      df2 <- cbind(df2, tmp)
      df2$label <- 1:nrow(df2)
      # df2$seurat_clusters <- factor(df2$seurat_clusters,levels = c(0,1,4,6,7,8,10))
      ggplot() +
        geom_point(data = df1, aes_string(x = a, y = b),
                   colour = "grey80", alpha = 0.2, size = 2.3) +
        geom_point(data = df2, aes_string(x = a, y = b, colour = lab),
                   alpha = 0.6, size = 2.3) +
        geom_text(data = df2, aes_string(x = a, y = b, label = "label"),
                  size = 2, check_overlap = FALSE) +
        scale_colour_d3("category20") +
        guides(colour = guide_legend(override.aes=list(size=4))) +
        rj.ftheme
    })
    p <- wrap_plots(plist, ncol = 2, guides = "collect")
    ggplot2::ggsave(paste0(dir2, "2_murp_label_",lab,".pdf"), p, device = cairo_pdf, width = width, height = ht)
    # ggplot2::ggsave(paste0(dir2, "2_murp_label_",lab,".png"), p, width = width, height = 8)
  }
}


#' PlotBinBif
#'
#' @description
#' Visualize bifurcation on the binary-tree
#'
#' @param object MGPfact object
#' @param plot_label logical value, whether to plot label
#' @export
PlotBinBif <- function(object,
                       plot_label = FALSE){

  pointColorValue <- colorRampPalette(pal_d3("category10")(10))(10)[c(8,1,4)]
  names(pointColorValue) <- c(0,1,2)

  bintree_all = GetBinTreeResult(object)
  P = getParams(object, "murp_number")
  L = getParams(object, "trajectory_number")
  Q = getParams(object, "murp_pc_number")

  plist <- lapply(1:L, function(i){

    layout = bintree_all[[i]]$layout
    graph = bintree_all[[i]]$graph

    p <- ggraph(graph, layout = layout) +
      geom_edge_diagonal(alpha = 0.7, width = 0.5, check_overlap = FALSE) +
      geom_node_point(aes(color = factor(.data$c_label, levels = 0:2) ), size = 6, alpha = 0.6) +
      scale_color_manual(values = pointColorValue, breaks = c(0:2), labels = paste0("Phase ",0:2), drop = F) +
      coord_flip() +
      scale_y_reverse() +
      labs(title = paste0("Trajectory ",i), color = "") +
      rj.graph.ftheme

    # 判断phase1是否在上面
    tm = bintree_all[[i]]$vertex[,paste0("C0_",i)]
    in1 = which(tm==1)
    in2 = which(tm==2)

    if( (length(in1)!=0) & (length(in2)!=0) ){
      if(layout[in1[1],1]<layout[in2[1],1]){
        p = p +
          scale_x_reverse()
      }
    }

    if(plot_label){
      p + geom_node_text(aes(label = .data$label), size = 4)
    }else{
      p
    }
  })
  p <- patchwork::wrap_plots(plist, nrow = 3, guides = "collect")
  ggsave(paste0("2_pseudotime/2.2_binarytree/1_binarytree_bif.pdf"), p,
         width = 12, height = 9)
}

#' PlotBinLabel
#'
#' @description
#' Visualize cell property labels on the binary-tree
#'
#' @param object MGPfact object
#' @param plot_label logical value, whether to plot label
#' @param labels cell property labels
#' @export
PlotBinLabel <- function(object,
                         plot_label = FALSE,
                         labels = NULL){


  sdf = GetMURPInfo(object)
  bintree_all = GetBinTreeResult(object)
  P = getParams(object, "murp_number")
  L = getParams(object, "trajectory_number")
  Q = getParams(object, "murp_pc_number")

  for(i in 1:length(labels) ){
    lab = labels[i]
    cat(lab, "\n")
    plist <- lapply(1:L, function(j){

      layout = bintree_all[[j]]$layout
      graph = bintree_all[[j]]$graph
      if(is.null(levels(sdf[,lab]))){
        levs = unique(sdf[,lab])
      }else{
        levs = levels(sdf[,lab])
      }
      p <- ggraph(graph, layout = layout) +
        geom_edge_diagonal(alpha = 0.7, width = 0.5, check_overlap = FALSE) +
        geom_node_point(aes(colour = factor(get(lab),levels = levs) ),
                        size = 6, alpha = 0.6) +
        scale_colour_d3("category20", breaks = levs, labels = levs) +
        coord_flip() +
        scale_y_reverse() +
        labs(title = paste0("L",j), colour = lab) +
        rj.graph.ftheme

      # 判断phase1是否在上面
      tm = bintree_all[[j]]$vertex[,paste0("C0_",j)]
      in1 = which(tm==1)
      in2 = which(tm==2)
      if( (length(in1)!=0) & (length(in2)!=0) ){
        if(layout[in1[1],1]<layout[in2[1],1]){
          p = p +
            scale_x_reverse()
        }
      }

      if(plot_label){
        p + geom_node_text(aes(label = .data$label), size = 4)
      }else{
        p
      }
    })
    p <- patchwork::wrap_plots(plist, nrow = L, guides = "collect")
    ggsave(paste0("2_pseudotime/2.2_binarytree/2_binarytree_",lab,".pdf"),
           device = cairo_pdf, p, width = 17, height = 9)
  }
  invisible()

}


#' PlotTbPse
#'
#' @description
#' Visualize pseudotime on the binary-tree
#'
#' @param object MGPfact object
#' @param plot_label logical value, whether to plot label
#' @export
PlotTbPse <- function(object,
                      plot_label = FALSE){

  sdf = GetMURPInfo(object)
  tbtree = GetTbTreeResult(object)
  P = getParams(object, "murp_number")
  L = getParams(object, "trajectory_number")
  Q = getParams(object, "murp_pc_number")

  p <- ggraph(tbtree$graph) +
    geom_edge_diagonal(alpha = 0.7, width = 0.5, check_overlap = FALSE) +
    geom_node_point(aes(colour = .data$T ), size = 10, alpha = 0.8) +
    scale_colour_material(name = "pse", palette = "deep-orange") +
    coord_flip() +
    scale_y_reverse() +
    rj.graph.ftheme
  if(plot_label){
    p <- p + geom_node_text(aes(label = .data$label), size = 4)
  }
  ggsave(paste0("2_pseudotime/2.3_tbtree/1_tbtree_label_", plot_label, "_pse.pdf"), p, width = 15, height = 5, limitsize = FALSE)
}

#' PlotTbLabel
#'
#' @description
#' Visualize cell property labels on the tb-tree
#'
#' @param object MGPfact object
#' @param labels cell property labels
#' @export
PlotTbLabel <- function(object,
                        labels = NULL){

  sdf = GetMURPInfo(object)
  tbtree = GetTbTreeResult(object)
  P = getParams(object, "murp_number")
  L = getParams(object, "trajectory_number")
  Q = getParams(object, "murp_pc_number")

  if(!file.exists("2_pseudotime/2.3_tbtree")){
    dir.create("2_pseudotime/2.3_tbtree", recursive = TRUE)
  }
  for(i in 1:length(labels)){
    lab = labels[i]
    cat(lab, "\n")
    if(is.null(levels(sdf[,lab]))){
      levs = unique(sdf[,lab])
    }else{
      levs = levels(sdf[,lab])
    }
    p <- ggraph(tbtree$graph) +
      geom_edge_diagonal(alpha = 0.7, width = 0.5, check_overlap = FALSE) +
      geom_node_point(aes(colour = factor(get(lab),levels = levs) ),
                      size = 8, alpha = 0.8) +
      coord_flip() +
      scale_y_reverse() +
      scale_colour_d3("category20") +
      labs(colour = lab) +
      rj.graph.ftheme
    ggsave(paste0("2_pseudotime/2.3_tbtree/2_tbtree_",lab,".pdf"),
           device = cairo_pdf, p, width = 15, height = 5, limitsize = FALSE)
  }
}


#' PlotTbBif
#'
#' @description
#' Visualize the branching of independent differentiation trajectories on the tb-tree
#'
#' @param object MGPfact object
#' @param plot_label logical value, whether to plot label
#' @export
PlotTbBif <- function(object,
                      plot_label = FALSE){

  pointColorValue <- colorRampPalette(pal_d3("category10")(10))(10)[c(8,1,4)]
  names(pointColorValue) <- c(0,1,2)

  sdf = GetMURPInfo(object)
  tbtree = GetTbTreeResult(object)
  P = getParams(object, "murp_number")
  L = getParams(object, "trajectory_number")
  Q = getParams(object, "murp_pc_number")

  layout = tbtree$layout
  graph = tbtree$graph
  plist <- lapply(1:L, function(i){

    c_label = paste0("C0_",i)
    p <- ggraph(graph, layout = layout) +
      geom_edge_diagonal(alpha = 0.7, width = 0.5, check_overlap = FALSE) +
      geom_node_point(aes(color = get(c_label) ), size = 6, alpha = 0.6) +
      scale_color_manual(values = pointColorValue) +
      coord_flip() +
      scale_y_reverse() +
      labs(title = paste0("L",i), color = "Bif") +
      rj.graph.ftheme

    if(plot_label){
      p + geom_node_text(aes(label = .data$label), size = 4)
    }else{
      p
    }
  })
  p <- patchwork::wrap_plots(plist, nrow = 3, guides = "collect")
  ggsave(paste0("2_pseudotime/2.3_tbtree/1_tbtree_bif.pdf"),
         device = cairo_pdf, p, width = 12, height = 9)
}


#' PlotTbMainPse
#'
#' @description
#' plot pseudotime on the main path of the consensus trajectory
#'
#' @param object MGPfact object
#' @param col color
#' @export
PlotTbMainPse <- function(object, col = "blue"){

  tbtree_all = GetTbTreeAllResult(object)

  graph <- tbtree_all$graph
  edge <- tbtree_all$edge
  vertex <- tbtree_all$vertex
  start <- tbtree_all$start
  fork_p <- tbtree_all$fork_p
  bb <- tbtree_all$bb

  ### 15. 主点 + 主路线
  p <- ggraph(graph, layout = "manual", x = bb$xy[, 1], y = bb$xy[, 2]) +
    geom_edge_link0(aes(filter = .data$weight==1), width = 0.5, alpha = 0.2) +
    geom_node_point(aes(filter = .data$alpha_murp!=0, size = .data$rat, colour = .data$size), alpha = 0.8, shape = 16) +
    geom_node_label(aes(filter = .data$name %in% start, label = "root"), fill = "#9bacb9", colour = "white", repel = FALSE) +
    geom_node_label(aes(filter = .data$name %in% fork_p, label = "fork"), fill = "#9bacb9", colour = "white", repel = FALSE) +
    # geom_node_point(aes(filter = name %in% fork_p), size = 6, colour = "red", alpha = 1, shape = 17) +
    # geom_node_point(aes(filter = name %in% start, size = rat), colour = "blue", alpha = 1, shape = 17) +
    # scale_colour_gradient(name = "PseudoT", low = scol[[col]][5], high = scol[[col]][10]) +
    scale_color_material(palette = col) +
    scale_size_continuous(name = "Cell Ratio") +
    scale_edge_width(range = c(0.05,0.6)) +
    labs(color = "Pseudotime") +
    guides(alpha="none") +
    rj.graph.ftheme
  ggsave(paste0("2_pseudotime/2.3_tbtree/3_main_pse.pdf"), p, width = 12, height = 5)
}

#' PlotTbMainLabel
#'
#' @description
#' Visualize cell property labels on the main path of the consensus trajectory
#'
#' @param object MGPfact object
#' @param labels cell property labels
#' @export
PlotTbMainLabel <- function(object, labels = NULL){

  sdf = GetMURPInfo(object)
  tbtree_all = GetTbTreeAllResult(object)

  graph <- tbtree_all$graph
  edge <- tbtree_all$edge
  vertex <- tbtree_all$vertex
  start <- tbtree_all$start
  fork_p <- tbtree_all$fork_p
  bb <- tbtree_all$bb

  ### 15. 主点 + 主路线 + celltype
  for(i in 1:length(labels)){
    lab = labels[i]
    cat(lab, "\n")
    if(is.null(levels(sdf[,lab]))){
      levs = unique(sdf[,lab])
    }else{
      levs = levels(sdf[,lab])
    }
    p <- ggraph(graph, layout = "manual", x = bb$xy[, 1], y = bb$xy[, 2]) +
      geom_edge_link0(aes(filter = .data$weight==1), width = 0.5, alpha = 0.1) +
      #geom_node_point(aes(filter = alpha_murp==0, size = rat), alpha = 0.7, colour = scol[[col]][4], shape = 16) +
      geom_node_point(aes(filter = .data$alpha_murp!=0, size = .data$rat,
                          colour = factor(get(lab),levels = levs) ),
                      alpha = 0.8, shape = 16) +
      geom_node_label(aes(filter = .data$name %in% start, label = "root"), repel = FALSE) +
      geom_node_label(aes(filter = .data$name %in% fork_p, label = "fork"), repel = FALSE) +
      # geom_node_point(aes(filter = name %in% fork_p), size = 8, colour = "black", alpha = 1, shape = 8) +
      # geom_node_point(aes(filter = name %in% start, size = rat), size = 8, colour = "black", alpha = 1, shape = 14) +
      # scale_colour_gradient(name = "PseudoT", low = scol[[col]][5], high = scol[[col]][10]) +
      scale_colour_d3("category20") +
      scale_size_continuous(range = c(1,10)) +
      scale_edge_width(range = c(0.05,0.6)) +
      guides(alpha="none", colour = guide_legend(override.aes=list(size=4))) +
      labs(colour = lab, size = "Cell Ratio") +
      rj.graph.ftheme
    ggsave(paste0("2_pseudotime/2.3_tbtree/3_main_",lab,".pdf"),
           device = cairo_pdf, p, width = 13, height = 5)
  }
}

#' PlotTbAllLabel
#'
#' @description
#' Visualize cell property labels on the consensus trajectory
#'
#' @param object MGPfact object
#' @param labels cell property labels
#' @export
PlotTbAllLabel <- function(object, labels = NULL){

  tbtree_all = GetTbTreeAllResult(object)
  sdf = GetMURPInfo(object)

  graph <- tbtree_all$graph
  edge <- tbtree_all$edge
  vertex <- tbtree_all$vertex
  start <- tbtree_all$start
  fork_p <- tbtree_all$fork_p
  bb <- tbtree_all$bb

  ### 主点 + 主路线 + celltype
  for(i in 1:length(labels)){
    lab = labels[i]
    cat(lab, "\n")
    if(is.null(levels(sdf[,lab]))){
      levs = unique(sdf[,lab])
    }else{
      levs = levels(sdf[,lab])
    }

    # cat(lab, "\n")
    p <- ggraph(graph, layout = "manual", x = bb$xy[, 1], y = bb$xy[, 2]) +
      geom_edge_link0(colour = "#202020", width = 0.1, alpha = 0.1) +
      geom_node_point(aes(colour = factor(get(lab),levels = levs),
                          size = .data$rat, alpha = .data$alpha)) +
      geom_node_label(aes(filter = .data$name %in% start, label = "root"), repel = FALSE) +
      geom_node_label(aes(filter = .data$name %in% fork_p, label = "fork"), repel = FALSE) +
      # geom_node_point(aes(filter = name %in% fork_p), size = 6, colour = "red", alpha = 1, shape = 17) +
      # geom_node_point(aes(filter = name %in% start, size = rat), colour = "blue", alpha = 1, shape = 17) +
      # scale_edge_color_gradient(low = scol[["grey"]][5], high = scol[["grey"]][10]) +
      scale_size_continuous(name = "Cell Ratio", range = c(1,5)) +
      scale_alpha_continuous(range = c(0.1,0.9)) +
      scale_colour_d3("category20") +
      labs(colour = lab) +
      guides(colour = guide_legend(override.aes=list(size=4)), alpha = "none") +
      rj.graph.ftheme
    ggsave(paste0("2_pseudotime/2.3_tbtree/4_all_",lab,".pdf"),
           device = cairo_pdf, p, width = 13, height = 5)
  }

}

#' rotate_point
#'
#' @description
#' Rotate coordinates counterclockwise.
#' @param x X-axis
#' @param y Y-axis
#' @param angle Rotation angle
#' @export
rotate_point <- function(x, y, angle) {
  rad <- angle * pi / 180  # 将角度转换为弧度
  new_x <- x * cos(rad) - y * sin(rad)
  new_y <- x * sin(rad) + y * cos(rad)
  return(c(new_x, new_y))
}

#' pie plot for every row of data.frame
#'
#' @description
#' Pie chart composition
#'
#' @param i Extract the attributes of the i-th column
#' @param tab DataFrame related to the property
#' @param levs facter levels
#' @param size point size
#' @param color point outline color
#' @param col_values point fill color
plot_pie <- function(i, tab, levs, size = 0.2, color = NA, col_values = NULL) {

  if(is.null(col_values)){
    col_values = setNames(pal_d3("category20")(length(levs))[seq_along(levs)],levs)
  }

  df = melt(tab[i,,drop = FALSE],"name")
  df$variable = factor(df$variable, levels = levs)
  if(is.na(color)){
    ggplot(df, aes_string(x = 1, fill = "variable", weight = "value")) +
      geom_bar(size = size, position = "stack") +
      coord_polar(theta = "y", start = 0) +
      scale_fill_manual(values = col_values, limits = levs, breaks = levs, labels = levs, drop=FALSE ) +
      theme_void() +
      theme(legend.position = "none")
  }else{
    ggplot(df, aes_string(x = 1, fill = "variable", weight = "value")) +
      geom_bar(size = size, position = "stack",color=color) +
      coord_polar(theta = "y", start = 0) +
      scale_fill_manual(values = col_values, limits = levs, breaks = levs, labels = levs, drop=FALSE ) +
      theme_void() +
      theme(legend.position = "none")
  }

}

#' PlotPieBinLabel
#'
#' @description
#' Visualize cell property labels on the binary-tree,
#' using the proportion of cell types in the MURP subset to represent the points.
#'
#' @param object MGPfact object
#' @param labels cell property labels
#' @param size point size
#' @param color Point outline color
#' @export
PlotPieBinLabel <- function(object,
                            labels = NULL,
                            size = 0.01,
                            color = "black"){

  sdf = GetMURPInfo(object)
  bintree_all = GetBinTreeResult(object)
  P = getParams(object, "murp_number")
  L = getParams(object, "trajectory_number")
  Q = getParams(object, "murp_pc_number")


  for(i in 1:length(labels) ){
    lab = labels[i]
    levs = unique(sdf[,lab])
    cat(lab, "\n")

    murp_c = paste0("T[",object$murp_cluster,"]")
    tab = table(object@MetaData[,lab],murp_c) %>% as.matrix
    # tab_r = apply(tab, 2, function(x)x/sum(x))
    tab = apply(tab, 2, function(x)x/sum(x)) %>% as.data.frame.matrix %>% t
    # colnames(tab) = gsub(" ","_",colnames(tab))
    # colnames(tab) = gsub("\\+","",colnames(tab))

    plist = list()
    for(j in 1:L){

      cat(j, "\n")

      edge =  bintree_all[[j]]$edge
      vertex = bintree_all[[j]]$vertex
      layout = bintree_all[[j]]$layout
      graph = bintree_all[[j]]$graph

      tab2 = data.frame(tab[rownames(vertex),], name = rownames(tab))
      colnames(tab2)[1:(ncol(tab2)-1)] = colnames(tab)

      ## 逆时针90度
      layout2 <- t(apply(layout, 1, function(row) rotate_point(row[1], row[2], 90)))
      dff <- as.data.frame(layout2)
      dff$pie <- lapply(1:nrow(tab2), plot_pie, tab = tab2, levs = levs,
                        size = size, color = color)

      ## plot
      ggraph(graph, layout = layout2) +
        geom_edge_diagonal(alpha = 0.7, width = 0.5, check_overlap = FALSE) +
        geom_node_point(aes(colour = factor(get(lab),levels = levs) ),
                        size = 6, alpha = 0) +
        geom_subview(data = dff, aes(x = V1, y = V2, subview = pie),
                     # width=1.5, height=1.5) +
                     width=(range(dff$V1)[2] - range(dff$V1)[1])/50,
                     height=(range(dff$V1)[2] - range(dff$V1)[1])/50) +
        scale_colour_d3("category20", breaks = levs, labels = levs) +
        scale_y_continuous(limits = c(-1,1)) +
        labs(title = paste0("Trajectory ",j), colour = lab) +
        guides(colour = guide_legend(override.aes=list(alpha = 1, size=4))) +
        rj.graph.ftheme -> p
      plist = append(plist, list(p))
    }
    p <- patchwork::wrap_plots(plist, nrow = L, guides = "collect")
    ggsave(paste0("2_pseudotime/2.2_binarytree/2_binarytree_pie_",lab,".pdf"), p, width = 17, height = 12)

  }

}

#' PlotPieTbLabel
#'
#' @description
#' Visualize cell property labels on the main path of the tb-tree,
#' using the proportion of cell types in the MURP subset to represent the points.
#'
#' @param object MGPfact object
#' @param labels cell property labels
#' @param size point size
#' @param color Point outline color
#' @param sub_width sub view width
#' @export
PlotPieTbLabel <- function(object,
                           labels = NULL,
                           size = 0.2,
                           color = "black",
                           sub_width = 1){

  sdf = GetMURPInfo(object)
  tbtree = GetTbTreeResult(object)
  P = getParams(object, "murp_number")
  L = getParams(object, "trajectory_number")
  Q = getParams(object, "murp_pc_number")

  edge =  tbtree$edge
  vertex = tbtree$vertex
  layout = tbtree$layout
  graph = tbtree$graph

  for(i in 1:length(labels) ){

    lab = labels[i]
    cat(lab, "\n")
    if(is.null(levels(sdf[,lab]))){
      levs = unique(sdf[,lab])
    }else{
      levs = levels(sdf[,lab])
    }

    ## 得到饼图数据
    murp_c = paste0("T[",object$murp_cluster,"]")
    tab = table(object@MetaData[,lab],murp_c) %>% as.matrix
    tab = apply(tab, 2, function(x)x/sum(x)) %>% as.data.frame.matrix %>% t
    tab2 = data.frame(tab[rownames(vertex),], name = rownames(tab))
    colnames(tab2)[1:(ncol(tab2)-1)] = colnames(tab)

    ## 逆时针90度
    layout2 <- t(apply(layout, 1, function(row) rotate_point(row[1], row[2], 90)))
    dff <- as.data.frame(layout2)
    dff$pie <- lapply(1:nrow(tab2), plot_pie, tab = tab2, levs = levs, size = size, color = color)

    ## plot
    ggraph(graph, layout = layout2) +
      geom_edge_diagonal(alpha = 0.7, width = 0.5, check_overlap = FALSE) +
      geom_node_point(aes(colour = factor(get(lab),levels = levs) ),
                      size = 6, alpha = 0) +
      geom_subview(aes(x = V1, y = V2, subview = pie), width=sub_width, height=sub_width, data = dff) +
      scale_colour_d3("category20", breaks = levs, labels = levs) +
      labs(colour = lab) +
      guides(colour = guide_legend(override.aes=list(alpha = 1, size=4))) +
      rj.graph.ftheme -> p
    ggsave(paste0("2_pseudotime/2.3_tbtree/2_tbtree_pie_",lab,".pdf"),
           p, width = 15, height = 5, limitsize = FALSE)

  }

}

#' PlotPieTbMainLabel
#'
#' @description
#' Visualize cell property labels on the main path of the consensus trajectory,
#' using the proportion of cell types in the MURP subset to represent the points.
#'
#' @param object MGPfact object
#' @param labels cell property labels
#' @param size point size
#' @param color Point outline color
#' @param sub_scale_factor Sub view's size scaling factor
#' @export
PlotPieConsensusMainLabel <- function(object,
                                      labels = NULL,
                                      size = 0.2,
                                      color = "black",
                                      sub_scale_factor = 0.7){

  sdf = GetMURPInfo(object)
  tbtree_all = GetTbTreeAllResult(object)

  graph <- tbtree_all$graph
  edge <- tbtree_all$edge
  vertex <- tbtree_all$vertex
  start <- tbtree_all$start
  fork_p <- tbtree_all$fork_p
  bb <- tbtree_all$bb

  ### 15. 主点 + 主路线 + celltype
  for(i in 1:length(labels)){
    lab = labels[i]
    cat(lab, "\n")
    if(is.null(levels(sdf[,lab]))){
      levs = unique(sdf[,lab])
    }else{
      levs = levels(sdf[,lab])
    }

    ## 得到饼图数据
    murp_c = paste0("T[",object$murp_cluster,"]")
    tab = table(object@MetaData[,lab],murp_c) %>% as.matrix
    tab = apply(tab, 2, function(x)x/sum(x)) %>% as.data.frame.matrix %>% t
    tab2 = data.frame(tab, name = rownames(tab))
    colnames(tab2)[1:(ncol(tab2)-1)] = colnames(tab)

    # og = table(object$cellnames,object$majortype)
    # og = og %>% as.data.frame.matrix
    # og2 = data.frame(og, name = rownames(og))
    # colnames(og2)[1:(ncol(og2)-1)] = colnames(og)
    #
    # tab2 = rbind(og2, tab2)
    # tab2 = tab2[rownames(vertex),]

    ## 逆时针90度
    # layout2 <- t(apply(bb$xy, 1, function(row) rotate_point(row[1], row[2], 90)))
    # dff <- as.data.frame(layout2)
    # dff$pie <- lapply(1:nrow(tab2), plot_pie, tab = tab2, levs)

    V(graph)$x = bb$xy[,1]
    V(graph)$y = bb$xy[,2]
    graph2 = subgraph(graph, vids = paste0("T[",1:nrow(sdf),"]"))

    ## 逆时针90度
    layout2 =  data.frame(V1 = V(graph2)$x, V2 = V(graph2)$y)
    # layout2 <- t(apply(layout_graph2, 1, function(row) rotate_point(row[1], row[2], 90)))
    dff <- as.data.frame(layout2)
    tab2 = tab2[V(graph2)$name,]
    dff$pie <- lapply(1:nrow(tab2), plot_pie, tab = tab2, levs = levs, size = size, color = color)
    dff$rat = V(graph2)$rat * 100
    dff$rat[dff$rat<0.5] = 1
    dff$rat[dff$rat>3] = 3

    ## plot
    dff$rat = dff$rat * sub_scale_factor
    ggraph(graph2, layout = "manual", x = V(graph2)$x, y = V(graph2)$y) +
      geom_edge_link0(width = 0.5, alpha = 0.1) +
      geom_node_point(aes(size = .data$rat, colour = factor(get(lab),levels = levs) ),
                      alpha = 0, shape = 16) +

      geom_subview(aes(x = V1, y = V2, subview = pie, width=rat, height=rat), data = dff) +

      geom_node_label(aes(filter = .data$name %in% start, label = "root"), repel = FALSE) +
      geom_node_label(aes(filter = .data$name %in% fork_p, label = "fork"), repel = FALSE) +
      scale_colour_d3("category20",breaks = levs, labels = levs) +
      # scale_size_continuous(range = c(1,10)) +
      scale_edge_width(range = c(0.05,0.6)) +
      guides(alpha="none",
             colour = guide_legend(override.aes=list(size=4,alpha=1)),
             size = guide_legend(override.aes=list(alpha = 0.6))) +
      labs(colour = lab, size = "Cell Ratio") +
      rj.graph.ftheme -> p
    ggsave(paste0("2_pseudotime/2.3_tbtree/3_main_pie_",lab,".pdf"), p, width = 13, height = 5)

  }
}

#' PlotPieTbAllLabel
#'
#' @description
#' Visualize cell property labels on the main path of the consensus trajectory,
#' using the proportion of cell types in the MURP subset to represent the points.
#'
#' @param object MGPfact object
#' @param labels cell property labels
#' @param size pie size
#' @param color Point outline color
#' @param sub_scale_factor Sub view's size scaling factor
#' @param cell_size point size
#' @param cell_alpha point alpha
#' @export
PlotPieConsensusAllLabel <- function(object,
                                     labels = NULL,
                                     size = 0.2,
                                     color = "black",
                                     sub_scale_factor = 0.7,
                                     cell_size = 1,
                                     cell_alpha = 0.5
                                     ){

  sdf = GetMURPInfo(object)
  tbtree_all = GetTbTreeAllResult(object)

  graph <- tbtree_all$graph
  edge <- tbtree_all$edge
  vertex <- tbtree_all$vertex
  start <- tbtree_all$start
  fork_p <- tbtree_all$fork_p
  bb <- tbtree_all$bb

  ### 15. 主点 + 主路线 + celltype
  for(i in 1:length(labels)){
    lab = labels[i]
    cat(lab, "\n")
    if(is.null(levels(sdf[,lab]))){
      levs = unique(sdf[,lab])
    }else{
      levs = levels(sdf[,lab])
    }

    ## 得到饼图数据
    murp_c = paste0("T[",object$murp_cluster,"]")
    tab = table(object@MetaData[,lab],murp_c) %>% as.matrix
    tab = apply(tab, 2, function(x)x/sum(x)) %>% as.data.frame.matrix %>% t
    tab2 = data.frame(tab, name = rownames(tab))
    colnames(tab2)[1:(ncol(tab2)-1)] = colnames(tab)

    # og = table(object$cellnames,object$majortype)
    # og = og %>% as.data.frame.matrix
    # og2 = data.frame(og, name = rownames(og))
    # colnames(og2)[1:(ncol(og2)-1)] = colnames(og)
    #
    # tab2 = rbind(og2, tab2)
    # tab2 = tab2[rownames(vertex),]

    ## 逆时针90度
    # layout2 <- t(apply(bb$xy, 1, function(row) rotate_point(row[1], row[2], 90)))
    # dff <- as.data.frame(layout2)
    # dff$pie <- lapply(1:nrow(tab2), plot_pie, tab = tab2, levs)

    V(graph)$x = bb$xy[,1]
    V(graph)$y = bb$xy[,2]
    graph2 = subgraph(graph, vids = paste0("T[",1:nrow(sdf),"]"))

    ## 逆时针90度
    layout2 =  data.frame(V1 = V(graph2)$x, V2 = V(graph2)$y)
    dff <- as.data.frame(layout2)
    tab2 = tab2[V(graph2)$name,]
    dff$pie <- lapply(1:nrow(tab2), plot_pie, tab = tab2, levs = levs, size = size, color = color)
    dff$rat = V(graph2)$rat * 100
    dff$rat[dff$rat<0.5] = 1
    dff$rat[dff$rat>3] = 3

    ## plot
    dff$rat = dff$rat * sub_scale_factor
    ggraph(graph, layout = "manual", x = bb$xy[, 1], y = bb$xy[, 2]) +
      geom_edge_link0(colour = "#202020", width = 0.3, alpha = 0.1) +
      geom_node_point(aes(filter = .data$alpha_murp!=1, colour = factor(get(lab),levels = levs)),
                      size = cell_size, alpha = cell_alpha ) +
      geom_node_point(aes(filter = .data$alpha_murp==1, colour = factor(get(lab),levels = levs),
                          size = .data$rat,alpha = 0)) +

      geom_subview(aes(x = V1, y = V2, subview = pie, width=rat, height=rat), data = dff) +

      geom_node_label(aes(filter = .data$name %in% start, label = "root"), repel = FALSE) +
      geom_node_label(aes(filter = .data$name %in% fork_p, label = "fork"), repel = FALSE) +
      scale_size_continuous(name = "Cell Ratio", range = c(1,5)) +
      scale_alpha_continuous(range = c(0.1,0.9)) +
      scale_colour_d3("category20") +
      labs(colour = lab) +
      guides(colour = guide_legend(override.aes=list(size=4, alpha = 0.6)),
             size = guide_legend(override.aes=list(alpha = 0.6)),
             alpha = "none") +
      rj.graph.ftheme -> p
    ggsave(paste0("2_pseudotime/2.3_tbtree/4_all_pie_",lab,".pdf"), p, width = 13, height = 5)

  }
}

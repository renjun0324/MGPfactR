
#' split_tb
#'
#' @description cut T (pseudotime) according to tb (bifurcation time)
#'
#' @param t pseudotime of murps
#' @param tb bifurcation time of different trajectory
#'
#' @return a matrix with Gene Number * Cell Number
#' @noRd
#'
split_tb <- function(t, tb){
  pos = lapply(1:(length(tb)+1), function(i){
    if(i==1){
      which(t < tb[i])
    }else if(i==(length(tb)+1)){
      which(t >= tb[i-1])
    }else {
      which((t >= tb[i-1]) & (t < tb[i]) )
    }
  })
  return(pos)
}


#' split_binary_tree
#'
#' @description
#' 对拆分之后的前后分组再根据分叉情况分组
#' @param df result of 'split_tb'
#' @param sep trajectory name
#'
#' @return list
#' @noRd
split_binary_tree <- function(df, sep = "L1"){

  final_result <- lapply(1:length(df), function(i){
    # print(i)
    if(i == 1){
      x = df[[1]] %>%
        rownames_to_column("c_names") %>%
        select("c_names", T) %>%
        group_split()
      names(x) = "g1"
    }else{
      c = sep
      x = df[[i]] %>% rownames_to_column("c_names") %>%
        select("c_names", T, c) %>%
        group_by(across(all_of(c)) ) %>%
        group_split()
      tmp = do.call(rbind,lapply(x, function(y)y[1,3:ncol(y)]))
      if(!is.null(tmp)){
        names(x) = apply(tmp,1,function(x)paste0("g",i,"_",paste0(x,collapse = "_")))
      }
    }
    x
  })

  return(final_result)
}

#' split_tb_tree
#'
#' @description  根据Tb进行分段之后，对每一个分段再进行不同分支的分组
#' 每一个分段的分组依据是当前tb以及之前tb的L的排列组合
#' @param df_list result of 'split_binary_tree'
#' @return list
#' @noRd
#'
split_tb_tree <- function(df_list){

  tree <- lapply(1:length(df_list), function(i){
    # print(i)
    df = df_list[[i]]

    if(nrow(df)==0){
      x = NULL
    }else if(i == 1){
      x = df_list[[1]] %>%
        rownames_to_column("c_names") %>%
        select("c_names", T) %>%
        group_split()
      names(x) = "g1"
    }else{
      # c = paste0("L",1:(i-1))
      c = tail(colnames(df), ncol(df)-1)
      x = df_list[[i]] %>% rownames_to_column("c_names") %>%
        select("c_names", T, c) %>%
        group_by(across(all_of(c)) ) %>%
        group_split()
      tmp = do.call(rbind,lapply(x, function(y)y[1,3:ncol(y)]))
      names(x) = apply(tmp,1,function(x) paste0("g",i,"_",paste0(x,collapse = "_")))
    }
    x
  })

  return(tree)
}


#' adj_Tb_tree / adj_binary_tree
#'
#' @description
#' adjacency matrix of binary tree
#' @param tb_group_list result of 'split_tb_tree'
#' @param P node number
#' @param point_names murp names
#'
#' @return matrix
#'
#' @noRd
adj_tb_tree <- function(tb_group_list, P, point_names = NULL){
  Group = tb_group_list

  # construct adj matrix
  adj_matrix <- matrix(0, nrow =  P, ncol =  P)
  if(is.null(point_names)){
    rownames(adj_matrix) <- paste0("c",1:P)
    colnames(adj_matrix) <- paste0("c",1:P)
  }else{
    rownames(adj_matrix) <- point_names
    colnames(adj_matrix) <- point_names
  }


  # 在按照tb分割后的每一段分别进行循环
  for (g in 1:length(Group)){

    G = Group[[g]]

    # cat(paste0("G",g,"\n"))
    if(is.null(G)){
      cat("G",g,"is null \n")
      Group[[g]] = NA
      next
    }else if(length(G)==0){
      cat("G",g,"is null \n")
      Group[[g]] = NA
      next
    }else if(nrow(G[[1]])==0){
      cat("G",g,"is null \n")
      Group[[g]] = NA
      next
    }else {
      cat("G",g,"is not null \n")
    }
    # if(nrow(as.data.frame(G))==0){
    #   cat("G",g,"is null \n")
    #   G = NULL
    #   next
    # }

    # 1. 在这一段中按照不同的簇分组，组内按照pseudotime连接
    cat("Loop in every layer \n")
    for(d in G){
      if(nrow(d)==0) next;
      n = d$c_names
      # cat(n, "\n")
      for(i in 1:nrow(d) ){
        #print(i)
        if(i==nrow(d)) break;
        adj_matrix[n[i], n[i+1]] = 1
        #adj_matrix[n[i+1], n[i]] = 1
      }
    }

    # 2. 为这一段中的每一个小组找到上家
    # g=1; g=2; g=...
    cat("Connect between layers \n")
    if(g==1){
      next;
    }else if(g==2){

      G2 = Group[[1]]

      # check if all tb < T
      # if(nrow(G2[[1]])==0){
      if(is.na(G2)|is.null(G2)){
        tmp = do.call(rbind,G)
        tmp = tmp[order(tmp$T),]
        min.name = tmp[1,]$c_names
        for(n in names(G)){
          x_c = G[[n]][1,]$"c_names"
          if(min.name!=x_c){
            adj_matrix[min.name, x_c] = 1
          }
        }
      }else{
        y_c = tail(G2[[1]],1)$c_names
        for(n in names(G)){
          x_c = G[[n]][1,]$"c_names"
          adj_matrix[y_c, x_c] = 1
        }
      }

    }else{

      G_n = str_split_fixed(names(G),"_",2)
      G_n[,2] = paste0("_", G_n[,2])

      for(n in G_n[,2]){
        cat("Group name: ", n, "\n")
        for(g2 in (g-1):1){
          cat("----Last Group: ", g2, "\n")
          G2 = Group[[g2]]
          aim = str_sub(n, 1, 2*(g2-1))
          ind = grep(aim, names(G2))
          if(length(ind)!=0){
            cat("----Find it!! \n")
            x_n = paste0(G_n[1,1],n)
            y_n = names(G2)[ind]
            x_c = G[[x_n]][1,]$"c_names"
            y_c = tail(G2[[y_n]],1)$c_names

            adj_matrix[y_c, x_c]=1
            break;
          }

          # 没找到上家，直接连接到第一个点
          if(g2==1){
            cat("----No find. \n ")
            # 生成原始的df
            # dd = lapply(tb_group_list, function(x) {
            #   if(!is.null(x)){
            #     do.call(rbind,x) %>% select(c_names, T) %>% data.frame
            #   } })
            # orig_df = do.call(rbind,dd)
            # start = orig_df[order(orig_df$T),][1,"c_names"]
            # x_n = paste0(G_n[1,1],n)
            # x_c = G[[x_n]][1,]$"c_names"
            # adj_matrix[start, x_c] = 1
          }

        }
      }
      ##
    }
  }
  return(adj_matrix)
}


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' GetBinTree
#'
#' @description
#' Get the deconstructed multiple trajectories,
#' and display them through a single-fork binary tree
#'
#' @param object MGPfact object
#' @param save Logical value, whether save the result in object
#'
#' @export
#'
GetBinTree <- function(object,
                       save = TRUE){

  sdf = GetMURPInfo(object)
  sdf = sdf[order(sdf$T, decreasing = FALSE), ]
  P = getParams(object, "murp_number")
  L = getParams(object, "trajectory_number")
  Q = getParams(object, "murp_pc_number")

  mamba = object@OptimResult$pse
  murp = object@MURP
  metadata = object@MetaData
  chains = mamba@chains

  ## 0. 获取bin_f和Tb
  bin_f <- sdf[, c("T", paste0("C0_",1:L))]
  Tb <- sdf[1,paste0("Tb_",1:L)] %>% as.matrix %>% as.vector
  names(Tb) <- paste0("Tb_", 1:L)

  ## 1. Segment bin_f according to different differentiation situations (tb)
  ## length(df_list) = length(Tb)
  tb_cut_list <- lapply(1:L, function(i){
    tb = paste0("Tb_", i)
    # cat("tb: ", tb, "\n")
    pos = split_tb(bin_f$T, Tb[tb])
    lapply(pos, function(x) bin_f[x,])
  })
  names(tb_cut_list) <- paste0("Tb_", 1:L)

  ## 2. group the (bin_f) of each trajectory according to the fork (T>tb)
  ## length(group_list) = length(Tb)
  tb_group_list <- lapply(1:L, function(i){
    tb = paste0("Tb_", i)
    l = paste0("C0_",i)
    # cat("tb: ", tb, " C: ", l, "\n")
    split_binary_tree(tb_cut_list[[tb]], sep = l)
  })

  ## 3. get the adjacency matrix and get the binary tree
  binary_tree_list <- lapply(tb_group_list, function(tb_group){
    # print(paste0("tree:",i))
    graph_from_adjacency_matrix(adj_tb_tree(tb_group, P, point_names =  rownames(sdf)) )
  })

  ## 4. get layout of each tree
  layout_list <- lapply(1:L, function(i){
    layout_as_tree(binary_tree_list[[i]],
                   root = rownames(bin_f)[1],
                   circular = FALSE,
                   mode = 'out',
                   flip.y = TRUE)
  })

  ## 5. get each bin tree
  centers <- murp$Recommended_K_cl$centers
  rownames(centers) <- paste0("T[",1:nrow(centers),"]")
  centers <- centers[rownames(bin_f),]
  bintree_all <- lapply(1:L, function(i){
    tree <- binary_tree_list[[i]]
    layout <- layout_list[[i]]
    edge <- as.data.frame(get.edgelist(tree), stringsAsFactors = FALSE)
    vertex <- data.frame( name = rownames(bin_f),
                          label = str_sub(rownames(bin_f), start = 3, end = -2),
                          c_label = bin_f[ ,paste0("C0_",i)],
                          pse = bin_f[ ,"T"],
                          centers,
                          sdf)
    newGraph <- graph_from_data_frame(edge, directed = TRUE, vertices = vertex)
    list(edge = edge, vertex = vertex, layout = layout, graph = newGraph)
  })
  names(bintree_all) <- paste0("trajectory",1:L)
  if(save){
    save(bintree_all, file = "2_pseudotime/bintree_all.rda")
  }
  object@Tree$bintree = bintree_all

  return(object)
}


#' GetTbTree
#'
#' @description
#' Get a multi-forked binary-tree by merging different trajectories
#'
#' @param object MGPfact object
#' @param save Logical value, whether save the result in object
#'
#' @export
#'
GetTbTreeDrq <- function(object,
                         save = TRUE){

  ### prepare
  sdf = GetMURPInfo(object)
  sdf = sdf[order(sdf$T, decreasing = FALSE), ]
  P = getParams(object, "murp_number")
  L = getParams(object, "trajectory_number")
  Q = getParams(object, "murp_pc_number")
  mamba = object@OptimResult$pse
  murp = object@MURP
  metadata = object@MetaData
  chains = mamba@chains

  ### function
  make_adj <- function(names){
    tmp = matrix(0, length(names), length(names)-1)
    diag(tmp) = 1
    tmp = cbind(matrix(0,length(names),1), tmp)
    dimnames(tmp) = list(names, names)
    return(tmp)
  }

  ###########################################
  ### 0. 获取tb_f和Tb
  tb_f <- sdf[, c("T", paste0("C0_",1:L))]
  Tb <- sdf[1,paste0("Tb_",1:L)] %>% as.matrix %>% as.vector
  names(Tb) <- colnames(sdf)[(2+L*2):(1+3*L)]

  ### 1. make a null matrix
  adj = matrix(0, P, P)
  dimnames(adj) = list(rownames(sdf), rownames(sdf))

  ## 2. 寻找跟节点在哪个TB处
  # for(i in 1:L){ # 选择第一个大于>min(T)的tb，相当于前面的都要剔除
  #   cat("root:", i, "\n")
  #   tb = Tb[i]
  #   if(tb<min(tb_f$T)){
  #     next
  #   } else{
  #     names = rownames(tb_f)[which(tb_f$T<=tb)] # root 节点前面的连接
  #     adj[names, names] = make_adj(names)
  #     end = tail(names, 1) ### end point
  #     root = i
  #     break
  #   }
  # }

  index = which(Tb > min(tb_f$T)) # 所有Tb都要保留
  if(length(index)==0){ # 没有tb > min(T), 一条直线
    adj <- matrix(0, nrow =  P, ncol =  P)
    rownames(adj) <- paste0("T[",1:P,"]")
    colnames(adj) <- paste0("T[",1:P,"]")
    for(i in 1:(nrow(tb_f)-1) ){
      adj[rownames(tb_f)[i],rownames(tb_f)[i+1]] = 1
    }
  }else if(length(index)==1){ # 只有1个tb > min(T), 一个二叉树
    pos = split_tb(tb_f$T, Tb[index])
    tb_cut = lapply(pos, function(x) tb_f[x,])
    # tb_group = split_binary_tree(tb_cut, sep = paste0("L",index))
    tb_group = split_binary_tree(tb_cut, sep = paste0("C0", str_sub(names(index), start= 3) ))
    adj = adj_tb_tree(tb_group, P, point_names =  rownames(sdf))
  }else{ # 两个及以上，无论如何以第一个tb作为第一个root
    tb = Tb[1]
    names = rownames(tb_f)[which(tb_f$T<=tb)] # root 节点前面的连接
    if(length(names)==0){
      names = rownames(tb_f)[1]
      first_t =  tb_f$T[1]
      end = rownames(tb_f)[1]
    }else{
      adj[names, names] = make_adj(names)
      end = tail(names, 1) ### end point
      first_t =  Tb[1]
    }
    # if(length(names)!=0){
    #   adj[names, names] = make_adj(names)
    #   end = tail(names, 1) ### end point
    #   first_t =  Tb[1]
    # }else{
    #   end = rownames(tb_f)[1]
    #   first_t =  tb_f$T[1]
    # }
    root = 1

    ## 3. 在root之后的分成c1, c2
    root_c = paste0("C0_",root)
    root_ind1 = rownames(tb_f)[which(tb_f$T > first_t & tb_f[,root_c]==1)]
    root_ind2 = rownames(tb_f)[which(tb_f$T > first_t & tb_f[,root_c]==2)]
    ind1_list = list( list(ind = root_ind1, end = end) ) # 把应该连接的上一个根节点也存起来
    ind2_list = list( list(ind = root_ind2, end = end) )
    adj[root_ind1, root_ind1] = make_adj(root_ind1) ## 第一个二叉树建立好
    adj[root_ind2, root_ind2] = make_adj(root_ind2)
    adj[end, root_ind1[1]] = 1
    adj[end, root_ind2[1]] = 1

    for(i in (root+1):L){

      # cat("adj", i, "\n")
      tb = Tb[i]
      c = paste0("C0_",i)
      new_ind1_list = list() # 创建空list准备放入
      new_ind2_list = list()

      # 上一个L中c1的后面
      for(j in 1:length(ind1_list)){

        # cat(j, "\n")
        ind = ind1_list[[j]]$ind
        if(is.null(ind)) next
        tmp = tb_f[ind,]
        median = rownames(tmp)[which(tmp$T<=tb)] # 中间的
        if(length(median)==0){
          new_end = end
        }else{
          new_end = tail(median,1) # 更新本次tb的end
        }

        upper = rownames(tmp)[which(tmp$T > tb)]
        adj[upper, ] = 0
        adj[ ,upper] = 0
        # adj[new_end, ] = 0

        # if(length(new_end)==0) new_end = end
        tmp2 = tmp[upper,,drop=FALSE]
        if(nrow(tmp2)==0) next
        new_id1 = rownames(tmp2)[which(tmp2[, c]==1)]
        new_id2 = rownames(tmp2)[which(tmp2[, c]==2)]
        # adj[new_end, new_id1[1]] = 1
        # adj[new_end, new_id2[1]] = 1

        # 把对应的新生成的c1,c2放在新的list里面
        if(length(new_id1)!=0){
          adj[new_end, new_id1[1]] = 1
          adj[new_id1,new_id1] = make_adj(new_id1)
          new_ind1_list = append(new_ind1_list,  list(list(ind = new_id1, end = new_end)) )
        }
        if(length(new_id2)!=0){
          adj[new_end, new_id2[1]] = 1
          adj[new_id2,new_id2] = make_adj(new_id2)
          new_ind2_list = append(new_ind2_list, list(list(ind = new_id2, end = new_end)) )
        }
      }

      # 上一个L中c2的后面
      for(j in 1:length(ind2_list)){
        ind = ind2_list[[j]]$ind
        if(is.null(ind)) next
        tmp = tb_f[ind,]
        median = rownames(tmp)[which(tmp$T<=tb)] # 中间的
        if(length(median)==0){
          new_end = end
        }else{
          new_end = tail(median,1) # 更新本次tb的end
        }

        upper = rownames(tmp)[which(tmp$T > tb)]
        adj[upper, ] = 0; adj[ ,upper] = 0;
        # adj[new_end, ] = 0

        # if(length(new_end)==0) new_end = end
        tmp2 = tmp[upper,]
        if(nrow(tmp2)==0) next
        new_id1 = rownames(tmp2)[which(tmp2[, c]==1)]
        new_id2 = rownames(tmp2)[which(tmp2[, c]==2)]
        # adj[new_end, new_id1[1]] = 1
        # adj[new_end, new_id2[1]] = 1

        # 把对应的新生成的c1,c2放在新的list里面
        if(length(new_id1)!=0){
          adj[new_end, new_id1[1]] = 1
          adj[new_id1,new_id1] = make_adj(new_id1)
          new_ind1_list = append(new_ind1_list,  list(list(ind = new_id1, end = new_end)) )
        }
        if(length(new_id2)!=0){
          adj[new_end, new_id2[1]] = 1
          adj[new_id2,new_id2] = make_adj(new_id2)
          new_ind2_list = append(new_ind2_list, list(list(ind = new_id2, end = new_end)) )
        }
      }
      ind1_list = new_ind1_list
      ind2_list = new_ind2_list
    }
  }

  y = adj
  # x = c(names, median, new_id1, new_id2)
  # y = adj[x,x]
  tmp <- graph_from_adjacency_matrix(y)
  edge <- as.data.frame(get.edgelist(tmp), stringsAsFactors = FALSE)
  vertex = data.frame(row.names = rownames(y),
                      name = rownames(y))
  g <- graph_from_data_frame(edge, directed = TRUE, vertices = vertex)
  layout <- layout_as_tree(g,
                           root = rownames(y)[1],
                           circular = FALSE,
                           mode = 'out',
                           flip.y = TRUE)
  ggraph(g, layout = layout) +
    geom_edge_diagonal(alpha = 0.7, width = 0.5, check_overlap = FALSE) +
    geom_node_point(color = "blue",size = 6, alpha = 0.6) +
    geom_node_text(aes(label = .data$name), size = 4) +
    coord_flip() +
    scale_y_reverse() +
    labs(title = paste0("L",i), color = "Bif") +
    rj.graph.ftheme

  tb_adj_matrix = adj

  ### 2. 画图准备, get tbtree list
  tmp <- graph_from_adjacency_matrix(tb_adj_matrix)
  edge <- as.data.frame(get.edgelist(tmp), stringsAsFactors = FALSE)
  rat <- table(murp$Recommended_K_cl$cluster)/length(murp$Recommended_K_cl$cluster)
  names(rat) <- paste0("T[", names(rat), "]")

  centers <- murp$Recommended_K_cl$centers
  rownames(centers) <- paste0("T[",1:nrow(centers),"]")
  centers <- centers[rownames(tb_f),]

  vertex = data.frame(row.names = rownames(tb_f),
                      name = rownames(tb_f),
                      label = str_sub(rownames(tb_f), start = 3, end = -2),
                      pse = tb_f$T,
                      grp = rownames(tb_f),
                      rat = as.vector(rat[rownames(tb_f)]),
                      sdf,
                      centers)

  ### 做一个C_all
  # vertex = tbtree$vertex
  # edge = tbtree$edge
  # layout = tbtree$layout
  tbx = unlist(vertex[1,paste0("Tb_",1:L),drop=TRUE])
  tbx = c(0,tbx,1)
  tmp = cut(vertex$T,tbx,labels = FALSE)
  tmp[which(tmp==1)] = "start"
  for(i in 2:length(tbx)){
    ind = which(tmp==i)
    l = i-1
    tmp[ind] = paste0("L",l,"_",vertex[ind,paste0("C_",l)])
  }
  vertex$C_all = tmp

  for(i in 1:L){ vertex[,paste0("C0_",i)] = as.factor(vertex[,paste0("C0_",i)]) }
  graph_tbtree <- graph_from_data_frame(edge, directed = TRUE, vertices = vertex)
  layout <- layout_as_tree(graph_tbtree,
                           root = rownames(tb_f)[1],
                           circular = FALSE,
                           mode = 'out',
                           flip.y = TRUE)
  tbtree <- list(edge = edge,
                 vertex = vertex,
                 graph = graph_tbtree,
                 layout = layout,
                 adj = tb_adj_matrix)

  if(save){
    save(tbtree, file = "2_pseudotime/tbtree.rda")
  }
  object@Tree$tbtree = tbtree

  return(object)
}

#' GetTbTree
#'
#' @description
#' get tbtree result of murp
#'
#' @param object MGPfact object
#' @param save Logical value, whether save tbtree result in mgpfact
#'
#' @export
#'
GetTbTree <- function(object,
                      save = TRUE){

  ### judge sdf
  sdf = GetMURPInfo(object)
  sdf = sdf[order(sdf$T, decreasing = FALSE), ]
  P = getParams(object, "murp_number")
  L = getParams(object, "trajectory_number")
  Q = getParams(object, "murp_pc_number")

  mamba = object@OptimResult$pse
  murp = object@MURP
  metadata = object@MetaData
  chains = mamba@chains

  ### 0. 获取tb_f和Tb
  tb_f <- sdf[, c("T", paste0("C0_",1:L))]
  Tb <- sdf[1,paste0("Tb_",1:L)] %>% as.matrix %>% as.vector
  names(Tb) <- paste0("Tb_",1:L)
  # names(Tb) <- colnames(sdf)[(2+L*2):(1+3*L)]


  ### 1. 获得Tb-tree的邻接矩阵tb_adj_matrix
  index = which(Tb > min(tb_f$T))
  if(length(index)==0){
    itb = length(Tb)
    pos = split_tb(tb_f$T, Tb[itb])
    tb_cut = lapply(pos, function(x) tb_f[x,])
    tb_group = split_binary_tree(tb_cut, sep = paste0("C0_", itb))
    adj_matrix = adj_tb_tree(tb_group, P, point_names =  rownames(sdf))

  }else if(length(index)==1){
    pos = split_tb(tb_f$T, Tb[index])
    tb_cut = lapply(pos, function(x) tb_f[x,])
    # tb_group = split_binary_tree(tb_cut, sep = paste0("L",index))
    tb_group = split_binary_tree(tb_cut, sep = paste0("C0", str_sub(names(index), start= 3) ))
    adj_matrix = adj_tb_tree(tb_group, P, point_names =  rownames(sdf))
  }else{
    pos <- split_tb(tb_f$T, Tb)
    tb_cut_list <- lapply(pos, function(x) tb_f[x,])
    tb_group_list <- split_tb_tree(tb_cut_list)
    adj_matrix <- adj_tb_tree(tb_group_list, P, point_names =  rownames(sdf))
  }
  tb_adj_matrix = adj_matrix

  ### 2. 画图准备, get tbtree list
  tmp <- graph_from_adjacency_matrix(tb_adj_matrix)
  edge <- as.data.frame(get.edgelist(tmp), stringsAsFactors = FALSE)
  rat <- table(murp$Recommended_K_cl$cluster)/length(murp$Recommended_K_cl$cluster)
  names(rat) <- paste0("T[", names(rat), "]")

  centers <- murp$Recommended_K_cl$centers
  rownames(centers) <- paste0("T[",1:nrow(centers),"]")
  centers <- centers[rownames(tb_f),]

  vertex = data.frame(row.names = rownames(tb_f),
                      name = rownames(tb_f),
                      label = str_sub(rownames(tb_f), start = 3, end = -2),
                      pse = tb_f$T,
                      grp = rownames(tb_f),
                      rat = as.vector(rat[rownames(tb_f)]),
                      sdf,
                      centers)

  ### 做一个C_all
  # vertex = tbtree$vertex
  # edge = tbtree$edge
  # layout = tbtree$layout
  tbx = unlist(vertex[1,paste0("Tb_",1:L),drop=TRUE])

  if(1 %in% tbx){ tbx = c(0,tbx) }else{ tbx = c(0,tbx,1) }
  tmp = cut(vertex$T,tbx,labels = FALSE)
  tmp[which(tmp==1)] = "start"
  for(i in 2:length(tbx)){
    ind = which(tmp==i)
    l = i-1
    tmp[ind] = paste0("L",l,"_",vertex[ind,paste0("C_",l)])
  }
  vertex$C_all = tmp

  for(i in 1:L){ vertex[,paste0("C0_",i)] = as.factor(vertex[,paste0("C0_",i)]) }
  graph_tbtree <- graph_from_data_frame(edge, directed = TRUE, vertices = vertex)
  layout <- layout_as_tree(graph_tbtree,
                           root = rownames(tb_f)[1],
                           circular = FALSE,
                           mode = 'out',
                           flip.y = TRUE)
  tbtree <- list(edge = edge,
                 vertex = vertex,
                 graph = graph_tbtree,
                 layout = layout,
                 adj = tb_adj_matrix)

  if(save){
    save(tbtree, file = "2_pseudotime/tbtree.rda")
  }
  object@Tree$tbtree = tbtree
  return(object)
}

#' GetTbTreeAllpoint
#'
#' @description
#' add all point in tbtree
#'
#' @param object MGPfact object
#' @param labels some attribute about cells
#' @param save Logical value, whether save tbtree result in mgpfact
#'
#' @export
#'
GetTbTreeAllpoint <- function(object,
                              save = TRUE,
                              labels = NULL){

  ### prepare
  sdf = GetMURPInfo(object)
  sdf = sdf[order(sdf$T, decreasing = FALSE), ]
  tbtree = GetTbTreeResult(object)
  P = getParams(object, "murp_number")
  L = getParams(object, "trajectory_number")
  Q = getParams(object, "murp_pc_number")

  mamba = object@OptimResult$pse
  murp = object@MURP
  metadata = object@MetaData
  chains = mamba@chains

  tb_f <- sdf[, c("T", paste0("C0_",1:L))]

  ### 1. 分别把每个簇的矩阵提取出来做成list
  sdata_cluster <- tapply(1:length(murp$Recommended_K_cl$cluster),
                          murp$Recommended_K_cl$cluster,
                          function(x,y){
                            exp = matrix(murp$rawdata[x,], ncol = ncol(murp$rawdata))
                            dimnames(exp) = list(rownames(murp$rawdata)[x], colnames(murp$rawdata))
                            exp })

  ### 2. 把每个中心点加入分别加入每个簇并命名
  sdata_cluster_centers = lapply(1:P, function(i){
    # print(i)
    exp = sdata_cluster[[i]]
    centers = matrix(murp$Recommended_K_cl$centers[i,], ncol = ncol(murp$rawdata))
    tmp = rbind(exp, centers)
    # rownames(tmp) = c(rownames(exp), rownames(murp$Recommended_K_cl$centers)[i])
    ci = rownames(murp$Recommended_K_cl$centers)[i]
    rownames(tmp) = c(rownames(exp),
                      paste0("T[", ci, "]"))
    tmp
  })

  ### 3. 分别把每个簇的距离矩阵计算出来
  cores = 3
  dist_list <- lapply(1:P, function(i){
    exp = sdata_cluster_centers[[i]]
    dist_m = parDist(x = exp, method = "euclidean", diag = TRUE, upper = TRUE, threads = cores-2)
    dist_m
  })

  ### 4. 针对每个簇建立最小生成树
  max_v = max(unlist(dist_list))
  min_v = 1e-9
  mst_list <- lapply(1:P, function(i){
    # print(i)
    dist_m = dist_list[[i]]
    if(max(dist_m)==0){
      dist_m = as.matrix(replace(dist_m, dist_m == 0, min_v))
    }else{
      dist_m = as.matrix(dist_m)
    }
    graph = graph.adjacency(dist_m, weighted=TRUE)
    mst = minimum.spanning.tree(graph)
    mst
  })

  ### 5. 获取所有簇的邻接矩阵
  adj_matrix_list <- lapply(1:P, function(i){
    # print(i)
    mst = mst_list[[i]]
    adj = as.matrix(as_adjacency_matrix(mst))
    adj
  })

  ### 6. 得到所有点（包括中心点）的邻接矩阵
  adj_final <- as.matrix(do.call(bdiag, adj_matrix_list))
  rownames(adj_final) <- unlist(lapply(adj_matrix_list, rownames))
  colnames(adj_final) <- rownames(adj_final)

  ### 7. 把tb-tree的结果加上
  cc = rownames(tbtree$adj)
  for(i in 1:P){
    for(j in 1:P){
      adj_final[cc[i],cc[j]] = tbtree$adj[cc[i],cc[j]]
    }
  }

  ### 8. 准备数据
  ### get edge
  tmp <- graph_from_adjacency_matrix(adj_final, mode = "directed")
  edge <- as.data.frame(get.edgelist(tmp), stringsAsFactors = FALSE)

  ### 每个簇的细胞个数/总细胞数
  cell_cluster_centers <- paste0("c", murp$Recommended_K_cl$cluster)
  rat <- table(murp$Recommended_K_cl$cluster)/length(murp$Recommended_K_cl$cluster)
  names(rat) <- paste0("T[", names(rat), "]")

  ### 点的gene属性
  centers <- murp$Recommended_K_cl$centers
  rownames(centers) <- paste0("T[",1:nrow(centers),"]")
  centers <- centers[rownames(tb_f),]
  sdata2 <- rbind(centers, murp$rawdata)

  ## 得到vertex
  vertex <- data.frame(row.names = c(rownames(tb_f), rownames(murp$rawdata)),
                       name = c(rownames(tb_f), rownames(murp$rawdata)),
                       pse = c(tb_f$T, tb_f[cell_cluster_centers,"T"]),
                       pse_na = c(tb_f$T, rep("",nrow(murp$rawdata))),
                       grp = c(rownames(tb_f), cell_cluster_centers),
                       grp_na = c(rownames(tb_f), rep("",nrow(murp$rawdata))),
                       rat = c(rat[rownames(tb_f)], rep(0,nrow(murp$rawdata))),
                       alpha = c(rep(0.8,nrow(tb_f)), rep(0.5,nrow(murp$rawdata))),
                       alpha_murp = c(rep(1,nrow(tb_f)), rep(0,nrow(murp$rawdata))),
                       size = c(tb_f$T, rep(0.01,nrow(murp$rawdata))),
                       size2 = c(tb_f$T, rep(0,nrow(murp$rawdata))),
                       sdata2)

  ### 计算每种细胞类型在每个MURP中所占的比例（中心点的比例设为0）
  if(!is.null(labels)){

    df2_list <- lapply(labels, function(lab){
      cat(lab, "\n")
      tmp <- table(object@MetaData[,c("murp_cluster",lab)]) %>% as.data.frame
      df <- dcast(tmp, murp_cluster~get(lab))
      df <- apply(df, 2, as.numeric)
      rownames(df) <- df[,"murp_cluster"]
      df <- df[,-1,drop=FALSE]

      if(ncol(df)==1){
        df[,1] = rep(1, nrow(df))
      }else{
        df <- apply(df,1,function(x){ x/sum(x) })
        df <- t(df)
      }

      rownames(df) <- paste0("T[",rownames(df),"]")
      df <- df[rownames(tb_f),,drop=FALSE]

      tmp <- matrix(0, nrow =  nrow(object@assay$data_matrix), ncol =  ncol(df))
      rownames(tmp) <- rownames(murp$rawdata)
      df2 <- rbind(df, tmp)
    })
    df2 <- do.call(cbind, df2_list) %>% data.frame
    if(nrow(df2)==0) df2 = NULL

    ### 点的celltype属性
    if(is.null(labels)){
      meta_tmp = NULL
    }else if(length(labels)==1){
      meta_tmp <- c(sdf[rownames(tb_f), labels],
                    metadata[rownames(murp$rawdata),labels]) %>% data.frame
      colnames(meta_tmp) <- labels
    }else{
      meta_tmp <- rbind(sdf[rownames(tb_f), labels],
                        metadata[rownames(murp$rawdata),labels])
    }

    ### 合并
    vertex <- data.frame(vertex,meta_tmp,df2)
  }

  ## graph_tbtree_all
  graph_tbtree_all <- graph_from_data_frame(edge, directed = FALSE, vertices = vertex)

  ### 9. save
  s_ind <- unique(intersect(grep("T\\[",edge[,1]),grep("T\\[",edge[,2])))
  b_ind <- setdiff(1:nrow(edge), s_ind)
  edge$weight <- 10
  edge$weight[b_ind] <- 0.1
  edge$weight[s_ind] <- 1
  new_edge <- edge
  new_graph <- graph_from_data_frame(new_edge, directed = FALSE, vertices = vertex)
  start = rownames(sdf)[which.min(sdf$T)]
  fork_p = names(which(degree(tbtree$graph)>=3))

  ### 10. backbone
  bb <- layout_as_backbone(new_graph, keep = 0.4, backbone = TRUE)
  tbtree_all <- list(edge = new_edge,
                     vertex = vertex,
                     graph = new_graph,
                     adj = adj_final,
                     bb = bb,
                     start = start,
                     fork_p = fork_p)

  if(save){
    save(tbtree_all, file = "2_pseudotime/tbtree_all.rda")
  }
  object@Tree$tbtree_all = tbtree_all
  return(object)
}

#' GetTreeResult
#' @description
#' extract bifurcation result
#' @param object MGPfact object
#' @export
#'
GetBinTreeResult <- function(object){
  r = object@Tree[["bintree"]]
  return(r)
}

#' GetTbTreeResult
#' @description
#' extract consensus tree result
#' @param object MGPfact object
#' @export
#'
GetTbTreeResult <- function(object){
  r = object@Tree[["tbtree"]]
  return(r)
}

#' GetTbTreeAllResult
#' @description
#' extract consensus tree result, which including all cells
#' @param object MGPfact object
#' @export
#'
GetTbTreeAllResult <- function(object){
  r = object@Tree[["tbtree_all"]]
  return(r)
}



# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                          functions related to murp
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' om
#'
#' @Description:
#' omega compute
#'
#' @param x: expression matrix, row is cell, col is gene
#' @export
#'
om <- function(x){
  r = -0.036 + 78.34 * (1/nrow(x)) + 1.85 * 1e-5 * ncol(x)
  return(r)
}

##### MURPDownsampling
#' @param name description
MURPDownsampling <- function(object,
                             omega = 0.5,
                             iter = 10,
                             seed = 723,
                             fast = T,
                             max_murp = 500,
                             cores = 1,
                             pca.center = FALSE,
                             pca.scale = FALSE,
                             plot = T){

  require(MURP)
  require(ggplot2)
  require(patchwork)

  object@MURP <- MURP(Data = object@assay$data_matrix,
                      cores = cores,
                      omega = omega,
                      iter = iter,
                      seed = seed,
                      fast = fast,
                      max_murp = max_murp)
  object@MetaData$murp_cluster <- object@MURP$Recommended_K_cl$cluster
  # object@MURP$centersPCA <- prcomp(object@MURP$Recommended_K_cl$centers, center = pca.center, scale. = pca.scale)
  object <- MURPPCA(object, pca.center = pca.center, pca.scale = pca.scale)

  if(plot){
    ggsave("1_murp/murp_bic_k.pdf", MURPNestedGridPlot(object@MURP), width = 3.5, height = 3.5)

    pl = PCANestedGridPlot(pca_result = object@MURP$centersPCA, sd_cutoff = 1, max_pc = 100)
    p = wrap_plots(pl, ncol = 3, guides = "collect")
    ggsave("1_murp/pca_scree_30pc.pdf",p, width = 10, height = 3.7)
  }

  command = GetCommand()
  object@Command$murp$MURPDownsampling = command
  return(object)
}

##### MURPDownsampling
#' @param name description
MURPPCA <- function(object,
                    pca.center = FALSE,
                    pca.scale = FALSE){

  object@MURP$centersPCA <- prcomp(object@MURP$Recommended_K_cl$centers, center = pca.center, scale. = pca.scale)

  command = GetCommand()
  object@Command$murp$MURPPCA = command
  return(object)
}
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                              useful function
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' %nin%
#'
#' @export
'%nin%' <- function(x,y)!('%in%'(x,y))

#' rmx
#' @description
#' Remove other variables in the environment except x
#'
#' @param x character vector
#' @export
#'
rmx <- function(x){
  # rm(list=setdiff(ls(),x), envir = .GlobalEnv)
  all = ls(envir = .GlobalEnv)
  rm(list=setdiff(all,x), envir = .GlobalEnv)
}

#'
#' RenameCluster
#'
#' @description
#' A replacement tool, such as replacing all "a" in the vector vector with "1"
#'
#' @param vector a vector to replace
#' @param replace the content to be replaced
#'
#' @usage
#' RenameCluster(c(1,1,1,2), c("1"="a","2"="b"))
#' "a" "a" "a" "b"
#'
#' @export
RenameCluster <- function(vector = NULL,
                          replace = NULL){

  orig_n <- names(replace)
  new_n <- replace

  for(i in orig_n ){
    vector[which(vector==i)] = new_n[i]
  }
  return(vector)
}

#'
#' enum
#'
#' @Description:
#' tiny pipline
#' unique %>% length
#' @export
#'
enum <- function(a){
  length(unique(a))
}

#' jac
#'
#' @Description:
#' compute jaccard similarity between two vector
#'
#' @param a
#' @param b
#' @export
#'
jac <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (round(intersection/union,4))
}

#' aeq
#'
#' @Description:
#' all.equal
#'
#' @param a
#' @param b
#' @export
#'
eql <- function(x, y, ...) all.equal(as.vector(x), as.vector(y), ...)

#' subString
#' @description
#' split string with sep, get idx substring
#'
#' @param strings
#' @param idx
#' @param sep
#'
#' @export
#'
subString <- function(strings, idx, sep = NA){

  strings = as.character(strings)
  if(is.na(sep)){
    res = as.character(lapply(strings, function(x) paste(strsplit(x, "")[[1]][idx], collapse = "")))
  } else{
    res = sapply(strsplit(strings, sep), function(x) x[idx])
  }
  return(res)
}

#' ExceptSet
#'
#' @description
#'
#' @param set_list
#' @export
ExceptSet <- function(set_list = NULL){

  sets = set_list
  x = names(set_list)
  if(is.null(x)){
    x = 1:length(set_list)
  }

  com = combn(x, 2, simplify = TRUE)
  result = apply(com, 2, function(i){
    a = setdiff(sets[[i[1]]], sets[[i[2]]])
    b = setdiff(sets[[i[2]]], sets[[i[1]]])
    r = list(a, b)
    names(r) = c(i[1],i[2])
    r
  })

  names = apply(com, 2, function(i){
    paste0(i[1],"_",i[2])
  })

  names(result) = names

  return(result)
}

#' IntersectSet
#'
#' @description
#'
#' @param set_list
#' @export
IntersectSet <- function(set_list = NULL, num=3){

  sets = set_list
  x = names(set_list)
  if(is.null(x)){
    x = 1:length(set_list)
  }

  com = combn(x, num, simplify = TRUE)
  result = apply(com, 2, function(i){

    tmp = lapply(seq_along(i), function(j){
      sets[[i[j]]]
      })
    c = Reduce(intersect, tmp)
    if(length(c)==0) c = "none"
    return(c)
  })

  names = apply(com, 2, function(i){
    paste0(i, collapse="_")
  })

  names(result) = names

  return(result)
}

#' #' get_blank_coord
#' #'
#' #' @export
#' get_blank_coord <- function(p){
#'
#'   l_dat = layer_data(p,1)
#'   sapply(c("x","y"), function(coord_tag){
#'
#'     ## get range of coordinates x/y
#'     rge = ggplot_build(p)$layout$panel_params[[1]][[paste0(coord_tag,".range")]]
#'
#'     # 插值
#'     coord_df <- data.frame(coord = seq(rge[1],rge[2], length.out=1000),
#'                            tag = "0") %>%
#'       rbind(data.frame(coord = l_dat[,coord_tag], tag = "1")) %>%
#'       arrange(coord)
#'
#'     # 判断连续的空白区域（0是插入的空白点，1是原先的点）
#'     # 用start命名
#'     rl = rle(coord_df$tag)
#'     rl_i = c(1) ## named as start index
#'     for(i in 1:(length(rl$lengths)-1) ){
#'       rl_i = append(rl_i, rl_i[i]+rl$lengths[i])
#'     }
#'     rl_df = data.frame(row.names = rl_i,
#'                        ind = rl_i,
#'                        v = rl$values,
#'                        l = rl$lengths) %>%
#'       filter(v=="0")
#'     mm = which(rl_df$l==max(rl_df$l))
#'     if(length(mm)>1) mi = tail(mm,1) else mi = mm
#'     ind = rl_df[mi,"ind"]
#'
#'     # 利用取出的短的均值作为可以放label的坐标
#'     crd = coord_df[ind:(ind+rl_df[as.character(ind),"l"]-1), "coord"]
#'     if(coord_tag=="x"){
#'       if((max(coord_df$coord)-tail(crd,1))<=0.1){
#'         star = crd[1]-0.15
#'       }
#'       if((crd[1]-min(coord_df$coord))<=0.1){
#'         star = tail(crd,1)
#'       }
#'     }else{
#'       star = prod(crd)^(1/length(crd))
#'     }
#'
#'     return(star)
#'   } ) -> crdd
#'
#'   return(crdd)
#' }


#' get_blank_coord2
#'
#' @export
get_blank_coord2 <- function(p, x_cut = NULL, y_cut = NULL){

  l_dat = layer_data(p,1)
  mat = l_dat %>%
    dplyr::select(x,y) %>%
    mutate(value = 1) %>%
    reshape2::acast(x~y, value.var = "value", fill = 0,
                    drop = FALSE, fun.aggregate = max)

  rge = range(l_dat$x)
  coord = seq(rge[1],rge[2], length.out=1000)
  unq_nodes = union(coord, rownames(mat))
  dim = length(unq_nodes)
  mat <- rbind(
    mat,
    matrix(
      rep(0, ncol(mat) * (dim - nrow(mat))),
      ncol = ncol(mat),
      dimnames = list(setdiff(unq_nodes,rownames(mat)),
                      colnames(mat))
    )
  )

  rge = range(l_dat$y)
  coord = seq(rge[1],rge[2], length.out=1000)
  unq_nodes = union(coord, colnames(mat))
  dim = length(unq_nodes)
  mat <- cbind(
    mat,
    matrix(
      rep(0, nrow(mat) * (dim - ncol(mat))),
      nrow = nrow(mat),
      dimnames = list(rownames(mat),
                      setdiff(unq_nodes,colnames(mat)))
    )
  )
  mat = mat[order(rownames(mat)%>%as.numeric), order(colnames(mat)%>%as.numeric)]

  if(!is.null(x_cut)){
    xl = rownames(mat) %>% as.numeric
    xi = which((xl>x_cut[1]) & (xl < x_cut[2]))
    mat = mat[xi,]
  }
  if(!is.null(y_cut)){
    yl = colnames(mat) %>% as.numeric
    yi = which((yl > y_cut[1]) & (yl < y_cut[2]))
    mat = mat[,yi]
  }
  # 定义动态规划矩阵和最大子矩阵大小
  dp_mat <- matrix(0, nrow = nrow(mat), ncol = ncol(mat))
  max_size <- 0
  max_row <- 0
  max_col <- 0

  # 填充动态规划矩阵
  for (i in 1:nrow(mat)) {
    for (j in 1:ncol(mat)) {
      if (mat[i,j] == 0) {
        if (i == 1 || j == 1) {
          dp_mat[i,j] <- 1
        } else {
          dp_mat[i,j] <- min(dp_mat[i-1,j], dp_mat[i,j-1], dp_mat[i-1,j-1]) + 1
        }
        if (dp_mat[i,j] > max_size) {
          max_size <- dp_mat[i,j]
          max_row <- i
          max_col <- j
        }
      }
    }
  }

  # 打印最大子矩阵
  if (max_size == 0) {
    cat("No all-zero submatrix found.\n")
  } else {
    row_idx <- max_row - max_size + 1
    col_idx <- max_col - max_size + 1
    sub_mat <- mat[row_idx:(row_idx+max_size-1), col_idx:(col_idx+max_size-1)]
    # cat(sprintf("Largest all-zero submatrix found with size %dx%d at position (%d,%d):\n",
    #     max_size, max_size, row_idx, col_idx))
    # print(sub_mat)
  }

  y_rge = colnames(mat)[col_idx:(col_idx+max_size-1) ] %>% as.numeric %>% range
  x_rge = rownames(mat)[row_idx:(row_idx+max_size-1) ] %>% as.numeric %>% range

  # r = c(x=mean(x_rge), y = mean(y_rge))
  r = list(x=x_rge, y = y_rge)
  return(r)
}

#' RunTime
#'
#' @param t1
#' @param t2
#' @param timeunits
#'
#' @export
#'
RunTime <- function(t1, t2,
                    timeunits='secs_mins') {

  time_secs = as.numeric(difftime(t2,t1,units='secs'))
  time_mins = as.numeric(difftime(t2,t1,units='mins'))
  time_hours = as.numeric(difftime(t2,t1,units='hours'))

  myRunTime = data.frame(t1 = t1,
                         t2 = t2,
                         exetime_secs = time_secs,
                         exetime_mins = time_mins,
                         exetime_hours = time_hours)
  return(myRunTime)

}
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                                plot function
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' quantile_breaks
#'
#' @export
#'
quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}


#' geom_rectriangle
#'
#' @export
geom_rectriangle <- function(mapping = NULL, data = NULL,
                             stat = "identity", position = "identity",
                             ...,
                             linejoin = "mitre",
                             na.rm = FALSE,
                             show.legend = NA,
                             inherit.aes = TRUE) {

  require(vctrs)
  require(grid)
  require(ggplot2)

  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomRectriangle,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      linejoin = linejoin,
      na.rm = na.rm,
      ...
    )
  )
}

#' GeomRectriangle
#' @importFrom ggproto ggplot2
#' @export
#'
require(ggplot2)
GeomRectriangle <- ggproto(
  "GeomRectriangle", Geom,
  default_aes = aes(r = 1, colour = "grey35", fill = NA, size = 0.25, linetype = 1,
                    alpha = NA,type = "upper"),
  required_aes = c("x", "y"),
  draw_panel = function(self, data, panel_params, coord, linejoin = "mitre",type = "upper") {
    aesthetics <- setdiff(names(data), c("x", "y"))

    polys <- lapply(split(data, seq_len(nrow(data))), function(row) {
      rectriangle <- point_to_rectriangle(row$x, row$y, row$r, row$type)
      aes <- new_data_frame(row[aesthetics])[rep(1, 4), ]
      GeomPolygon$draw_panel(cbind(rectriangle, aes), panel_params, coord)
    })

    ggplot2:::ggname("geom_rectriangle", do.call("grobTree", polys))
  },
  draw_key = draw_key_polygon
)

#' point_to_rectriangle
#'
#' @export
point_to_rectriangle <- function(x, y, r, type = type) {
  r <- 0.5 * sign(r) * sqrt(abs(r))
  #r0 = 0.5
  xmin <- - r + x
  xmax <- r + x
  ymin <- - r + y
  ymax <- r + y
  if(type == "upper"){
    df = new_data_frame(list(
      y = c(ymax, ymax, ymin, ymax),
      x = c(xmin, xmax, xmin, xmin)
    ))
  }else if(type == "lower"){
    df = new_data_frame(list(
      y = c(ymax, ymin, ymin, ymax),
      x = c(xmax, xmax, xmin, xmax)
    ))
  }
  df
}

#' GeomSplitViolin
#'
#' @export
GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin,
                           draw_group = function(self, data, ..., draw_quantiles = NULL) {
                             data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
                             grp <- data[1, "group"]
                             newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
                             newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
                             newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])

                             if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                               stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
                                                                         1))
                               quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
                               aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
                               aesthetics$alpha <- rep(1, nrow(quantiles))
                               both <- cbind(quantiles, aesthetics)
                               quantile_grob <- GeomPath$draw_panel(both, ...)
                               ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
                             }
                             else {
                               ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
                             }
                           })

#' geom_split_violin
#'
#' @export
geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ...,
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE,
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin,
        position = position, show.legend = show.legend, inherit.aes = inherit.aes,
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                                Nebulosa
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' get_dens <- function(data, dens, method) {
#'   require(ks)
#'   if (method == "ks") {
#'     ix <- findInterval(data[, 1], dens$eval.points[[1]])
#'     iy <- findInterval(data[, 2], dens$eval.points[[2]])
#'     ii <- cbind(ix, iy)
#'     z <- dens$estimate[ii]
#'   } else if (method == "wkde") {
#'     ix <- findInterval(data[, 1], dens$x)
#'     iy <- findInterval(data[, 2], dens$y)
#'     ii <- cbind(ix, iy)
#'     z <- dens$z[ii]
#'   }
#'   z
#' }
#'
#' #' @title Estimate weighted kernel density
#' #' @author Jose Alquicira-Hernandez
#' #' @param w Vector with weights for each observation
#' #' @param x Matrix with dimensions where to calculate the density from. Only
#' #' the first two dimensions will be used
#' #' @param method Kernel density estimation method:
#' #' \itemize{
#' #' \item \code{ks}: Computes density using the \code{kde} function from the
#' #'  \code{ks} package.
#' #' \item \code{wkde}: Computes density using a modified version of the
#' #'  \code{kde2d} function from the \code{MASS}
#' #' package to allow weights. Bandwidth selection from the \code{ks} package
#' #'  is used instead.
#' #' }
#' #' @param adjust Numeric value to adjust to bandwidth. Default: 1. Not available
#' #'  for \code{ks} method
#' #' @param map Whether to map densities to individual observations
#' #' @return If \code{map} is \code{TRUE}, a vector with corresponding densities
#' #'  for each observation is returned. Otherwise,
#' #' a list with the density estimates from the selected method is returned.
#' #' @importFrom ks kde hpi
#' #' @examples
#' #'
#' #' dens <- Nebulosa:::calculate_density(iris[, 3], iris[, 1:2], method = "wkde")
#' calculate_density <- function(w, x, method, adjust = 1, map = TRUE) {
#'   if (method == "ks") {
#'     dens <- kde(x[, c(1, 2)],
#'                 w = w / sum(w) * length(w)
#'     )
#'   } else if (method == "wkde") {
#'     dens <- wkde2d(
#'       x = x[, 1],
#'       y = x[, 2],
#'       w = w / sum(w) * length(w),
#'       adjust = adjust
#'     )
#'   }
#'
#'   if (map) {
#'     get_dens(x, dens, method)
#'   } else {
#'     dens
#'   }
#' }



#' # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' #' GetMapClusterLabel
#' #'
#' #' @Description:
#' #' Map MURP clustering results to all cells
#' #'
#' #' @param mic mic results
#' #' @param cluster cluster object of kmeans
#' #'
#' GetMapClusterLabel <- function(orig_cluster = NULL,
#'                                second_cluster = NULL){
#'   tmp1 <- orig_cluster
#'   tmp2 <- second_cluster
#'
#'   for(i in 1:length(tmp1)){
#'     for(j in 1:length(tmp2)){
#'       # cat("i:",i,", j:", j, "\n")
#'       if(as.character(tmp1[i])==names(tmp2)[j]){
#'         # cat("yes \n")
#'         tmp1[i]=tmp2[j]
#'       }
#'     }
#'   }
#'   return(tmp1)
#' }

#' # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' #' Get z-scored data (urd)
#' #'
#' #' Returns z-scored data, calculated on the log-transformed, normalized data. Eliminates
#' #' any genes that produce NA values, which can occur if they have standard deviation 0.
#' #' @param data row is genes
#' #' @param genes (Character vector) Genes to return in z-scored matrix (default, all genes)
#' #'
#' #' @return Matrix of z-scored data (Warning: memory hog as data will no longer be sparse.)
#' #'
#' #' @export
#' GetZcoreData <- function(data, genes=NULL) {
#'   # Z-scored data
#'   if (is.null(genes)) genes <- rownames(data)
#'   data.mean <- apply(data[genes,], 1, mean)
#'   data.sd <- apply(data[genes,], 1, sd)
#'   z.data <- as.matrix(sweep(sweep(data[genes,], 1, data.mean, "-"), 1, data.sd, "/"))
#'   z.data <- z.data[complete.cases(z.data),]
#'   return(z.data)
#' }

#' # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' #' Normalize
#' #'
#' #' @param count row is gene, col is cell
#' #' @param base default is 10, x=log(y,10)
#' #' @return a matrix with Gene Number * Cell Number
#' #'
#' #' @examples
#' #'
#' #' data(cts)
#' #' normed <- Normalize(count = cts, base = 10)
#' #'
#' Normalize <- function(count = NULL, base = 10){
#'   bc_sums <- colSums(count)
#'   median_sum <- median(bc_sums)
#'   new_matrix <- sweep(count, 2, median_sum/bc_sums, "*")
#'   data <- log(1 + new_matrix, base = base)
#'
#'   return(data)
#' }



#' # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' #' GetFinalData
#' #'
#' #' @param selectG a vector include genes'name that you choosed
#' #' @param geneMeta compute from mkGeneMeta
#' #' @param normed compute from Normalize
#' #' @param zScore locical value, make the expression of each gene a standard normal distribution,
#' #' x-mean/sd, default FALSE
#' #' @return a matrix with Cell Number * Gene Number
#' #' @examples
#' #'
#' #'  data(cts)
#' #'  data(selectG)
#' #'  normed <- Normalize(count = cts, base = 10)
#' #'  geneMeta <- mkGeneMeta(exprs = normed)
#' #'  # selectG is character vector
#' #'  sdata <- GetFinalData(selectG = selectG, geneMeta = geneMeta, normed = normed)
#' #'
#' GetFinalData <- function(selectG = NULL,
#'                          geneMeta = geneMeta,
#'                          normed = normed,
#'                          zScore = FALSE){
#'
#'   geneMeta$Select <- 'N'
#'   geneMeta[selectG,'Select'] <- 'Y'
#'
#'   sel_data <- t(normed[rownames(geneMeta)[which(geneMeta$Select == 'Y')],])
#'   sel_data <- as.matrix(sel_data)
#'
#'   if(zScore){
#'     final.gene.mean <- apply(sel_data, 2, mean)
#'     final.gene.sd <- apply(sel_data, 2, sd)
#'     sdata <- apply(sel_data, 1, function(x){(x-final.gene.mean)/final.gene.sd})
#'   }else{
#'     sdata <- as.matrix(sel_data)
#'   }
#'
#'   # check max(expr) of every gene
#'   max <- apply(sdata,2,max)
#'   sdata <- sdata[, which(max!=0)]
#'
#'   return(sdata)
#'
#' }


#' # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' #' mkGeneMeta
#' #'
#' #' @param exprs row is gene, col is cell
#' #' @param presentCall.cutoff proportion of cells greater than the expression threshold,
#' #' this parameter is the expression threshold you set, which can be a value or a vector
#' #' default: seq(0,1,0.5)
#' #'
#' #' @importFrom pbapply pbsapply
#' #'
#' #' @return attribute matrix of all genes
#' #' @examples
#' #'
#' #'  data(cts)
#' #'  normed <- Normalize(count = cts, base = 10)
#' #'  geneMeta <- mkGeneMeta(exprs = normed)
#' #'
#' mkGeneMeta <- function(exprs = NULL,
#'                        presentCall.cutoff = seq(0.,1,.5)){
#'
#'   require(pbapply)
#'   require(matrixStats)
#'
#'   cat('Calculating Mean\n')
#'   # Mean <- apply(exprs, 1, mean)
#'   Mean <- rowMeans2(as.matrix(exprs))
#'   cat('Calculating Sd and Var\n')
#'   # Sd <- apply(exprs, 1, sd)
#'   Sd <- rowSds(as.matrix(exprs))
#'   Var <- Sd**2
#'   cat('Calculating CV\n')
#'   CV <- Sd/Mean
#'
#'   # cat('Calculating PresentCall\n')
#'   # if (length(presentCall.cutoff) == 0) {
#'   #   presentCall.cutoff = c(0,1)
#'   # }
#'   # pC <- pbsapply(presentCall.cutoff, simplify = T,
#'   #                function(cutoff) {
#'   #                  apply(exprs,1,function(x) {sum(x > (cutoff + min(x)))/length(x)})
#'   #                })
#'   # colnames(pC) <- paste('pC',presentCall.cutoff,sep = '_')
#'
#'   df <- data.frame(Mean = Mean,
#'                    Sd  = Sd,
#'                    Var = Var,
#'                    CV  = CV)
#'   # df <- cbind(df,pC)
#'   rownames(df) <- rownames(exprs)
#'
#'   return(df)
#' }

#' #' scaleFactor
#' #'
#' #' @Description:
#' #' Reduce data to a fixed interval
#' #'
#' #' Arguments:
#' #' @param x orig data, a vector
#' #' @param min min value of interval
#' #' @param max max value of interval
#' #'
#' #' @export
#' scaleFactor <- function(x, min = 0, max = 1){
#'   k <- (max - min)/(max(x) - min(x) )
#'   new <- min + k * (x - min(x))
#'
#'   return(new)
#' }


#' ObjectFunction
#'
#' @description
#' The objective function for calculating trajectory scores
#'
#' @param par Initial values for the parameters to be estimated
#' @param Yx Matrix required for calculating trajectory scores
#' @param S_list a list, covariance matrix for all trajectories
#' @param debug_file log files
#' @export
ObjectFunction = function(par, Yx, S_list, debug_file) {

  L = nrow(Yx)
  P = ncol(Yx)
  R = as.SO3(par[1:L])

  s = par[L+1]
  W = matrix(R, nrow = L, ncol = L) %*% Yx[1:L, ]
  E_tmp = t(matrix(par[-c(1:(L+1))],nrow=P))

  lik <- sapply(1:L, function(l){
    cov_matrix = S_list[[l]]
    E = E_tmp[l,]
    Ws = W[l, ] + E
    PR = dmvnorm(Ws, mean = rep(0, P), sigma = cov_matrix, log = TRUE)
    PE = dnorm(E, mean = 0, sd = s, log = TRUE)
    return(-PR - sum(PE))
  })

  Ps = dinvgamma(s^2, shape = 0.01, scale = 0.01, log = TRUE)

  write.table(paste0("liklyhood",sum(lik) - Ps),
              file = paste0(getwd(),"/3_tracking/3.1_optim_result/",debug_file),
              quote = FALSE, append = TRUE,
              row.names = FALSE, col.names = FALSE, sep = "\n")

  return(sum(lik) - Ps)
}

#' RunningmodMGPtrack
#'
#' @description calculating trajectory scores
#' @param object MGPfact object
#' @param iter optim iterations
#' @param cores number of threads
#' @param seed random seed
#' @export
#'
RunningmodMGPtrack <- function(object,
                               iter = 100,
                               cores = 1,
                               seed = 723){
  sdf = GetMURPInfo(object)
  S_list = object@GPR$murp_cov$cov_l
  P = getParams(object, "murp_number")
  L = getParams(object, "trajectory_number")
  Q = getParams(object, "murp_pc_number")
  Yx = t(object@MURP$centersPCA$x[,1:L])

  ## ordering
  ord = order(sdf$T)
  Yx = Yx[,ord]
  for(l in 1:L){S_list[[l]] = S_list[[l]][ord,ord] }

  ## updating S_list
  new_S_list <- lapply(S_list, function(S){
    s = svd(S)
    eig_vals = s$d
    eig_vals[which(eig_vals<1)] = 0
    s$u %*% diag(s$d) %*% t(s$v)
  })

  #### parallel
  # iter = 100
  # cores = 1
  # debug_file = paste0("3_tracking/3.1_optim_result/debug_iter_",iter,".txt")
  debug_file = paste0("debug_iter_",iter,".txt")
  CL <- makeCluster(cores)
  setDefaultCluster(cl = CL)
  clusterEvalQ(CL,library(invgamma))
  clusterEvalQ(CL,library(dplyr))
  clusterEvalQ(CL,library(mvtnorm))
  clusterEvalQ(CL,library(rotations))
  system(paste0("rm ",debug_file))

  t1 = Sys.time()

  set.seed(seed)
  init_value = rnorm(L+1+P*L, mean = 0, sd = 1)
  lower = rep(-Inf, L+1+P*L)
  lower[L+1] = 1e-8

  set.seed(seed)
  m <- optimParallel(par = init_value,
                     fn = ObjectFunction,
                     Yx = Yx,
                     S_list = S_list,
                     debug_file = debug_file,
                     method = "L-BFGS-B",
                     control = list(maxit=iter),
                     parallel = list(loginfo = TRUE),
                     lower = lower,
                     upper = rep(Inf, L+1+P*L))
  t2 = Sys.time()
  setDefaultCluster(cl = NULL)
  stopCluster(CL)
  time = RunTime(t1,t2)
  optim_result = list(optim = m,
                      init = init_value,
                      execution_time = time)
  save(optim_result, file = paste0("3_tracking/3.1_optim_result/iter",iter,"_optim.rda"))

  ## import track optim result to object
  object = assignSettings(object, "track_optim_iterations", iter)
  object = ImportTrackResult(object, optim_result = optim_result, hasinits = TRUE)

  # command
  command = GetCommand()
  object@Command$track$RunningmodMGPtrack = command

  return(object)
}

#' ImportOptimR
#'
#' @description import trajectory scores in mgpfact object
#' @param object MGPfact object
#' @param optim_result tracking result from optim
#' @param hasinits Logical value, whether has initial values
#' @export
ImportTrackResult <- function(object,
                              optim_result = NULL,
                              hasinits = TRUE){

  if(is.null(optim_result)){
    iter = getParams(object, "track_optim_iterations")
    load(paste0("3_tracking/3.1_optim_result/iter",iter,"_optim.rda"))
  }
  m = optim_result$optim
  init_value = optim_result$init

  P = getParams(object, "murp_number")
  L = getParams(object, "trajectory_number")
  Q = getParams(object, "murp_pc_number")

  values1 = array(m$loginfo[,2:(L+1+L*P+1)], dim = c(nrow(m$loginfo), L+1+L*P, 1))
  names = lapply(1:L, function(i) paste0("E[", i, ", ",1:P,"]")) %>% unlist
  iter = nrow(m$loginfo)
  burnin = 0
  chains = 1

  if(hasinits){
    inits = list(rot_seed = init_value[1:3], s2 = 1,
                 E = matrix(init_value[5:length(init_value)], ncol= L, nrow= P))

    OptimResult <- CreateTrackOptimResultObject(sample_value = values1,
                                                dimnames = list(paste0("iter",1:dim(values1)[1]),
                                                                names = c(paste0("rot_seed[",1:L,"]"), "s", names),
                                                                paste0("chain",1:dim(values1)[3])),
                                                burnin = burnin,
                                                chains = chains,
                                                iter = iter,
                                                hasinits = hasinits,
                                                inits = inits)
  }else{
    OptimResult <- CreateTrackOptimResultObject(sample_value = values1,
                                                dimnames = list(paste0("iter",1:dim(values1)[1]),
                                                                names = c(paste0("rot_seed[",1:L,"]"), "s", names),
                                                                paste0("chain",1:dim(values1)[3])),
                                                burnin = burnin,
                                                chains = chains,
                                                iter = iter,
                                                hasinits = hasinits )
  }

  save(OptimResult, file = "3_tracking/3.1_optim_result/OptimResult.rda")
  object@OptimResult$track = OptimResult
  assignSettings(object, "track_optim_iterations", object@OptimResult$track@iter)

  return(object)
}

#' GetTrackSdf
#'
#' @description
#' Organize the results obtained from 'RunningmodMGPtrack'
#' @param object MGPfact object
#' @param iter_range Retain the results of parameter optimization
#' within a specified iteration range
#' @param save logical value, whether to save locally
#'
#' @export
#'
GetTrackSdfFromOptim <- function(object = NULL,
                                 iter_range = NULL,
                                 save = TRUE){

  if(is.null(iter_range)){ iter_range = c(object@OptimResult$track@iter,
                                          object@OptimResult$track@iter)}
  pse_sdf = GetMURPInfo(object)
  pse_sdf = pse_sdf[,1:which(colnames(pse_sdf)=="unified_direction")]
  L = getParams(object, "trajectory_number")
  P = getParams(object, "murp_number")

  ## 1. get meanvalue
  param_meanvalue <- GetAllParamMeanChain(object = object,
                                          iter_range = iter_range,
                                          aspect = "track")
  R = matrix(as.SO3(param_meanvalue$rot_seed[1,]),ncol = L, nrow = L)

  ## 2. ordering
  ord = order(pse_sdf$T)
  Yx = pse_sdf[,paste0("PC_",1:L)] %>% t %>% as.matrix
  Yx = Yx[,ord]

  W = as.matrix(R) %*% Yx
  E = matrix(apply(param_meanvalue$E, 2, mean), ncol= L) %>% t
  Ws = t(W + E)
  colnames(Ws) = paste0("W_",1:L)

  ## 3. create track_sdf
  track_sdf <- data.frame(pse_sdf[ord,], Ws)
  tmplist <- lapply(1:L, function(l){

    tb = paste0("Tb_", l)
    c = paste0("C0_", l)
    w = paste0("W_", l)
    t = "T"
    df = track_sdf[, c(t, w, c, tb)]

    df$C_1 = NA
    df$C_2 = NA
    df$C_1[which(df[,c]==0 | df[,c]==1)] = 1
    df$C_2[which(df[,c]==0 | df[,c]==2)] = 2
    df$C_1 = factor(df$C_1, levels = c(0,1,2))
    df$C_2 = factor(df$C_2, levels = c(0,1,2))

    df$W_1_1 = as.numeric(as.character(df$C_1))
    df$W_1_2 = as.numeric(as.character(df$C_2))
    df[which(df$W_1_1==1), "W_1_1"] = df[which(df$W_1_1==1), w]
    df[which(df$W_1_2==2), "W_1_2"] = df[which(df$W_1_2==2), w]

    df = df[,c("C_1","C_2","W_1_1","W_1_2")]
    colnames(df) = c(paste0("C0_",l,"_",1:2),
                     paste0("W_",l,"_",1:2) )
    df
  })
  tmp = data.frame(do.call(cbind, tmplist), Ws)
  rownames(tmp) = rownames(track_sdf)
  tmp = tmp[paste0("T[",1:P,"]"),]
  object = AddMURPMetadata(object, tmp) # add metadata
  object = assignSettings(object,"weight_iter_range",iter_range) ## set iter_range for weight
  object = ComputeGeneWeight(object)

  ## 4. create orig_track_sdf
  # orig_sdf = ct@MetaData
  # w_names <- colnames(track_sdf)[grep("^W", colnames(track_sdf))]
  # tmp <- lapply(w_names,  function(i){
  #   GetPseudoTt(CellNum = getParams(object, "cell_number"),
  #               pred_t = track_sdf[,i],
  #               Cluster_Label = object@MURP$Recommended_K_cl$cluster) %>% as.numeric
  # })

  w_names <- colnames(tmp)[grep("^W", colnames(tmp))]
  orig_tmp <- lapply(w_names,  function(i){
    MapMURPLabelToAll(vecc = tmp[,i],
                      orig = object@MURP$Recommended_K_cl$cluster) %>% as.numeric
  })
  w  = do.call(cbind, orig_tmp) %>% data.frame
  colnames(w) = w_names
  object = AddMetadata(object, w) # assign metadata

  ## 5. save
  if(save){
    sdf = GetMURPInfo(object)
    orig_sdf = object@MetaData
    save(sdf, file = paste0("3_tracking/sdf_",iter_range[1],"_",iter_range[2],".rda"))
    save(orig_sdf, file = paste0("3_tracking/orig_sdf_",iter_range[1],"_",iter_range[2],".rda"))
  }

  ## 6. command
  command = GetCommand()
  object@Command$track$GetTrackSdfFromOptim = command

  return(object)
}

#' ComputeWeightGene
#'
#' @description Calculate gene weights across different trajectories
#' @param object MGPfact object
#' @export
#'
ComputeGeneWeight <- function(object = NULL){

  # if(is.null(iter_range)){ iter_range = getParams(object,"weight_iter_range") }
  iter_range = getParams(object,"weight_iter_range")
  P = getParams(object, "murp_number")
  L = getParams(object, "trajectory_number")
  Q = getParams(object, "murp_pc_number")

  ## 1. get meanvalue and track_sdf
  # load(paste0("3_tracking/sdf_",iter_range[1], "_",iter_range[2],".rda"))
  # load(paste0("3_tracking/param_meanvalue_",iter_range[1], "_", iter_range[2],".rda"))
  param_meanvalue <- GetAllParamMeanChain(object = object,
                                          iter_range = iter_range,
                                          aspect = "track")

  ## 2. get R
  if("R" %in% names(param_meanvalue)){
    ## 为了满足旋转矩阵的特性，这部分需要每条chain的R结果单独得到weight
    weight_list <- lapply(1:nrow(param_meanvalue$R),function(i){
      R <- matrix(param_meanvalue$R[i,], ncol= L) %>% as.data.frame
      rownames(R) <- paste0("R_", 1:L)
      R <- as.matrix(R)
      U <- object@MURP$centersPCA$rotation[,1:L]
      weight <- U %*% R
      return(weight)
    })
    weight = Reduce("+",weight_list)/nrow(param_meanvalue$R)
  }else{
    R = matrix(as.SO3(param_meanvalue$rot_seed[1,]),ncol = L, nrow = L)
    U <- object@MURP$centersPCA$rotation[,1:L]
    weight <- U %*% R
  }

  ## 3. gene weight
  # U <- object@MURP$centersPCA$rotation[,1:L]
  # weight <- U %*% R

  ## 4. gene weight order
  gene_weight_order <- apply(abs(weight),2,function(x){
    rownames(weight)[order(x,decreasing = TRUE)] })

  ## 5. save to object
  result = list(weight = weight)
  names = paste0("iter_", iter_range[1], "_", iter_range[2])
  gene_weight <- object@OptimResult$track@gene_weight
  if(names %in% names(gene_weight)){
    gene_weight[[names]] = result
  }else{
    len = length(gene_weight)
    gene_weight[[len+1]] = result
    names(gene_weight)[len+1] <- paste0("iter_", iter_range[1], "_", iter_range[2])
  }

  object@OptimResult$track@gene_weight = gene_weight

  ## 6. save
  # save(weight, file = paste0("3_tracking/3.3_gene_weight/gene_weight_",iter_range[1],"_",iter_range[2],".rda"))
  # write.csv(gene_weight_order, file = paste0("3_tracking/3.3_gene_weight/gene_weight_order_", iter_range[1],"_",iter_range[2],".csv") )
  # command = GetCommand()
  # object@Command$track$ComputeWeightGene = command

  return(object)
}


#' WeightGeneFilter
#'
#' @description Filter gene weights
#'
#' @param object MGPfact object
#' @param iter_range Retain the results of parameter optimization
#' within a specified iteration range
#' @param method filtering methods
#' c("top_ratio","quantile","confidence","top_number","weight_cut")
#' @param q_seq parameters for the quantile method: quantile intervals
#' @param top_ratio parameters for the top_ratio: select genes with the highest weights based on proportion
#' @param top_number patameters for the top_number: select genes with the highest weights based on quantity
#' @param weight_cut patameters for the weight_cut: select genes with the highest weights based on weight
#' @param conf.int patameters for the confidence: select genes with the highest weights based on confidence
#' @export
#'
WeightGeneFilter<- function(object,
                            iter_range = NULL,
                            method = c("top_ratio","top_number","weight_cut"),
                            conf.int = 0.95,
                            q_seq = 0.25,
                            top_ratio = 0.1,
                            top_number = 500,
                            weight_cut = 0.03){

  if(is.null(iter_range)){ iter_range = getParams(object,"weight_iter_range") }
  index = paste0("iter_", iter_range[1], "_", iter_range[2])
  weight = object@OptimResult$track@gene_weight[[index]]$weight
  L = getParams(object, "trajectory_number")

  if(method=="top_ratio"){
    object = assignSettings(object, paste0("weight_",method), top_ratio)
  }
  # if(method=="quantile"){
  #   object = assignSettings(object, paste0("weight_",method), q_seq)
  # }
  if(method=="top_number"){
    object = assignSettings(object, paste0("weight_",method), top_number)
  }
  # if(method=="confidence"){
  #   object = assignSettings(object, paste0("weight_",method), conf.int)
  # }
  if(method=="weight_cut"){
    object = assignSettings(object, paste0(method), weight_cut)
  }

  # set filename
  # if(method == "confidence"){
  #   # file_name = paste0("confidence_", getParams(object, "weight_conf.int"))
  #   file_name = paste0("confidence_", conf.int)
  # }
  # if(method == "quantile"){
  #   # file_name = paste0("quantile_", getParams(object, "weight_q_seq"))
  #   file_name = paste0("quantile_", q_seq)
  # }
  if(method == "top_ratio"){
    # file_name = paste0("top_ratio_", getParams(object, "weight_top_ratio"))
    file_name = paste0("top_ratio_", top_ratio)
  }
  if(method == "top_number"){
    # file_name = paste0("top_ratio_", getParams(object, "weight_top_ratio"))
    file_name = paste0("top_number_", top_number)
  }
  if(method == "weight_cut"){
    # file_name = paste0("top_ratio_", getParams(object, "weight_top_ratio"))
    file_name = paste0("weight_cut_", weight_cut)
  }

  ## 0. function
  rf = function(w,
                method = method,
                q_seq = q_seq,
                conf.int = conf.int,
                top_ratio = top_ratio,
                weight_cut = weight_cut){

    if(method == "confidence"){
      t = t.test(w, conf.level=conf.int)
      q2 = t$conf.int[1]
      q4 = t$conf.int[2]
    }
    if(method == "quantile"){
      q2_seq = q_seq
      q4_seq = 1-q_seq
      q2 = qnorm(q2_seq, mean=mean(w),sd=sd(w))
      q4 = qnorm(q4_seq, mean=mean(w),sd=sd(w))
    }
    if(method == "top_ratio"){
      gene_num = round(length(w) * top_ratio)
      w2 = w[order(abs(w), decreasing=T)]
      tmp_g = names(w2)[1:gene_num]

      status = names(w)
      status[which(w<0)] = 1
      status[which(w>0)] = 2
      status[which(names(w) %nin% tmp_g)] <- 3

      # x_label
      s1 = w[which(w<0)]
      s2 = w[which(w>0)]

      g1 = names(w)[status == 1]
      g2 = names(w)[status == 2]
      x1 = s1[ which(s1>max(s1[g1])) ] %>% min
      x2 = s2[ which(s2<min(s2[g2])) ] %>% max
      q2 = x1 + (abs(x1) - abs(max(s1[g1])) )/2
      q4 = x2 + (min(s2[g2]) - x2)/2
    }
    if(method == "top_number"){
      gene_num = top_number
      w2 = w[order(abs(w), decreasing=T)]
      tmp_g = names(w2)[1:gene_num]

      status = names(w)
      status[which(w<0)] = 1
      status[which(w>0)] = 2
      status[which(names(w) %nin% tmp_g)] <- 3

      # x_label
      s1 = w[which(w<0)]
      s2 = w[which(w>0)]

      g1 = names(w)[status == 1]
      g2 = names(w)[status == 2]
      x1 = s1[ which(s1>max(s1[g1])) ] %>% min
      x2 = s2[ which(s2<min(s2[g2])) ] %>% max
      q2 = x1 + (abs(x1) - abs(max(s1[g1])) )/2
      q4 = x2 + (min(s2[g2]) - x2)/2
    }
    if(method == "weight_cut"){
      tmp_g = names(w)[which(abs(w)>weight_cut)]

      status = names(w)
      status[which(w<0)] = 1
      status[which(w>0)] = 2
      status[which(names(w) %nin% tmp_g)] <- 3

      # x_label
      s1 = w[which(w<0)]
      s2 = w[which(w>0)]

      g1 = names(w)[status == 1]
      g2 = names(w)[status == 2]
      x1 = s1[ which(s1>max(s1[g1])) ] %>% min
      x2 = s2[ which(s2<min(s2[g2])) ] %>% max
      q2 = x1 + (abs(x1) - abs(max(s1[g1])) )/2
      q4 = x2 + (min(s2[g2]) - x2)/2
    }
    x = list(q2 = q2, q4 = q4)
    return(x)
  }

  ## 1. PLOT: qq + CI
  df <- weight %>% data.frame
  colnames(df) <- paste0("L", 1:L)
  plist <- lapply(1:L, function(i){
    # cat(i, "\n")
    c = paste0("L",i)
    w = weight[,i]

    # cat("q2: ", q2, " q4: ", q4, "\n")
    qq = rf(w, method = method, q_seq = q_seq,  conf.int = conf.int, top_ratio = top_ratio, weight_cut = weight_cut)
    q2 = qq$q2; q4 = qq$q4

    # plot
    ggplot(df, aes(sample = get(c))) +
      geom_qq(alpha = 0.7) +
      geom_qq_line(fullrange = TRUE) +
      geom_hline(yintercept = q2, linetype = "dashed", color = "blue") +
      geom_hline(yintercept = q4, linetype = "dashed", color = "blue") +
      labs(x = "Theoretical Quantiles", y = "Sample Quantiles") +
      guides(colour = FALSE, fill = FALSE) +
      rj.ftheme
  })
  p = wrap_plots(plist, ncol = L)

  ggsave(paste0("3_tracking/3.3_gene_weight/weight_qq_ci_",
                iter_range[1], "_", iter_range[2], "_", file_name, ".pdf"), p,
         width = 12, height = ceiling(L/3)*3.5, limitsize = FALSE)

  ## PLOT: hist + CI
  plist <- lapply(1:L, function(i){
    # cat(i, "\n")
    c = paste0("L",i)
    w = weight[,i]

    qq = rf(w, method = method, q_seq = q_seq,  conf.int = conf.int, top_ratio = top_ratio,weight_cut = weight_cut)
    q2 = qq$q2; q4 = qq$q4

    ggplot(df, aes(x = get(c)),  color = "grey", fill = "grey") +
      geom_histogram(position="identity", binwidth = 0.002, alpha = 0.9) +
      geom_vline(xintercept = q2, linetype = "dashed", color = "blue") +
      geom_vline(xintercept = q4, linetype = "dashed", color = "blue") +
      labs(x = "Gene Weight", y = "Frequency") +
      guides(colour = FALSE, fill = FALSE) +
      rj.ftheme
  })
  p = wrap_plots(plist, ncol = L)
  ggsave(paste0("3_tracking/3.3_gene_weight/weight_hist_ci_",
                iter_range[1], "_", iter_range[2], "_", file_name, ".pdf"), p,
         width = 12, height = ceiling(L/3)*3.5, limitsize = FALSE)


  # 利用所有轨迹的weight计算
  # x = as.vector(weight)
  # qt = quantile(x, probs = seq(0, 1, q_seq))
  # q2 = qt[2]
  # q4 = tail(qt,2)[1]
  # cat("q2: ", q2, " q4: ", q4, "\n")
  weight_filter <- lapply(1:L, function(i){

    cat("Trajectory ", i, "\n")
    w = weight[,i]

    # 直接利用quantile计算
    # qt = quantile(x, probs = seq(0, 1, q_seq))
    # q2 = qt[2]
    # q4 = tail(qt,2)[1]
    # cat("q2: ", q2, " q4: ", q4, "\n")

    # 利用qnorm计算
    # 如qnorm(0.05)=-1.64,即x<=-1.64小于0.05
    # q2 = qnorm(q2_seq,mean=mean(w),sd=sd(w))
    # q4 = qnorm(q4_seq,mean=mean(w),sd=sd(w))
    # cat("q2: ", q2, " q4: ", q4, "\n")

    # 利用置信计算
    qq = rf(w, method = method, q_seq = q_seq,  conf.int = conf.int, top_ratio = top_ratio,weight_cut = weight_cut)
    q2 = qq$q2; q4 = qq$q4

    g = rownames(weight)
    ind = which(w <= q2 | w >= q4)
    g[ind]
  })

  save(weight_filter,
       file = paste0("3_tracking/3.3_gene_weight/weight_filter_qq_",
                     iter_range[1], "_", iter_range[2],
                     "_", file_name, ".rda") )

  # save
  tmp = object@OptimResult$track@gene_weight[[index]]
  if(file_name %in% names(tmp)){
    tmp[[file_name]] = weight_filter
  }else{
    len = length(tmp)
    tmp[[len+1]] = weight_filter
    names(tmp)[len+1] <- file_name
  }
  object@OptimResult$track@gene_weight[[index]] = tmp

  return(object)
}

#' GetWeightGene
#'
#' @param object MGPfact object
#' @param method filtering methods
#' c("top_ratio","quantile","confidence","top_number","weight_cut")
#' @param trajectory trajectory index
#'
#' @export
#'
GetWeightGene <- function(object,
                          method = c("weight_cut","top_ratio","top_number"),
                          trajectory){

  track = getTrack(object)

  # conf.int = getParams(object, "weight_conf.int")
  # confidence = getParams(object, "weight_confidence")
  # q_seq = getParams(object, "weight_q_seq")
  # top_ratio = getParams(object, "weight_top_ratio")

  iter_range = getParams(object, "weight_iter_range")
  index = paste0("iter_", iter_range[1], "_", iter_range[2])
  t = track@gene_weight[[index]]

  # if(method == "confidence"){
  #   index2 = paste0("confidence_", getParams(object, "weight_conf.int"))
  # }
  # if(method == "quantile"){
  #   index2 = paste0("quantile_", getParams(object, "weight_q_seq"))
  # }
  if(method == "top_ratio"){
    index2 = paste0("top_ratio_",getParams(object, "weight_top_ratio"))
  }
  if(method == "top_number"){
    index2 = paste0("top_number_",getParams(object, "weight_top_number"))
  }
  if(method == "weight_cut"){
    index2 = paste0("weight_cut_",getParams(object, "weight_cut"))
  }
  genes = t[[index2]][[trajectory]]
  return(genes)
}

#' GetGeneWeight
#' @description Obtain gene weights for each trajectory
#' @param object MGPfact object
#'
#' @export
#'
GetGeneWeight <- function(object){
  track = getTrack(object)
  index1 = paste0("iter_", getParams(object, "weight_iter_range")[1],
                  "_", getParams(object, "weight_iter_range")[2])
  # if(getParams(object, "weight_confidence")){
  #   index_name = paste0("filter_confidence_", getParams(ct, "weight_conf.int") )
  # }else{
  #   index_name = paste0("filter_qseq_", getParams(ct, "weight_q_seq") )
  # }
  w = track@gene_weight[[index1]][["weight"]]
  return(w)
}

#' WriteWeightGene
#'
#' @param object MGPfact object
#' @param method filtering methods
#' c("top_ratio","quantile","confidence","top_number","weight_cut")
#' @export
WriteWeightGene <- function(object,
                            method = c("weight_cut","top_ratio","top_number")){
  gw = GetGeneWeight(object)

  num_trajectories = getParams(ct, "trajectory_number")

  wt <- lapply(1:num_trajectories, function(i){
    x = GetWeightGene(ct, method = method, trajectory = i)
    w = gw[x, i] %>% abs
    w = sort(w, decreasing = TRUE)
    w = data.frame(gene = names(w), weight = w %>% round(3))
    return(w)
  })

  max_rows = max(sapply(wt, nrow))
  wtm = matrix("", ncol = num_trajectories * 2, nrow = max_rows)
  colnames(wtm) = unlist(lapply(1:num_trajectories, function(i) paste0(c("Trajectory ", "Weight "), i)))

  for(i in 1:num_trajectories){
    wtm[1:nrow(wt[[i]]), (i * 2 - 1):(i * 2)] = wt[[i]][1:nrow(wt[[i]]), 1:2] %>% as.matrix
  }

  iter_range = object@Settings@settings$weight_iter_range
  output_file = paste0("4_differential_genes/gene_weight_", method, "_",
                       paste0(iter_range, collapse = "_"), ".csv")
  write.csv(wtm, row.names = FALSE, file = output_file)
  return(wtm)
}

#' TrajPlot
#' @description
#' smooth trajectory visualization
#' @param object MGPfact object
#' @param col the colors of different branches
#' @param save logical value, whether to save pdf
#' @param box_fill box color
#' @param title plot title
#' @param g_title
#' @param tb logical value, whether label bifurcation line
#' @param tb_pred bifurcation point
#' @param pointSize point size
#' @param pointAlpha point alpha
#' @param pointLabel logical value, whether label point
#' @param pointLabelsize label size
#' @param rug logical value, whether to plot rug
#' @param rugAlpha rug alpha
#' @param legend locical value, whether to add legend
#' @param legend_title legend title
#' @param lm whether add lm line
#' @param lineMethod linear model methods
#' @param lineType line type
#' @param lineSize line size
#' @param lineAlpha line alpha

#' @export
#'
TrajPlot <- function(object ,
                     box_fill = NULL,
                     plot_title = TRUE,
                     g_title = NULL,
                     col = NULL,
                     plot_tb = TRUE,
                     tb_pred = NULL,
                     pointSize = 5,
                     pointAlpha = 0.4,
                     pointLabel = FALSE,
                     pointLabelsize = 3,
                     rug = TRUE,
                     rugAlpha = 0.3,
                     legend = TRUE,
                     legend_title = NULL,
                     lm = FALSE,
                     lineMethod = loess,
                     lineType = "solid",
                     lineSize = 1,
                     lineAlpha = 0.3,
                     se = F,
                     span = 0.75,
                     formula = y~x,
                     family = NULL,
                     vlineType = "dashed",
                     ...){

  coln = ifelse(is.null(col), "C", col)
  L = getParams(object, "trajectory_number")

  plist = list()
  for(indexL in 1:L){
    df = GetMURPInfo(object)
    tb_pred = df[1,paste0("Tb_", indexL)]
    c = paste0("C0_", indexL)
    w = paste0("W_", indexL)
    # if(is.null(c)) { c = paste0("C0_", indexL)}
    # if(is.null(w)) { w = paste0("W_", indexL)}
    df[,paste0(c,"_1")] = as.factor(df[,paste0(c,"_1")])
    df[,paste0(c,"_2")] = as.factor(df[,paste0(c,"_2")])
    df[,c] = as.factor(df[,c])
    df[,paste0(c,"_1")] = as.factor(df[,paste0(c,"_1")])
    df[,paste0(c,"_2")] = as.factor(df[,paste0(c,"_2")])
    t = "T"
    if(coln=="C"){
      col = c
    }
    if(is.null(box_fill)){
      if(col==c){
        box_fill = colorRampPalette(pal_d3("category10")(10))(10)[c(8,1,4)]
        names(box_fill) = c(0,1,2)
      } else{
        box_fill = colorRampPalette(pal_d3("category20")(20))(20)
      }
    }

    p = ggplot(df) +
      geom_point(alpha = pointAlpha, size = pointSize,
                 # aes(x = get(t), y = get(w), colour = factor(get(col)))
                 aes_string(x = t, y = w, colour = col ))

    # whether tb
    if(plot_tb){
      p = p +
        geom_vline(xintercept = tb_pred, colour = "#990000", linetype = "dashed")
    }

    # rug
    if(rug){
      p = p +
        geom_rug(alpha = rugAlpha,
                 aes_string(x = t, y = paste0(w,"_1"), colour = paste0(c,"_1")) ) +
        geom_rug(alpha = rugAlpha,
                 aes_string(x = t, y = paste0(w,"_2"), colour = paste0(c,"_2")) )
    }

    # loess
    if(lm){
      p = p +
        geom_line( stat = "smooth",
                   aes_string(x = t, y = paste0(w,"_1"), colour = paste0(c,"_1")),
                   na.rm = TRUE,
                   method = lineMethod,
                   linetype = lineType,
                   size = lineSize,
                   alpha = lineAlpha,
                   se = se,
                   span = span,
                   family = family,
                   formula = formula) +
        geom_line( stat = "smooth",
                   aes_string(x = t, y = paste0(w,"_2"), colour = paste0(c,"_2")),
                   na.rm = TRUE,
                   method = lineMethod,
                   linetype = lineType,
                   size = lineSize,
                   alpha = lineAlpha,
                   se = se,
                   span = span,
                   family = family,
                   formula = formula)

    }

    # label
    if(pointLabel){
      p <- p +
        geom_text(aes_string(x = t, y = w, label = "names"),
                  color = "black",
                  size = pointLabelsize)
    }

    # final
    if(!plot_title){
      p = p +
        scale_colour_manual(values = box_fill) +
        labs(title = paste0(" "), x = "PseudoT", y = "Score", color = legend_title) +
        rj.ftheme
    }else{
      if(is.null(g_title)){
        p = p +
          scale_colour_manual(values = box_fill, drop = F) +
          labs(title = paste0("Trajectory ",indexL), x = "PseudoT", y = "Score", color = legend_title) +
          rj.ftheme
      }else {
        p = p +
          scale_colour_manual(values = box_fill) +
          labs(title = paste0(g_title," in trajectory ",indexL), x = "PseudoT", y = "Score", color = legend_title) +
          rj.ftheme
      }
    }

    # legend
    if(!legend){
      p = p + guides(colour = "none")
    }
    plist = append(plist, list(p))
  }
  pl <- wrap_plots(plist, ncol = L, guides = "collect")
  if(save){
    ggsave(paste0("3_tracking/3.2_trajectory/3_t_score_mean_",coln,".pdf"),
           pl, width = L*4, height = 4)
  }else{
    return(pl)
  }
}


# --------------

#' TrajGenePlot
#'
#' @description
#' The expression trend of any gene on a smooth trajectory
#'
#' @param df dataframe contain plot information
#' @param pointColorValue point color
#' @param traj trajectory index
#' @param gene gene name
#' @param title plot title
#' @param g_title
#' @param col the colors of different branches
#' @param c branching column
#' @param w trajectory score column
#' @param tb logical value, whether label bifurcation line
#' @param tb_pred bifurcation point
#' @param pointSize point size
#' @param pointAlpha point alpha
#' @param pointLabel logical value, whether label point
#' @param pointLabelsize label size
#' @param rug logical value, whether to plot rug
#' @param rugAlpha rug alpha
#' @param legend locical value, whether to add legend
#' @param legend_title legend title
#' @param lm whether add lm line
#' @param lineMethod linear model methods
#' @param lineType line type
#' @param lineSize line size
#' @param lineAlpha line alpha
#' @export
#'
# TrajGenePlot <- function(df = NULL,
#                          pointColorValue = NULL,
#                          expr_methods = "mean",
#                          traj = 1,
#                          gene = NULL,
#                          c = NULL,
#                          w = NULL,
#                          tb = TRUE,
#                          tb_pred = NULL,
#                          pointSize = 4,
#                          pointAlpha = 0.5,
#                          pointLabel = FALSE,
#                          pointLabelsize = 3,
#                          rug = TRUE,
#                          rugAlpha = 0.3,
#                          legend = FALSE,
#                          lm = TRUE,
#                          lineMethod = "loess",
#                          lineType = "solid",
#                          lineSize = 1,
#                          lineAlpha = 0.8,
#                          se = F,
#                          span = 0.75,
#                          formula = y~x,
#                          vlineType = "dashed",
#                          ...){
#
#   require(ggplot2)
#   require(cowplot)
#   require(ggsci)
#
#   l = traj
#   df = df
#   t = "T"
#   df[,paste0(c,"_1")] = as.factor(df[,paste0(c,"_1")])
#   df[,paste0(c,"_2")] = as.factor(df[,paste0(c,"_2")])
#   df[,c] = as.factor(df[,c])
#   col = gene
#   # if(is.null(col)){
#   #   col = c
#   # }
#   # if(is.null(pointColorValue)){
#   #   if(col==c){
#   #     pointColorValue = colorRampPalette(pal_d3("category10")(10))(10)[c(8,1,4)]
#   #     names(pointColorValue) = c(0,1,2)
#   #   } else{
#   #     pointColorValue = colorRampPalette(pal_d3("category20")(20))(20)
#   #   }
#   # }
#
#   p = ggplot(df) +
#     geom_point(alpha = pointAlpha, size = pointSize,
#                aes_string(x = t, y = w, colour = col )) +
#     scale_colour_gradientn(colours = pointColorValue) +
#     labs(title = paste0(gsub("\\.","-",col)," in Trajectory",l),
#          # subtitle = paste0("SRCC between pseudot and ", col, ": ",
#          #                   round(cor(df[,col],df$T,method = "spearman"),4)),
#          subtitle = paste0("SRCC: ",
#                            round(df[1, paste0(col,"_cor")],4)),
#          x = "PseudoT", y = "Score",
#          colour = expr_methods) +
#     rj.ftheme
#
#   # whether tb
#   if(tb){
#     p = p +
#       geom_vline(xintercept = tb_pred, colour = "#990000", linetype = "dashed")
#   }
#
#   # rug
#   if(rug){
#     p = p +
#       geom_rug(alpha = rugAlpha,
#                aes_string(x = t, y = paste0(w,"_1"), colour = paste0(c,"_1")) ) +
#       geom_rug(alpha = rugAlpha,
#                aes_string(x = t, y = paste0(w,"_2"), colour = paste0(c,"_2")) )
#   }
#
#   # loess
#   if(lm){
#     p = p +
#       geom_line( stat = "smooth",
#                  data = df,
#                  aes_string(x = t, y = paste0(w,"_1")),
#                  colour = "#D62728",
#                  na.rm = TRUE,
#                  method = lineMethod,
#                  linetype = lineType,
#                  size = lineSize,
#                  alpha = lineAlpha,
#                  se = se,
#                  span = span,
#                  formula = formula) +
#       geom_line( stat = "smooth",
#                  data = df,
#                  aes_string(x = t, y = paste0(w,"_2")),
#                  colour = "#1F77B4",
#                  na.rm = TRUE,
#                  method = lineMethod,
#                  linetype = lineType,
#                  size = lineSize,
#                  alpha = lineAlpha,
#                  se = se,
#                  span = span,
#                  formula = formula)
#   }
#
#   # label
#   if(pointLabel){
#     p <- p +
#       geom_text(aes_string(x = t, y = w, label = "names"),
#                 color = "black",
#                 size = pointLabelsize)
#   }
#
#   # legend
#   if(!legend){
#     p = p + guides(colour = "none")
#   }
#
#   return(p)
# }
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' GetTrackSdf
#'
#' @description
#'
#' @param object MGPfact object
#' @param iter_range
#' @param pse_sdf
#'
#' @export
#'
#' GetTrackSdf <- function(object = NULL,
#'                         iter_range = NULL,
#'                         save = TRUE){
#'
#'   pse_sdf = GetMURPInfo(object)
#'   L = getParams(object, "trajectory_number")
#'
#'   ## 1. get meanvalue
#'   param_meanvalue <- GetAllParamMeanChain(mamba = object@OptimResult$track,
#'                                           murp = object@MURP,
#'                                           iter_range = iter_range,
#'                                           aspect = "track")
#'
#'   W = matrix(apply(param_meanvalue$W, 2, mean), ncol= L) %>% as.data.frame
#'   colnames(W) = paste0("W_", 1:L)
#'
#'   ## 2. create track_sdf
#'   track_sdf <- data.frame(pse_sdf, W)
#'   tmplist <- lapply(1:L, function(l){
#'
#'     tb = paste0("Tb_", l)
#'     c = paste0("C0_", l)
#'     w = paste0("W_", l)
#'     t = "T"
#'     df = track_sdf[, c(t, w, c, tb)]
#'
#'     df$C_1 = NA
#'     df$C_2 = NA
#'     df$C_1[which(df[,c]==0 | df[,c]==1)] = 1
#'     df$C_2[which(df[,c]==0 | df[,c]==2)] = 2
#'     df$C_1 = factor(df$C_1, levels = c(0,1,2))
#'     df$C_2 = factor(df$C_2, levels = c(0,1,2))
#'
#'     df$W_1_1 = as.numeric(as.character(df$C_1))
#'     df$W_1_2 = as.numeric(as.character(df$C_2))
#'     df[which(df$W_1_1==1), "W_1_1"] = df[which(df$W_1_1==1), w]
#'     df[which(df$W_1_2==2), "W_1_2"] = df[which(df$W_1_2==2), w]
#'
#'     df = df[,c("C_1","C_2","W_1_1","W_1_2")]
#'     colnames(df) = c(paste0("C0_",l,"_",1:2),
#'                      paste0("W_",l,"_",1:2) )
#'     df
#'   })
#'   tmp = do.call(cbind, tmplist) %>% as.data.frame
#'   track_sdf = cbind(track_sdf, tmp)
#'   track_sdf$names = rownames(track_sdf)
#'   object = AddMURPMetadata(object, track_sdf) # add metadata
#'   object = assignSettings(object,"weight_iter_range",iter_range) ## set iter_range for weight
#'
#'   ## 3. create orig_track_sdf
#'   w_names <- colnames(track_sdf)[grep("^W", colnames(track_sdf))]
#'   tmp <- lapply(w_names,  function(i){
#'     MapMURPLabelToAll(vecc = track_sdf[,i],
#'                       orig = object@MURP$Recommended_K_cl$cluster) %>% as.numeric
#'   })
#'   w  = do.call(cbind, tmp) %>% data.frame
#'   colnames(w) = w_names
#'   object = AddMetadata(object, w) # assign metadata
#'
#'   ## 4. save
#'   if(save){
#'     sdf = GetMURPInfo(object)
#'     orig_sdf = object@MetaData
#'     save(sdf, file = paste0("3_tracking/sdf_",iter_range[1],"_",iter_range[2],".rda"))
#'     save(orig_sdf, file = paste0("3_tracking/orig_sdf_",iter_range[1],"_",iter_range[2],".rda"))
#'   }
#'
#'   command = GetCommand()
#'   object@MURP$Command$GetTrackSdf = command
#'
#'   return(object)
#' }

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' TrajArrowPlow
#'
#' @Arguments
#' @param track_df
#' @param pointColorValue
#' @param indexL
#' @param col
#' @param c
#' @param w
#' @param tb
#' @param rug
#' @param lm
#'
# TrajArrow3Plot <- function(track_sdf = NULL,
#                            p = NULL,
#                            c = NULL,
#                            w = NULL,
#                            arrowSize = 0.15){
#
#   df = track_sdf
#   t = "T"
#   # get index of every branch
#   c0 = which(df[,c]==0)
#   c1 = which(df[,c]==1)
#   c2 = which(df[,c]==2)
#
#   # get coordinate of start and tb_s
#   if(length(c0)==0){
#     ind = which(df$T==min(df$T))
#     start = c(T = df[ind, t], W = df[ind, w])
#     tb_s = c(T = df[ind, t], W = df[ind, w])
#   }else if(length(c0)==nrow(df)){
#     ind = which(df$T==min(df$T))
#     start = c(T = df[ind, t], W = df[ind, w])
#     ind = which(df$T==max(df$T))
#     tb_s = c(T = df[ind, t], W = df[ind, w])
#   }else{
#     start = c(T = mean(df[c0, t]), W = mean(df[c0, w]) )
#     ind = which.min(abs(df$T-tb))
#     tb_s = c(T = df[ind, t], W = df[ind, w])
#   }
#
#   # get coordinate of c1_s
#   if(length(c1)==0){
#     c1_s = NULL
#   }else{
#     c1_s = c(T = mean(df[c1, t]), W = mean(df[c1, w]))
#   }
#
#   # get coordinate of c2_s
#   if(length(c2)==0){
#     c2_s = NULL
#   }else{
#     c2_s = c(T = mean(df[c2, t]), W = mean(df[c2, w]))
#   }
#
#   # plot first segment
#   # tb < min(T)
#   # tb > min(T) & c1_s/c2_s = NULL
#   # tb > min(T) & c1_s/c2_s != NULL
#   if(identical(start,tb_s)){
#     seg1 = c(start, c1_s)
#     seg2 = c(start, c2_s)
#   }else{
#     seg0 = c(start, tb_s)
#     seg1 = c(tb_s, c1_s)
#     seg2 = c(tb_s, c2_s)
#
#     seg0 = data.frame(t(seg0),color=factor(0, levels=c(0,1,2)))
#     colnames(seg0) = c("x1", "y1", "x2", "y2", "color")
#
#     if(is.null(c1_s)&is.null(c2_s)){
#       p = p + geom_segment(data = seg0, size = 1,
#                            aes(x = x1, y = y1, xend = x2, yend = y2, colour = as.factor(color)),
#                            arrow = arrow(length = unit(arrowSize, "inches"), end = "last"))
#     }else{
#       p = p + geom_segment(data = seg0, size = 1,
#                            aes(x = x1, y = y1, xend = x2, yend = y2, colour = as.factor(color)) )
#     }
#
#   }
#
#   # plot second and third segment
#   if(is.null(c1_s)&is.null(c2_s)){
#     p = p
#   }else if(!is.null(c1_s) & is.null(c2_s)){
#     seg1 = data.frame(t(seg1), color = factor(1,levels=c(0,1,2)))
#     colnames(seg1) = c("x1", "y1", "x2", "y2", "color")
#     p = p + geom_segment(data = seg1, size = 1,
#                          aes(x = x1, y = y1, xend = x2, yend = y2, colour = as.factor(color)),
#                          arrow = arrow(length = unit(arrowSize, "inches"), end = "last"))
#   }else if(is.null(c1_s) & !is.null(c2_s)){
#     seg2 = data.frame(t(seg2), color = factor(2,levels=c(0,1,2)))
#     colnames(seg2) = c("x1", "y1", "x2", "y2", "color")
#     p = p + geom_segment(data = seg2, size = 1,
#                          aes(x = x1, y = y1, xend = x2, yend = y2, colour = as.factor(color)),
#                          arrow = arrow(length = unit(arrowSize, "inches"), end = "last"))
#   }else{
#     linedf = t(data.frame(seg1, seg2)) %>% as.data.frame
#     linedf$colour = c(1,2)
#     colnames(linedf) = c("x1", "y1", "x2", "y2", "color")
#     linedf = data.frame(linedf)
#
#     p = p + geom_segment(data = linedf, size = 1,
#                          aes(x = x1, y = y1, xend = x2, yend = y2, colour = as.factor(color)),
#                          arrow = arrow(length = unit(arrowSize, "inches"), end = "last") )
#   }
#
#   return(p)
# }

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' TrajArrow2Plot
#'
#' @Arguments
#' @param track_df
#' @param pointColorValue
#' @param indexL
#' @param col
#' @param c
#' @param w
#' @param tb
#' @param rug
#' @param lm
#'
# TrajArrow2Plot <- function(track_sdf = NULL,
#                            p = NULL,
#                            c = NULL,
#                            w = NULL,
#                            arrowSize = 0.15,
#                            arrowAlpha = 0.5){
#
#   df = track_sdf
#   t = "T"
#   # get index of every branch
#   c0 = which(df[,c]==0)
#   c1 = which(df[,c]==1)
#   c2 = which(df[,c]==2)
#
#   # get coordinate of start and tb_s
#   if(length(c0)==0){
#     ind = which(df$T==min(df$T))
#     start = c(T = df[ind, t], W = df[ind, w])
#   }else if(length(c0)==nrow(df)){
#     ind = which(df$T==min(df$T))
#     start = c(T = df[ind, t], W = df[ind, w])
#   }else{
#     start = c(T = mean(df[c0, t]), W = mean(df[c0, w]) )
#   }
#
#   # get coordinate of c1_s
#   if(length(c1)==0){
#     c1_s = NULL
#   }else{
#     c1_s = c(T = mean(df[c1, t]), W = mean(df[c1, w]))
#   }
#
#   # get coordinate of c2_s
#   if(length(c2)==0){
#     c2_s = NULL
#   }else{
#     c2_s = c(T = mean(df[c2, t]), W = mean(df[c2, w]))
#   }
#
#   # plot first segment
#   # tb < min(T)
#   # tb > min(T) & c1_s/c2_s = NULL
#   # tb > min(T) & c1_s/c2_s != NULL
#   if(length(c0)==nrow(df)){
#     ind = which(df$T==max(df$T))
#     end = c(T = df[ind, t], W = df[ind, w])
#     seg0 = c(start, end)
#     seg0 = data.frame(t(seg0),color=factor(0, levels=c(0,1,2)))
#     colnames(seg0) = c("x1", "y1", "x2", "y2", "color")
#     p = p + geom_segment(data = seg0, size = 1, alpha = arrowAlpha,
#                          aes(x = x1, y = y1, xend = x2, yend = y2, colour = as.factor(color)),
#                          arrow = arrow(length = unit(arrowSize, "inches"), end = "last"))
#   }else{
#     # seg0 = c(start, tb_s)
#     seg1 = c(start, c1_s)
#     seg2 = c(start, c2_s)
#
#     # seg0 = data.frame(t(seg0),color=factor(0, levels=c(0,1,2)))
#     # colnames(seg0) = c("x1", "y1", "x2", "y2", "color")
#     # if(is.null(c1_s)&is.null(c2_s)){
#     #   p = p + geom_segment(data = seg0, size = 1,
#     #                        aes(x = x1, y = y1, xend = x2, yend = y2, colour = as.factor(color)),
#     #                        arrow = arrow(length = unit(arrowSize, "inches"), end = "last"))
#     # }else{
#     #   p = p + geom_segment(data = seg0, size = 1,
#     #                        aes(x = x1, y = y1, xend = x2, yend = y2, colour = as.factor(color)) )
#     # }
#
#   }
#
#   # plot second and third segment
#   if(is.null(c1_s)&is.null(c2_s)){
#     p = p
#   }else if(!is.null(c1_s) & is.null(c2_s)){
#     seg1 = data.frame(t(seg1), color = factor(1,levels=c(0,1,2)))
#     colnames(seg1) = c("x1", "y1", "x2", "y2", "color")
#     p = p + geom_segment(data = seg1, size = 1, alpha = arrowAlpha,
#                          aes(x = x1, y = y1, xend = x2, yend = y2, colour = as.factor(color)),
#                          arrow = arrow(length = unit(arrowSize, "inches"), end = "last"))
#   }else if(is.null(c1_s) & !is.null(c2_s)){
#     seg2 = data.frame(t(seg2), color = factor(2,levels=c(0,1,2)))
#     colnames(seg2) = c("x1", "y1", "x2", "y2", "color")
#     p = p + geom_segment(data = seg2, size = 1, alpha = arrowAlpha,
#                          aes(x = x1, y = y1, xend = x2, yend = y2, colour = as.factor(color)),
#                          arrow = arrow(length = unit(arrowSize, "inches"), end = "last"))
#   }else{
#     linedf = t(data.frame(seg1, seg2)) %>% as.data.frame
#     linedf$colour = c(1,2)
#     colnames(linedf) = c("x1", "y1", "x2", "y2", "color")
#     linedf = data.frame(linedf)
#
#     p = p + geom_segment(data = linedf, size = 1, alpha = arrowAlpha,
#                          aes(x = x1, y = y1, xend = x2, yend = y2, colour = as.factor(color)),
#                          arrow = arrow(length = unit(arrowSize, "inches"), end = "last") )
#   }
#   return(p)
# }
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' CreateKalmanSdf
#'
#' @param M (Numeric) Number of rows in input data, pc.gene
#' @export
#'
# CreateKalmanSdf <- function(L = NULL,
#                            pse_sdf = NULL,
#                            traj_w = NULL){
#
#   l = L
#   c = paste0("C0_", l)
#
#   track_sdf = data.frame(pse_sdf, traj_w)
#   cdf = track_sdf[, c, drop = FALSE]
#   cdf$C_1 = NA
#   cdf$C_2 = NA
#   cdf$C_1[which(cdf[,c]==0 | cdf[,c]==1)] = 1
#   cdf$C_2[which(cdf[,c]==0 | cdf[,c]==2)] = 2
#   cdf$C_1 = factor(cdf$C_1, levels = c(0,1,2))
#   cdf$C_2 = factor(cdf$C_2, levels = c(0,1,2))
#   cdf = cdf[,c("C_1", "C_2")]
#
#   # colnames(df) = c(paste0("C0_",l,"_",1:2))
#   # track_sdf = cbind(track_sdf, df)
#   # cdf = track_sdf[,c(paste0("C0_",l,"_",1:2))]
#   # col
#
#   tmplist <- lapply(1:ncol(traj_w), function(i){
#
#     w = paste0("W_", i)
#     wdf = track_sdf[, w, drop = FALSE]
#
#     wdf$W_1_1 = cdf$C_1 %>% as.character %>% as.numeric
#     wdf$W_1_2 = cdf$C_2 %>% as.character %>% as.numeric
#     wdf[which(wdf$W_1_1==1), "W_1_1"] = wdf[which(wdf$W_1_1==1), w]
#     wdf[which(wdf$W_1_2==2), "W_1_2"] = wdf[which(wdf$W_1_2==2), w]
#
#     wdf = wdf[,c("W_1_1","W_1_2")]
#     colnames(wdf) = c(paste0("W_",i,"_",1:2) )
#
#     wdf
#   })
#
#   tmp = do.call(cbind, tmplist) %>% as.data.frame
#   colnames(cdf) = c(paste0("C0_",l,"_",1:2))
#   track_sdf = data.frame(names = rownames(track_sdf),
#                          track_sdf, tmp, cdf)
#
#   return(track_sdf)
# }


#' is.julia.setup
#'
#' @description
#' Determine if the Julia environment is connected.
#'
#' @export
#'
is.julia.setup <- function(){
  x = tryCatch({julia_command("a=1");julia_exists("b")}, error = function(e){ "error"})
  if(x[1]=="error"){
    return(FALSE)
  }else{
    return(TRUE)
  }
}

#' SaveMURPDatToJulia
#'
#' @description
#' extract murp data to julia environment
#'
#' @param object MGPfact object
#' @param murp_pc_number the number of principal components of
#' the expression matrix after MURP downsampling
#' @export
#'
SaveMURPDatToJulia <- function(object,
                               murp_pc_number = 3){

  # murp_matrix = ct@MURP$Recommended_K_cl$centers
  # save(murp_matrix, file = "1_murp/murp_matrix.rda")
  Q = murp_pc_number
  pca = object@MURP$centersPCA$x
  murp_matrix_pca = pca[,1:Q,drop=FALSE]
  save(murp_matrix_pca, file = "1_murp/murp_matrix_pca.rda")
}

#' RunningmodMGPpseudoT
#'
#' @description
#' trajectory parameters estimation through MCMC sampling
#'
#' @param object MGPfact object
#' @param julia_home Julia's bin path
#' @param seed random seed, default is 723
#' @param cores number of threads
#'
#' @export
#'
RunningmodMGPpseudoT <- function(object,
                                 julia_home,
                                 seed = 723,
                                 cores = 1){
  if(!is.julia.setup()){
    julia_setup(JULIA_HOME = julia_home)
  }

  cat("/// parameters setting \n")
  # list(getParams(object, "trajectory_number"),
  #      getParams(object, "murp_pc_number"),
  #      getParams(object, "pse_optim_iterations"),
  #      getParams(object, "start_murp"),
  #      getParams(object, "chains_number")) -> x
  # save(x, file = "~/x.rda")

  cmd = list(paste0("trajectory_number=", getParams(object, "trajectory_number")),
             paste0("pc_number=", getParams(object, "murp_pc_number")),
             paste0("iterations=", getParams(object, "pse_optim_iterations")),
             paste0("start_id=[", paste0(getParams(object, "start_murp"),collapse = ";"),"]"),
             paste0("chains_number=", getParams(object, "chains_number")))
  invisible(sapply(cmd, julia_command))
  # cmd = paste0("trajectory_number=", getParams(ct, "trajectory_number"),
  #              "; pc_number=", getParams(ct, "murp_pc_number"),
  #              "; iterations=", getParams(ct, "pse_optim_iterations"),
  #              "; start_id=", getParams(ct, "start_murp"),
  #              "; chains_number=", getParams(ct, "chains_number"))
  # cmd = paste0("trajectory_number=", getParams(ct, "trajectory_number"),
  #              "\n pc_number=", getParams(ct, "murp_pc_number"),
  #              "\n iterations=", getParams(ct, "pse_optim_iterations"),
  #              "\n start_id=", getParams(ct, "start_murp"),
  #              "\n chains_number=", getParams(ct, "chains_number"))
  # julia_command(paste0("trajectory_number=", getParams(ct, "trajectory_number"),
  #                      "; pc_number=", getParams(ct, "murp_pc_number"),
  #                      "; iterations=", getParams(ct, "pse_optim_iterations"),
  #                      "; start_id=", getParams(ct, "start_murp"),
  #                      "; chains_number=", getParams(ct, "chains_number")))

  julia_command("using Distributed")
  julia_command(paste0("addprocs(",cores,")"))

  cat("/// load data \n")

  cmd = '@everywhere using MGPfact, Mamba, RData, Distributions
data_path="1_murp/murp_matrix_pca.rda"
yx = RData.load(data_path)
yx = yx["murp_matrix_pca"][:,1:pc_number]'
#   cmd = '@everywhere using MGPfact, Mamba, RData, Distributions
# data_path="1_murp/murp_matrix_pca.rda"
# yx = RData.load(data_path)
# yx = yx["murp_matrix_pca"][:,1:pc_number]'
  cmd = gsub("\n", ";", cmd)
  julia_command(cmd)

  cat("/// running model \n")

  julia_command("using Random")
  julia_command(paste0("Random.seed!(",seed,")"))

  cat("/// running model \n")

  cmd = 'model = MGPfact.modMGPpseudoT()
SC,inits,scheme = MGPfact.Initialize(yx, pc_number, trajectory_number, start_id, iterations, chains_number)
setinputs!(model, SC)
setinits!(model, inits)
setsamplers!(model, scheme)
@time sim = mcmc(model, SC, inits, iterations, burnin = 0, chains = chains_number)'
  cmd = gsub("\n", ";", cmd)
  julia_command(cmd)

  cat("/// saving \n")

  cmd = 'using JLD2
write(string("2_pseudotime/2.1_julia_result/iter",sim.model.iter,"_bi",sim.model.burnin,".jls"), sim)
@save string("2_pseudotime/2.1_julia_result/inits.jld2") inits'
  cmd = gsub("\n", ";", cmd)
  julia_command(cmd)

  object = ImportPseResult(object = object, init = TRUE)
  return(object)
}

#' ReadPseSim
#'
#' @description
#' read and convert the object optimized by julia into R
#'
#' @param sim julia result
#' @param julia_home Julia's bin path
#' @param init Logical value, whether to import the initial value.
#' it can be viewed in the directory: 2_pseudotime/2.1_julia_result/inits.jld2
#'
#' @export
#'
ReadPseSim <- function(sim = NULL,
                       julia_home = NULL,
                       init = FALSE){
  if(!is.julia.setup()){
    julia_setup(JULIA_HOME = julia_home)
  }

  ## load sim.jls
  cmd = 'using Distributed, RData
nc = 10
addprocs(nc)
@everywhere using Pkg
@everywhere using MGPfact, Mamba, RData, Distributions, KernelFunctions'
  cmd=paste0(cmd,'\nsim1=read("2_pseudotime/2.1_julia_result/',sim,'",ModelChains)')
  cmd = gsub("\n", ";", cmd)
  julia_command(cmd)

  ## load inits
  if(init){
    julia_command("using JLD2")
    julia_command('@load "2_pseudotime/2.1_julia_result/inits.jld2"')
  }
}

#' ImportPseSimToR
#'
#' @description
#' Import the pseudotime object optimized by julia into R
#'
#' @param object MGPfact object
#' @param init
#' Logical value, whether to import the initial value.
#' it can be viewed in the directory: 2_pseudotime/2.1_julia_result/inits.jld2
#' @param julia_home julia bin path
#' @export
#'
ImportPseResult <- function(object,
                            init = FALSE,
                            julia_home = FALSE){

  if(!is.julia.setup()){
    julia_setup(JULIA_HOME = julia_home)
  }

  if(init){
    cmd = 'chains = sim.chains
names2 = sim.names
values1 = sim.value
burnin = sim.model.burnin
iter =  sim.model.iter
hasinits = sim.model.hasinits
@rput chains burnin iter inits names2 values1'
    cmd = gsub("\n", ";", cmd)
    julia_command(cmd)
    OptimResult <- CreatePseOptimResultObject(sample_value = values1,
                                              dimnames = list(paste0("iter",1:dim(values1)[1]),
                                                              names=names2,
                                                              paste0("chain",1:dim(values1)[3])),
                                              burnin = burnin,
                                              chains = chains,
                                              iter = iter,
                                              inits = as.list(inits),
                                              hasinits = TRUE)
  }else{
    cmd = 'chains = sim.chains
names2 = sim.names
values1 = sim.value
burnin = sim.model.burnin
iter =  sim.model.iter
hasinits = sim.model.hasinits
@rput chains burnin iter names2 values1'
    cmd = gsub("\n", ";", cmd)
    julia_command(cmd)
    OptimResult <- CreatePseOptimResultObject(sample_value = values1,
                                              dimnames = list(paste0("iter",1:dim(values1)[1]),
                                                              names=names2,
                                                              paste0("chain",1:dim(values1)[3])),
                                              burnin = burnin,
                                              chains = chains,
                                              iter = iter,
                                              # inits = as.list(inits),
                                              hasinits = FALSE)
  }

  object@OptimResult$pse = OptimResult
  return(object)
}

#' GetSigmaFromJulia
#'
#' @description get the covariance matrix
#'
#' @param object MGPfact object
#' @param init initial values
#' @param load logical value, whether to load julia's environment
#' @param julia_home Julia's bin path
#'
#' @export
#'
GetSigmaFromJulia <- function(object,
                              load = FALSE,
                              julia_home = NULL,
                              init = T){

  if(!is.julia.setup()){
    julia_setup(JULIA_HOME = julia_home)
    julia_command("using MGPfact, Mamba, RData, Distributions")
  }

  ## 1. load
  if(load){
    ReadPseSim(sim = paste0("iter",getParams(object, "pse_optim_iterations"),"_bi0.jls"),
               julia_home = julia_home,
               init = init)
#    cmd = 'using Distributed, RData
# nc = 10
# addprocs(nc)
# @everywhere using Pkg
# @everywhere using MGPfact, Mamba, RData, Distributions, KernelFunctions'
# cmd = gsub("\n", ";", cmd)
# julia_command(cmd)
  }

  ## 2. prepare in julia
  sdf = GetMURPInfo(object=object)
  cat(dim(sdf))
  julia_assign("sdf",sdf)
  # julia_command("using RCall; @rget sdf")

  julia_command("using KernelFunctions")
  cmd = 'P = size(sdf)[1]
L = size(filter(x->occursin("Tb", string(x)), names(sdf)))[1]
Q = size(filter(x->occursin("PC", string(x)), names(sdf)))[1]
Tb = Vector(sdf[1,filter(x->occursin("Tb", string(x)), names(sdf))])
lambda1 = Vector(sdf[1,filter(x->occursin("lambda1", string(x)), names(sdf))])
lambda2 = Vector(sdf[1,filter(x->occursin("lambda2", string(x)), names(sdf))])
lambda3 = Vector(sdf[1,filter(x->occursin("lambda3", string(x)), names(sdf))])
lambda = hcat(lambda1,lambda2,lambda3)
T = sdf.T
X = sdf.X
m_t = sdf.mt[1]
s2_t = sdf.s2t[1]
s2_x = sdf.s2x[1]
rho = sdf.rho[1]'
  cmd = gsub("\n", ";", cmd)
  julia_command(cmd)

  julia_command('tmp = sdf[:,filter(x->occursin(r"^C_", string(x)), names(sdf))]')
  julia_command("C = convert(Array{Int},Matrix(tmp)')")

  ## 3. get cov matrix
  cmd = 'sigma = 1.0
delta = 1.0
P = Int(P)
L = Int(L)
S_list = [( tb = Tb[i];
            cl = C[i,:];
            le = lambda[i,:];
            k = MGPfact.track.bifurcateKernel(tb, le, sigma, 0.0, delta);
            x = [MGPfact.track.Trejactory(T[j], cl[j], X[j]) for j in 1:P];
            ss = KernelFunctions.kernelmatrix(k, x);
            Matrix(Hermitian(ss .+ 1E-6 .* Matrix(I, P, P)))) for i in 1:L ]
S = sum(S_list)'
  cmd = gsub("\n", ";", cmd)
  cmd = gsub(";;", ";", cmd)
  julia_command(cmd)

  ## 4. get 100 points sigma
  cat("2 \n")
  cmd = 'P_ = 100
T1_ = collect(1:(P_) )/(P_)
T2_ = collect(1:(P_) )/(P_)
T_ = vec(vcat(T1_,T2_))
C1_ = fill(1, L, P_)
C2_ = fill(2, L, P_)
C_ = hcat(C1_,C2_)
sort_index = sortperm(T_)
T_ = T_[sort_index]
C_ = C_[:,sort_index]
P_ = P_*2
L = Int(L)
X_ = rand(truncated(Normal(0, sqrt(s2_x[1])), 0.0, 1.0), P_)
S_list_ = [( tb = Tb[i];
             le = lambda[i,:];
             k = MGPfact.track.bifurcateKernel(tb, le, sigma, 0.0, delta);
             x = [MGPfact.track.Trejactory(vcat(T,T_)[j], hcat(C,C_)[i,j], vcat(X, X_)[j]) for j in 1:(P+P_)];
             ss = KernelFunctions.kernelmatrix(k, x);
             Matrix(Hermitian(ss .+ 1E-6 .* Matrix(I, P+P_, P+P_)))) for i in 1:L ]
S_ = sum(S_list_)'
  cmd = gsub("\n", ";", cmd)
  cmd = gsub(";;", ";", cmd)
  julia_command(cmd)

  ## 5. save in julia
  cat("3 \n")
  cmd = 'using JLD2
  @save "2_pseudotime/2.4_cov_matrix/s_.jld2" S_ S_list_ C_ T_ X_ P_ rho m_t s2_t s2_x Tb lambda L Q P
  @save "2_pseudotime/2.4_cov_matrix/s.jld2" S S_list C T X rho m_t s2_t s2_x Tb lambda L Q P
  @rput S_ S_list_ C_ T_ X_ P_ rho m_t s2_t s2_x Tb lambda L Q P
  @rput S S_list C T X rho m_t s2_t s2_x Tb lambda L Q P'
  cmd = gsub("\n", ";", cmd)
  julia_command(cmd)

  ## 6. save in R
  S_list_ = as.list(S_list_)
  S_list = as.list(S_list)
  save(S_,S_list_,C_,T_,X_,P_,rho,m_t,s2_t,s2_x,Tb,lambda,L,Q,P, file = "2_pseudotime/2.4_cov_matrix/s_.rda")
  save(S,S_list,C,T,X,rho,m_t,s2_t,s2_x,Tb,lambda,L,Q,P, file = "2_pseudotime/2.4_cov_matrix/s.rda")

  object@GPR$reg_cov = list(cov = S_, cov_l = S_list_,
                            C = C_, T = T_, X = X_, P = P_,
                            rho = rho, m_t = m_t, s2_t = s2_t, s2_x = s2_x, Tb = Tb, lambda = lambda)
  object@GPR$murp_cov = list(cov = S, cov_l = S_list,
                             C = C, T = T, X = X,
                             rho = rho, m_t = m_t, s2_t = s2_t, s2_x = s2_x, Tb = Tb, lambda = lambda)
  return(object)
}

#' ComputeLogpdf
#'
#' @description
#' compute probability density values for different parameters
#' @param object MGPfact object
#' @param import_sim logical value, whether import result from julia environment
#' @param julia_home Julia's bin path
#' @export
ComputeLogpdf <- function(object,
                          import_sim = FALSE,
                          julia_home=NULL ){

  if(import_sim){
    ReadPseSim(sim =  paste0("iter",getParams(object, "pse_optim_iterations"),"_bi0.jls"),
               julia_home = julia_home, init = TRUE)
  }

  cmd = 'using JLD2
chain = size(sim1)[3]
iterations = size(sim1.value)[1]
burnin = size(sim1)[1] - size(sim1.value)[1]
sigma = sim1.model.nodes[:sigma]
delta = sim1.model.nodes[:delta]
P = sim1.model.nodes[:P]
L = sim1.model.nodes[:L]
Q = sim1.model.nodes[:Q]'
  cmd = gsub("\n", ";", cmd)
  julia_command(cmd)

  cat("all \n")
  cmd='value = logpdf(sim1)
l_all = zeros(iterations, chain)
for i in 1:iterations l_all[i, :] .= value.value[i,1,:] end'
  cmd = gsub("\n", ";", cmd)
  julia_command(cmd)

  cat("rho & pi \n")
  cmd = 'rho = hcat([sim1[burnin:(burnin+iterations), "rho", :].value ]...)
a = hcat([[sim1[burnin:(burnin+iterations), "a[$i, $j]", :].value for i in 1:L] for j in 1:P]...)
p = hcat([[sim1[burnin:(burnin+iterations), "p[$i, $j]", :].value for i in 1:L] for j in 1:P]...)
C = hcat([[sim1[burnin:(burnin+iterations), "C[$i, $j]", :].value for i in 1:L] for j in 1:P]...)
l_rho = zeros(iterations, chain)
for i in 1:iterations l_rho[i, :] .= logpdf.(truncated(Beta(1.0, 1.0), 1E-8, 0.99999999), vcat(rho[i,:,1:chain]...)) end
l_p = zeros(iterations, chain)
for i in 1:iterations l_p[i, :] .= [sum(hcat([vcat([logpdf(Beta(a[l,j][i,1,ch]*(1-rho[i, 1, ch]), a[l,j][i,1,ch]*rho[i, 1, ch]), p[l,j][i,1,ch]) for l in 1:L]...) for j in 1:P]...)) for ch in 1:chain] end'
  cmd = gsub("\n", ";", cmd)
  julia_command(cmd)

  cat("Tb \n")
  cmd = 'Tb = hcat([sim1[burnin:(burnin+iterations), "Tb[$i]", :].value for i in 1:L ]...)
l_tb = zeros(iterations, chain)
for i in 1:iterations l_tb[i, :] .= [sum([logpdf(truncated(Gamma(l, 1/L), 0.0, 1.0),Tb[i, l, ch]) for l in 1:L]) for ch in 1:chain] end'
  cmd = gsub("\n", ";", cmd)
  julia_command(cmd)

  cat("lambda \n")
  cmd = 'lambda = hcat([[sim1[burnin:(burnin+iterations), "lambda[$i, $j]", :].value for i in 1:L] for j in 1:3]...)
lambda1 = hcat(lambda[:, 1]...)
lambda2 = hcat(lambda[:, 2]...)
lambda3 = hcat(lambda[:, 3]...)
l_l1 = zeros(iterations, chain)
for i in 1:iterations l_l1[i, :] .= [sum(logpdf(sim1.model.nodes[:lambda1], lambda1[i,:,ch])) for ch in 1:chain] end
l_l2 = zeros(iterations, chain)
for i in 1:iterations l_l2[i, :] .= [sum(logpdf(sim1.model.nodes[:lambda2], lambda2[i,:,ch])) for ch in 1:chain] end
l_l3 = zeros(iterations, chain)
for i in 1:iterations l_l3[i, :] .= [sum(logpdf(sim1.model.nodes[:lambda3], lambda3[i,:,ch])) for ch in 1:chain] end'
  cmd = gsub("\n", ";", cmd)
  julia_command(cmd)

  cat("T \n")
  cmd = 'm_t = hcat([sim1[burnin:(burnin+iterations), "m_t", :].value ]...)
l_mt = zeros(iterations, chain)
for i in 1:iterations l_mt[i, :] .= logpdf.(sim1.model.nodes[:m_t], vcat(m_t[i,:,1:chain]...)) end
s2_t = hcat([sim1[burnin:(burnin+iterations), "s2_t", :].value ]...)
l_s2t = zeros(iterations, chain)
for i in 1:iterations l_s2t[i, :] .= logpdf.(sim1.model.nodes[:s2_t], vcat(s2_t[i,:,1:chain]...)) end'
  cmd = gsub("\n", ";", cmd)
  julia_command(cmd)

  cmd = 'T = hcat([sim1[burnin:(burnin+iterations), "T[$i]", :].value for i in 1:P]...)
l_t = zeros(iterations, chain)
rt = sim1.model.nodes[:rt]'
  cmd = gsub("\n", ";", cmd)
  julia_command(cmd)

  cmd = 'for i in 1:iterations
    tmp = zeros(P,chain)
    for ch in 1:chain
      for i2 in 1:P
        if i2 in rt
          tmp[i2,ch]=logpdf(truncated(Normal(0, 0.01), 0.0, 1.0), T[i,i2,ch])
        else
          m_t2 = m_t[i,1,ch]
          s2_t2 = s2_t[i,1,ch]
              tmp[i2,ch]=logpdf(truncated(Normal(m_t2, sqrt(s2_t2)), 0.0, 1.0), T[i,i2,ch])
          end
      end
  end
  l_t[i,:] .=  vcat(sum(tmp, dims=1)...)
end'
  julia_command(cmd)

  cat("X \n")
  cmd = 's2_x = hcat([sim1[burnin:(burnin+iterations), "s2_x", :].value ]...)
l_s2x = zeros(iterations, chain)
for i in 1:iterations l_s2x[i, :] .= logpdf.(sim1.model.nodes[:s2_x], vcat(s2_x[i,:,1:chain]...)) end
X = hcat([sim1[burnin:(burnin+iterations), "X[$i]", :].value for i in 1:P]...)
l_x = zeros(iterations, chain)
for i in 1:iterations l_x[i, :] .= [sum(logpdf(truncated(Normal(0, sqrt(s2_x[i,1,ch])), -1.0, 1.0), X[i, :, ch])) for ch in 1:chain] end'
  cmd = gsub("\n", ";", cmd)
  julia_command(cmd)

  # cat("9 \n")
  # cmd='S = [[
  #        Matrix(Hermitian(
  #         sum([(
  #             tb = Tb[i, l, ch];
  #             cl = hcat(C[l,:]...)[i,:,ch];
  #             le = [lambda1[i,l,ch], lambda2[i,l,ch], lambda3[i,l,ch]];
  #             k = MGPfact.pseudot.bifurcateKernel(tb, le, sigma, 0.0, delta);
  #             gp = [MGPfact.pseudot.Trejactory(T[i, p, ch], cl[p], X[i, p, ch]) for p in 1:P];
  #             KernelFunctions.kernelmatrix(k, gp)
  #             )
  #         for l in 1:L]) .+ 1E-6 .* Matrix(I, P, P) ))
  # for i in 1:5] for ch in 1:10]'
  #     julia_command(cmd)

  cmd='@rput l_all l_rho l_tb l_mt l_s2t l_s2x l_l1 l_l2 l_l3 l_t l_x l_p chain iterations
@save "2_pseudotime/2.1_julia_result/logpdf.jld2" l_all l_rho l_tb l_mt l_s2t l_s2x l_l1 l_l2 l_l3 l_t l_x l_p chain iterations'
  cmd = gsub("\n", ";", cmd)
  julia_command(cmd)

  logpdf = list(l_all = l_all,
                l_p = l_p,
                l_rho = l_rho,
                l_tb = l_tb,
                l_l1 = l_l1, l_l2 = l_l2, l_l3 = l_l3,
                l_x = l_x, l_s2x = l_s2x,
                l_t = l_t, l_mt = l_mt, l_s2t = l_s2t)

  mamba=object@OptimResult$pse
  mamba@logpdf$logpdf = logpdf

  chain = length(mamba@chains)
  iterations = mamba@iter
  logpdf_cor = lapply(c("spearman","pearson"), function(m){

    lapply(names(logpdf), function(x){
      lp = logpdf[[x]]
      apply(lp, 2, function(y) cor(y, 1:iterations, method = m) )
    }) -> tmp
    cor_df = do.call(rbind,tmp)
    dimnames(cor_df) = list(names(logpdf),paste0("chain",1:chain))
    return(cor_df)
  })
  names(logpdf_cor) = c("spearman","pearson")

  corr_spearman = logpdf_cor[[1]]
  corr_pearson = logpdf_cor[[1]]

  write.csv(corr_spearman, file = "2_pseudotime/2.1_julia_result/logpdf_iter_cor_spearman.csv")
  write.csv(corr_pearson, file = "2_pseudotime/2.1_julia_result/logpdf_iter_cor_pearson.csv")

  save(corr_spearman, corr_pearson, file = "2_pseudotime/2.1_julia_result/logpdf_iter_cor.rda")
  mamba@logpdf$cor = list(pearson = corr_pearson,
                          spearman = corr_spearman)
  object@OptimResult$pse=mamba
  return(object)
}

# -------------
#' #' # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' #' SaveMURPDatToJulia
#' #'
#' #' @description
#' #' extract murp data to julia environment
#' #' @param object mgpfact object
#' #' @param murp_pc_number the number of principal components of the expression matrix after MURP downsampling
#' #' @export
#' #'
#' SavePCAMURPDatToJulia <- function(object,
#'                                   murp_pc_number = 5){
#'
#'   # murp_matrix = ct@MURP$Recommended_K_cl$centers
#'   # save(murp_matrix, file = "1_murp/murp_matrix.rda")
#'   Q = murp_pc_number
#'   murp_matrix_pca = ct@MURP$Recommended_K_cl$centers[,1:Q]
#'   save(murp_matrix_pca, file = "1_murp/murp_matrix_pca.rda")
#' }

#' # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' #' ImportLogpdf
#' #'
#' #' @description
#' #'
#' ImportLogpdf <- function(){
#'   julia_setup(JULIA_HOME = "/public/home/renjun/tool/julia-1.6.6/bin")
#'
#'   julia_command("using JLD2, RCall")
#'   julia_command('@load "2_pseudotime/2.1_julia_result/logpdf.jld2"')
#'   julia_command('@rput l l_rho l_tb l_mt l_s2t l_s2x l_l1 l_l2 l_l3 l_t l_x l_p chain iterations')
#' }

#' # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' #' ImportTrackOptimToR
#' #'
#' #' @description
#' #'
#' #' @param object import tracking object from julia to R
#' #'
#' ImportTrackOptimToR <- function(object){
#'
#'   P = getParams(object, "murp_number")
#'   L = getParams(object, "trajectory_number")
#'   Q = getParams(object, "murp_pc_number")
#'
#'   ## load data from julia
#'   julia_setup(JULIA_HOME = "/public/home/renjun/tool/julia-1.6.6/bin")
#'
#'   cmd1 = 'sample_value = map(x -> r.trace[x].metadata["x"], 1:length(r.trace))'
#'   cmd2 = "sample_value = hcat(sample_value...)'"
#'   cmd3 = 'init_value = Optim.initial_state(r)
#' iterations = length(r.trace)
#' @rput sample_value init_value iterations'
#'   julia_command(cmd1)
#'   julia_command(cmd2)
#'   julia_command(gsub("\n", ";", cmd3))
#'
#'   ## set iteration in settings
#'   object = assignSettings(object, "track_optim_iterations", iterations)
#'
#'   ## create object in R
#'   names <- lapply(1:L, function(i) paste0("W[", i, ", ",1:P,"]")) %>% unlist
#'   burnin = 0
#'   chains = 1
#'   hasinits = TRUE
#'   inits = list(W = init_value[1:(length(init_value)-1)],
#'                s2 = tail(init_value,1))
#'   values1 = array(sample_value, dim = c(iterations, L*P+1, 1))
#'   OptimResult <- CreateTrackOptimResultObject(sample_value = values1,
#'                                               dimnames = list(paste0("iter", 1:iterations),
#'                                                               names = c(names,"s2"),
#'                                                               paste0("chain",1:dim(values1)[3])),
#'                                               burnin = burnin,
#'                                               chains = chains,
#'                                               iter = iterations,
#'                                               hasinits = hasinits,
#'                                               inits = inits)
#'   object@OptimResult$track = OptimResult
#'   return(object)
#' }

#' # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' #' ReadTrackSim
#' #'
#' #' @description
#' #' Import the object optimized by julia into R
#' #'
#' #' @param sim julia result
#' #' @param init
#' #' Logical value, whether to import the initial value.
#' #' it can be viewed in the directory: 2_pseudotime/2.1_julia_result/inits.jld2
#' #'
#' ReadTrackSim <- function(sim = NULL,
#'                          init = FALSE){
#'   julia_setup(JULIA_HOME = "/public/home/renjun/tool/julia-1.6.6/bin")
#'
#'   ## load sim.jls
#'   cmd = 'using Distributed, RData
#' nc = 10
#' addprocs(nc)
#' @everywhere using Pkg
#' @everywhere Pkg.activate(".")
#' @everywhere using StatsBase, Distributions, KernelFunctions, RandomMatrices
#' @everywhere using LinearAlgebra, Mamba, DataFrames, CSV
#' @everywhere include("/share/data6/tmp/renjun/CellTrekResult/CellTrek/merge_MURP_PCA/0_0_celltrek.jl")
#' @everywhere using .pseudot
#' @everywhere using .track'
#'   cmd=paste0(cmd,'\nsim1=read("3_tracking/3.1_optim_result/',sim,'",ModelChains)')
#'   cmd = gsub("\n", ";", cmd)
#'   julia_command(cmd)
#'
#'   ## load inits
#'   if(init){
#'     julia_command("using JLD2")
#'     julia_command('@load "3_tracking/3.1_optim_result/inits.jld2"')
#'   }
#' }

#' # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' #' ImportTrackSimToR
#' #'
#' #' @description
#' #'
#' ImportTrackSimToR <- function(object, init = FALSE){
#'
#'   julia_setup(JULIA_HOME = "/public/home/renjun/tool/julia-1.6.6/bin")
#'
#'   if(init){
#'     cmd = 'sim = sim1
#' chains = sim.chains
#' names2 = sim.names
#' values1 = sim.value
#' burnin = sim.model.burnin
#' iter =  sim.model.iter
#' hasinits = sim.model.hasinits
#' @rput chains burnin iter inits names2 values1'
#'     cmd = gsub("\n", ";", cmd)
#'     julia_command(cmd)
#'     OptimResult <- CreateTrackOptimResultObject(sample_value = values1,
#'                                                 dimnames = list(paste0("iter",1:dim(values1)[1]),
#'                                                                 names=names2,
#'                                                                 paste0("chain",1:dim(values1)[3])),
#'                                                 burnin = burnin,
#'                                                 chains = chains,
#'                                                 iter = iter,
#'                                                 inits = as.list(inits),
#'                                                 hasinits = TRUE)
#'
#'   }else{
#'     cmd = 'sim = sim1
#' chains = sim.chains
#' names2 = sim.names
#' values1 = sim.value
#' burnin = sim.model.burnin
#' iter =  sim.model.iter
#' hasinits = sim.model.hasinits
#' @rput chains burnin iter names2 values1'
#'     cmd = gsub("\n", ";", cmd)
#'     julia_command(cmd)
#'     OptimResult <- CreateTrackOptimResultObject(sample_value = values1,
#'                                                 dimnames = list(paste0("iter",1:dim(values1)[1]),
#'                                                                 names=names2,
#'                                                                 paste0("chain",1:dim(values1)[3])),
#'                                                 burnin = burnin,
#'                                                 chains = chains,
#'                                                 iter = iter,
#'                                                 # inits = as.list(inits),
#'                                                 hasinits = FALSE)
#'   }
#'
#'   save(OptimeResult, file = "3_tracking/3.1_optim_result/OptimeResult.rda")
#'   object@OptimResult$track = OptimResult
#'   return(object)
#' }

#' # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' #' ReadTrackSim
#' #'
#' #' @description
#' #'
#' ReadTrackOptim <- function(sim = NULL){
#'   julia_setup(JULIA_HOME = "/public/home/renjun/tool/julia-1.6.6/bin")
#'
#'   ## load sim.jls
#'   cmd = 'using Optim, JLD2'
#'   cmd = paste0(cmd,'\n@load "3_tracking/3.1_optim_result/',sim,'"')
#'   cmd = gsub("\n", ";", cmd)
#'   julia_command(cmd)
#' }

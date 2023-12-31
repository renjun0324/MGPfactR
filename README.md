
<!-- <br><img src="test/logo.png" align="right" width="300"/> -->
# MGPfactR

A model-based, unsupervised manifold learning method that factors complex cellular trajectories into interpretable bifurcating Gaussian processes of transcription. The complete functionality of MGPfact is accessible in [MGPfactR](https://github.com/renjun0324/MGPfactR), enabling the discovery of specific biological determinants of cell fate.

## installation

### (1) Install julia environment

```shell
JULIA_VERSION=1.6.6
sudo wget https://julialang-s3.julialang.org/bin/linux/x64/$(echo $JULIA_VERSION | cut -d. -f 1-2)/julia-$JULIA_VERSION-linux-x86_64.tar.gz \
    && tar -xvzf julia-$JULIA_VERSION-linux-x86_64.tar.gz -C ./ 
echo 'export PATH=$PATH:'"$(pwd)"/julia-$JULIA_VERSION/bin >> ~/.bashrc
rm julia-$JULIA_VERSION-linux-x86_64.tar.gz
```

### (2) Install MGPfact.jl and its dependencies

```shell
julia -e 'import Pkg; Pkg.add(url="https://github.com/renjun0324/MGPfact.jl")'
julia -e 'import Pkg; Pkg.add(["Mamba", "RData", "JLD2", "Distributions", "KernelFunctions"])'
julia -e 'import Pkg; Pkg.add(Pkg.PackageSpec(name="RCall", version="0.13.15"))'
julia -e 'import Pkg; Pkg.add(Pkg.PackageSpec(name="Suppressor", version="0.2.6"))'
```

```julia
# Test whether MGPfact.jl can be loaded
using MGPfact
using Mamba, RData, JLD2
```

### (3) Install MGPfactR packages

```shell
R -e 'devtools::install_version("JuliaCall","0.16")'
R -e 'devtools::install_github("renjun0324/MURP@v0.6.5")'
R -e 'devtools::install_cran(c("dplyr", "purrr", "stringr","JuliaCall", "pbmcapply", "doParallel", "reshape", "reshape2", "igraph", "graphlayouts","oaqc","parallelDist"))'
```

```r
# Test whether MGPfactR can be loaded
library(MGPfactR)

# Test if the R environment can be linked with the Julia environment
library(JuliaCall)
julia_home = gsub("/julia$","",system("which julia", intern = T))
julia_setup(JULIA_HOME=julia_home)
```

## quick start
For the detailed usage process of MGPfact, please click [here]()
```r
data(fibroblast_reprogramming_treutlein)
data = fibroblast_reprogramming_treutlein
counts = data$counts
cell_info = data$cell_info
rownames(cell_info) = cell_info$cell_id

expression = LogNormalize(t(counts)) %>% t

# create object
ct <- CreateMGPfactObject(data_matrix = expression, MetaData = cell_info)

# downsampling
ct = MURPDownsampling(ct, omega = 0.9, iter = 10, seed = 723, fast = T, cores = 1,
                      pca.center = FALSE, pca.scale = FALSE, plot = T, max_murp = 20)
ct = GetMURPMapLabel(ct, labels = "time_point")

# initialize parameters
SaveMURPDatToJulia(ct, murp_pc_number = 3)
ct = SetSettings(ct, murp_pc_number = 3, trajectory_number = 3, pse_optim_iterations = 100, start_murp = 999)

# forcast pseudotime
ct = RunningmodMGPpseudoT(ct, julia_home = julia_home, cores = 1)

# trahectory construction
ct <- GetIterCor(object = ct, iteration_list = list(c(1,getParams(ct,"pse_optim_iterations"))))
param_meanvalue <- GetAllParamMeanChain(object = ct, aspect = "pse", 
                                        iter_range = c(getParams(ct,"pse_optim_iterations"), 
                                                       getParams(ct,"pse_optim_iterations")))
ct <- GetPredT(object = ct, chains = 1:getParams(ct, "chains_number"), adjust = F, filter_chain = F, mean_th = 0)
ct <- GetPseSdf(ct, unified_direction = FALSE, rm_adjust_chain = TRUE, param_meanvalue = param_meanvalue)
ct <- GetBinTree(object = ct)
ct <- GetTbTree(object = ct)
ct <- GetTbTreeAllpoint(object = ct, save = T, labels = getParams(ct,"label"))
PlotPieBinLabel(ct, labels = getParams(ct,"label"))
PlotPieTbLabel(ct, labels = getParams(ct,"label"))
PlotPieConsensusMainLabel(ct, labels = getParams(ct,"label"))
PlotPieConsensusAllLabel(ct, labels = getParams(ct,"label"),size = 0.005)

```








<!-- <br><img src="test/logo.png" align="right" width="300"/> -->
# MGPfactR

A model-based, unsupervised manifold learning method that factors complex cellular trajectories into interpretable bifurcating Gaussian processes of transcription. The complete functionality of MGPfact is accessible in [MGPfactR](https://github.com/renjun0324/MGPfactR), enabling the discovery of specific biological determinants of cell fate.

## installation

### (1) Install julia environment
```shell
ENV JULIA_VERSION 1.6.6
wget https://julialang-s3.julialang.org/bin/linux/x64/$(echo $JULIA_VERSION | cut -d. -f 1-2)/julia-$JULIA_VERSION-linux-x86_64.tar.gz \
    && tar -xvzf julia-$JULIA_VERSION-linux-x86_64.tar.gz -C /usr/local/ \
    && ln -s /usr/local/julia-$JULIA_VERSION/bin/julia /usr/local/bin/julia \
    && rm julia-$JULIA_VERSION-linux-x86_64.tar.gz
```

### (2) Install MGPfact.jl and its dependencies
```shell
julia -e 'import Pkg; Pkg.add(url="https://github.com/renjun0324/MGPfact.jl")'
julia -e 'import Pkg; Pkg.add(["Mamba", "RData", "JLD2", "Distributions"])'
julia -e 'import Pkg; Pkg.add(Pkg.PackageSpec(name="RCall", version="0.13.15"))'
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
julia_home = gsub("/julia","",system("which julia", intern = T))
julia_setup(JULIA_HOME=julia_home)
```

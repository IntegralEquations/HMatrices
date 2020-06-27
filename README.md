# HMatrices.jl

*A package for assembling and doing algebra with hierarchical matrices with a focus on boundary integral equations* 

![Lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)
![CI](https://github.com/IntegralEquations/HMatrices/workflows/CI/badge.svg?branch=master)

## Installation
Install from the Pkg REPL:
```
pkg> add https://github.com/IntegralEquations/HMatrices
```
## Basic usage

A simple examples which uses the default parameters would be as follows:
```julia
    using HMatrices, Clusters, Plots, LinearAlgebra
    # create some random points
    pts         = rand(2,2000)
    # make a cluster tree
    clustertree = ClusterTree(pts)
    # then a block cluster tree
    blocktree   = BlockTree(clustertree,clustertree)
    # create your matrix
    f(x,y)      = exp(-norm(x-y)^2)
    M           = [f(x,y) for x in clustertree.data, y in clustertree.data]
    # compress it
    H           = HMatrix(M,blocktree)
```
Many of the steps above accept keyword arguments or functors for modifying their default behavior.

Often one cannot assemble the full matrix. In this case the `LazyMatrix` type is useful:
```julia
    L = LazyMatrix(f,clustertree.data,clustertree.data)
```
This is just like the matrix we build `M`, but it computes the entries *on demand* and does not store them.

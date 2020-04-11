# HMatrices.jl

![Lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)
![CI](https://github.com/IntegralEquations/HMatrices/workflows/CI/badge.svg?branch=master)
[![codecov.io](http://codecov.io/github/IntegralEquations/HMatrices.jl/coverage.svg?branch=master)](http://codecov.io/github/IntegralEquations/HMatrices.jl?branch=master)

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

## Integration with `BoundaryIntegralEquations`

Below is an example of how `HMatrix` is used together with e.g. `BoundaryIntegralEquations` to solve Helmholtz equations in 3d over a sphere:
```julia
using HMatrices, Clusters, Plots, LinearAlgebra, ParametricSurfaces, BoundaryIntegralEquations, Kernels, IterativeSolvers
geo    = Sphere()
refine!(geo)
qsize  = (10,10) #number of quadrature nodes in each direction
quad   = TensorQuadrature(qsize,geo)
pts    = quad.nodes |> vec
# build a cluster tree
clt    = ClusterTree(pts)
# then a block cluster tree
blocktree   = BlockTree(clt,clt)
# use the cluster tree numbering for the quadrature
permute!(quad,clt.perm)
# define the operators
pde = Helmholtz(k=1)
dim = 3
G      = SingleLayerKernel{dim}(pde)
dGdn   = DoubleLayerKernel{dim}(pde)
BoundaryIntegralEquations.kerneltype(dGdn) = DoubleLayer()
σ   = Density{ComplexF64}(quad)
𝒟  = IntegralPotential(dGdn,quad)x
D   = IntegralOperator(dGdn,quad,quad)
# chose a representation 
u   = 𝒟[σ]
# exact solution for reference
uₑ(x) =  G((0.1,0.2,-0.1),x)
# trace of exact solution
γ₀u = γ₀(uₑ,quad)
# solve for density σ
gmres!(σ,I/2 + D,γ₀u;verbose=true,maxiter=100,restart=1000)
# check the result on some point
xobs = (10,10,10)
rel_error = norm(u(xobs) - uₑ(xobs))/norm(uₑ(xobs))
println("relative error without compression: $(100*rel_error) %")
# now compress the matrix and do the same
H         = HMatrix(D,blocktree)
σ         = Density{ComplexF64}(quad)
gmres!(σ,I/2 + H,γ₀u;verbose=true,maxiter=100,restart=1000)
rel_error = norm(u(xobs) - uₑ(xobs))/norm(uₑ(xobs))
println("relative error with compression: $(100*rel_error) %")
```
 




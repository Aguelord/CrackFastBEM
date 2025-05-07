module CrackFastBEM

@info "Loading CrackFastBEM..."
@info "Loading dependencies packages..."

using Inti
using Gmsh
using HMatrices
using LinearAlgebra
using LinearMaps
using StaticArrays
using SparseArrays
using IterativeSolvers
using Optim
using Statistics
using NearestNeighbors
using Meshes
using Plots

# utils
include("utils.jl")

# Mesh API
include("mesh_api.jl")

# pre-processing
include("pre_processor.jl")
include("weight_functions.jl")

# solver
include("solver.jl")

# post-processing, visualization
include("graphic.jl")

# api
include("api.jl")

@info "CrackFastBEM successfully loaded!"

end

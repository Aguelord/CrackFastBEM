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
using GLMakie

# utils
include("utils.jl")
include("simple_crack_problems.jl")

# External API
include("mesh_api.jl")
include("inti_api.jl")

# pre-processing
include("pre_processor.jl")
include("weight_functions.jl")

# solver
include("solver.jl")

# post-processing, visualization
include("graphic.jl")

# Internal API
include("api.jl")

@info "CrackFastBEM successfully loaded!"

end

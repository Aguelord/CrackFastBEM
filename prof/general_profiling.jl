using Profile
using Pkg

Pkg.activate(".")


using FlameGraphs
using CrackBEM
using Inti
using Gmsh
using StaticArrays
using LinearAlgebra
using Statistics

meshsize = 0.5
power = 3

radius_x = 2.0
radius_y = 1.0

gmsh.initialize(String[], true)
gmsh.option.setNumber("Mesh.MeshSizeMin", meshsize)
gmsh.option.setNumber("Mesh.MeshSizeMax", meshsize)
# gmsh.option.setNumber("Mesh.Algorithm", 11)
gmsh.model.occ.addDisk(0.0, 0.0, 0.0, radius_x, radius_y)
gmsh.model.occ.synchronize()
gmsh.model.mesh.generate(2)
gmsh.model.mesh.recombine()
gmsh.model.mesh.setOrder(2)
msh = Inti.import_mesh(; dim = 3)
# gmsh.fltk.run()
gmsh.finalize()

Γ = Inti.Domain(e -> Inti.geometric_dimension(e) == 2, msh)

crack_mesh = msh[Γ]

μ = 1
ν = 0.15

material = CrackBEM.Material(
	"A-dimensional material",
	2 * μ * (1 + ν),
	ν,
	0.0,
	1.0,
)

boundary_mesh = CrackBEM.empty_mesh()
dirichlet_boundary_condition = Dict{Int, SVector{3, Float64}}()
neumann_boundary_condition = Dict{Int, SVector{3, Float64}}()
σ_inf = 1.0
traction_crack_condition = Dict(node_id => -SVector(0.0, 0.0, σ_inf) for (node_id, node) in enumerate(crack_mesh.nodes))
SIFs = Dict{Int, SVector{3, Float64}}()

p = CrackBEM.CrackProblem(
	"Test problem",
	material,
	crack_mesh,
	boundary_mesh,
	dirichlet_boundary_condition,
	neumann_boundary_condition,
	traction_crack_condition,
	SIFs,
)
@profview CrackBEM.solve_static!(p, meshsize; qorder = 4)

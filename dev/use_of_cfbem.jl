using Inti, Gmsh
using Revise
using CrackFastBEM
using StaticArrays
using Elliptic
using Statistics

using Pkg
Pkg.activate(".")

function main()
	meshsize = 0.2
	power = 3

	radius_x = 1.0
	radius_y = 1.0

	gmsh.initialize(String[], true)
	gmsh.option.setNumber("Mesh.MeshSizeMin", meshsize)
	gmsh.option.setNumber("Mesh.MeshSizeMax", meshsize)
	# gmsh.option.setNumber("Mesh.Algorithm", 11)
	gmsh.model.occ.addDisk(0.0, 0.0, 0.0, radius_x, radius_y)
	gmsh.model.occ.synchronize()
	gmsh.model.mesh.generate(2)
	gmsh.model.mesh.recombine()
	gmsh.model.mesh.setOrder(1)
	msh = Inti.import_mesh(; dim = 3)
	# gmsh.fltk.run()
	gmsh.finalize()

	Γ = Inti.Domain(e -> Inti.geometric_dimension(e) == 2, msh)

	crack_mesh = msh[Γ]

	μ = 1
	ν = 0.15

	material = CrackFastBEM.Material(
		"A-dimensional material",
		2 * μ * (1 + ν),
		ν,
		0.0,
		1.0,
	)

	boundary_mesh = CrackFastBEM.empty_mesh()
	dirichlet_boundary_condition = Dict{Int, SVector{3, Float64}}()
	neumann_boundary_condition = Dict{Int, SVector{3, Float64}}()
	σ_inf = 1.0
	traction_crack_condition = Dict(node_id => -SVector(0.0, 0.0, σ_inf) for (node_id, node) in enumerate(crack_mesh.nodes))
	SIFs = Dict{Int, SVector{3, Float64}}()

	p = CrackFastBEM.CrackProblem(
		"Test problem",
		material,
		crack_mesh,
		boundary_mesh,
		dirichlet_boundary_condition,
		neumann_boundary_condition,
		traction_crack_condition,
		SIFs,
	)
	weight_func = CrackFastBEM.create_elliptic_weight_function(radius_x, radius_y; d_inf = 2 * meshsize, d_sup = 4 * meshsize)
	CrackFastBEM.solve_static!(p; qorder = 2, weight_function = weight_func)
	fig = CrackFastBEM.plot_sifs(p)
	# display(fig)
	SIFs = p.SIFs
	front_nodes = p.crack_front_node
	front_nodes_idxs = p.crack_front_node_id
	function mode_I_sif(θ)
		a = radius_x
		b = radius_y
		k = sqrt(1 - b^2 / a^2)
		E_k = Elliptic.E(k^2)
		tempA = σ_inf * sqrt(π) / (E_k) * sqrt(b / a)
		tempB = (a^4 * sin(θ)^2 + b^4 * cos(θ)^2) / (a^2 * sin(θ)^2 + b^2 * cos(θ)^2)
		return tempA * tempB^(1 / 4)
	end
	analytical_SIFs = CrackFastBEM.nodal_scalar_field_in_XY_plane_from_polar_function(front_nodes, front_nodes_idxs, mode_I_sif)
	error_sifs_mode_I = Dict{Int, Float64}()
	for node_id in front_nodes_idxs
		error_sifs_mode_I[node_id] = abs(SIFs[node_id][1] - analytical_SIFs[node_id]) / abs(analytical_SIFs[node_id])
	end
	mean_error_sifs_mode_I = mean(values(error_sifs_mode_I))
	@info "Mean error on the SIF's mode I = $(mean_error_sifs_mode_I)"
	fig = CrackFastBEM.plot_sifs(p)
	display(fig)
end

main()

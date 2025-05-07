@enum ProblemStatus begin
	Undefined
	IllDefined
	WellDefined
	Unsolved
	Solved
end

abstract type AbstractFrame end

struct CrackKernel{T} <: Inti.AbstractKernel{T}
	K::Inti.AbstractKernel{T}
	w::Function
end

function Base.show(io::IO, ck::CrackKernel)
	print(io, "Weighted kernel corresponding to the kernel: \n")
	print(io, ck.K, "\n")
end

function (ck::CrackKernel{T})(source, target)::T where {T}
	return ck.K(source, target) * ck.w(source, target)
end

Inti.singularity_order(ck::CrackKernel) = Inti.singularity_order(ck.K)

struct Frame{T} <: AbstractFrame
	O::T
	n::T
	ν::T
	τ::T
end

function Base.show(io::IO, f::Frame)
	print(io, "Frame: \n")
	print(io, "O: ", f.O, "\n")
	print(io, "n: ", f.n, "\n")
	print(io, "ν: ", f.ν, "\n")
	print(io, "τ: ", f.τ, "\n")
end

abstract type AbstractMaterial end
abstract type AbstractCrackProblem end

struct Material <: AbstractMaterial
	name::String
	young_modulus::Float64
	poisson_ratio::Float64
	thermal_expansion::Float64
	density::Float64
end

function Base.getproperty(mat::Material, s::Symbol)
	if s == :lambda
		return mat.young_modulus * mat.poisson_ratio / ((1 + mat.poisson_ratio) * (1 - 2 * mat.poisson_ratio))
	elseif s == :mu
		return mat.young_modulus / (2 * (1 + mat.poisson_ratio))
	else
		return getfield(mat, s)
	end
end

mutable struct CrackProblem <: AbstractCrackProblem
	name::String
	material::Material
	crack_mesh::Inti.Mesh{3, Float64}
	boundary_mesh::Inti.Mesh{3, Float64}
	dirichlet_boundary_condition::Dict{Int, SVector{3, Float64}} ### Dirichlet boundary condition at each node mesh
	neumann_boundary_condition::Dict{Int, SVector{3, Float64}} ### Neumann boundary condition at each node mesh
	traction_crack_condition::Dict{Int, SVector{3, Float64}} ### Traction at each node mesh of the crack
	SIFs::Dict{Int, SVector{3, Float64}}
	status::Enum

	function CrackProblem(
		name::String,
		material::Material,
		crack_mesh::Inti.Mesh{3, Float64},
		boundary_mesh::Inti.Mesh{3, Float64},
		dirichlet_boundary_condition::Dict{Int, SVector{3, Float64}},
		neumann_boundary_condition::Dict{Int, SVector{3, Float64}},
		traction_crack_condition::Dict{Int, SVector{3, Float64}},
		SIFs::Dict{Int, SVector{3, Float64}},
	)
		p = new(
			name,
			material,
			crack_mesh,
			boundary_mesh,
			dirichlet_boundary_condition,
			neumann_boundary_condition,
			traction_crack_condition,
			SIFs,
			Undefined,
		)
		p.status = check_validity(p)
		return p
	end
end

function check_validity(p::CrackProblem)
	if (length(p.crack_mesh.nodes) != length(p.traction_crack_condition))
		status = IllDefined
		@warn "Crack problem is not well defined."
	else
		status = WellDefined
	end
	return status
end

function Base.getproperty(p::CrackProblem, s::Symbol)
	if s == :crack_front_node_id
		return get_border_nodes_idxs(p.crack_mesh)
	elseif s == :crack_front_node
		return get_border_nodes(p.crack_mesh)
	elseif s == :crack_front_edges
		return ordered_border_edges(p.crack_mesh)
	else
		return getfield(p, s)
	end
end

function Base.show(io::IO, p::CrackProblem)
	print(io, "CrackProblem: ", p.name, "\n")
	print(io, "Material: ", p.material.name, "\n")
	print(io, "Crack mesh: ", p.crack_mesh, "\n")
	print(io, "Boundary mesh: ", p.boundary_mesh.nodes, "\n")
	print(io, "Dirichlet boundary condition: ", length(p.dirichlet_boundary_condition), " nodes\n")
	print(io, "Neumann boundary condition: ", length(p.neumann_boundary_condition), " nodes\n")
	print(io, "Traction crack condition: ", length(p.traction_crack_condition), " nodes\n")
	print(io, "SIFs: ", length(p.SIFs), " nodes\n")
	print(io, "Status: ", p.status, "\n")
end

get_border_nodes_idxs(p::CrackProblem) = p.crack_front_node_id
get_border_nodes(p::CrackProblem) = p.crack_front_node
get_border_edges(p::CrackProblem) = p.crack_front_edges

function find_edges_from_node_id(crack_front_edges::Dict{DataType, Dict{Vector{Int}, Int}}, node_id::Int)
	edges = Dict{DataType, Vector{Vector{Int}}}()
	for (E, d) in crack_front_edges
		edges_list = Vector{Vector{Int}}()
		for (edge, el_id) in d
			if node_id in edge
				edges_list = push!(edges_list, edge)
			end
		end
		if length(edges_list) > 0
			edges[E] = edges_list
		end
	end
	return edges
end

function frenet_frames(p)
	frames = Dict{Int, Frame}()
	crack_front_node_id = p.crack_front_node_id
	crack_front_edges = p.crack_front_edges
	mesh = p.crack_mesh
	for node_id in crack_front_node_id
		# call from the mesh_api file
		frames[node_id] = frenet_frame_from_node_id(mesh, crack_front_edges, node_id)
	end
	return frames
end

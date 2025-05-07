function solve!(u, T₀, δT, t)
	L_ = LinearMap{Float64}(3 * size(T₀, 1)) do y, x
		σ = reinterpret(SVector{3, Float64}, x)
		μ = reinterpret(SVector{3, Float64}, y)
		mul!(μ, T₀, σ)
		mul!(μ, δT, σ, 1, 1)
		return y
	end
	u_ = reinterpret(Float64, u)
	t_ = reinterpret(Float64, t)
	u_, gmres_hist = gmres!(u_, L_, t_, restart = 1000, maxiter = 1000, log = true)
	@info "gmres : $gmres_hist"
	return u
end

solve(T₀, δT, t) = solve!(zero(t), T₀, δT, t)

function mat_frenet2cart(frame::Frame)
	τ = frame.τ
	ν = frame.ν
	n = frame.n
	P = SMatrix{3, 3, Float64}(
		n[1], n[2], n[3],
		τ[1], τ[2], τ[3],
		ν[1], ν[2], ν[3],
	)
	return P
end

function mat_cart2frenet(frame::Frame)
	return mat_frenet2cart(frame)'
end

function cart2frenet(frame::Frame, x::SVector{3, Float64})
	P = mat_cart2frenet(frame)
	return P * x
end

function precompute_sifs(μ, ν)::SVector{3, Float64}
	k1 = μ / (4 * (1 - ν)) * sqrt(2π)
	k2 = μ / 4 * sqrt(2π)
	return SVector(k1, k1, k2)
end

function compute_sifs(p::CrackProblem, solution_nvals::AbstractVector; weighted::Bool = true)
	crack_front_edges = p.crack_front_edges
	μ = p.material.mu
	ν = p.material.poisson_ratio
	pre_sifs = precompute_sifs(μ, ν)
	front_nodes_idxs = p.crack_front_node_id
	frenet_frs = frenet_frames(p)
	crack_front_nvals = Dict(node_id => solution_nvals[node_id] for node_id in front_nodes_idxs)
	SIFs = Dict{Int, SVector{3, Float64}}()
	for (node_id, frenet_frame) in pairs(frenet_frs)
		φ_frenet = cart2frenet(frenet_frame, crack_front_nvals[node_id])
		if weighted
			SIFs[node_id] = pre_sifs .* φ_frenet
		else
			throw("SIFs computation not implemented yet for non-weighted solutions")
		end
	end
	return SIFs
end

"""
function nodal_scalar_field(mesh::Inti.Mesh, node_idxs::Vector{Int}, F::Function)::Vector{Float64}
	From a function F called like F(θ), where θ ∈ [0, 2π], compute the value of F at each node of the mesh.
"""
function nodal_scalar_field_in_XY_plane_from_polar_function(nodes::Vector{SVector{3, Float64}}, nodes_idxs::Vector{Int}, F::Function)::Dict{Int, Float64}
	res = Dict{Int, Float64}()
	for (node, node_id) in zip(nodes, nodes_idxs)
		theta = atan(node[2], node[1])
		if theta < 0
			theta += 2 * π
		end
		res[node_id] = F(theta)
	end
	return res
end

function solve_static!(p::CrackProblem; qorder = 2, weight_function = CrackFastBEM.create_weight_function(p, rtol = 0.01, threshold_inf = 0.1, threshold_sup = 0.2))
	@info "Solving static problem..."
	@info "$(p)"
	# Extract material properties
	λ = p.material.lambda
	μ = p.material.mu
	E = p.material.young_modulus
	ν = p.material.poisson_ratio

	# Extract meshes
	crack_mesh = p.crack_mesh
	boundary_mesh = p.boundary_mesh

	# Create quadratures
	Q_crack = Inti.Quadrature(crack_mesh; qorder)
	Q_boundary = Inti.Quadrature(boundary_mesh; qorder)

	# Create elastic operator
	op = Inti.Elastostatic(; λ, μ, dim = 3)
	elastic_kernel = Inti.HyperSingularKernel(op)

	# Create the kernel function for the crack problem, witch is a weight function × the elastic kernel

	front_nodes = p.crack_front_node
	front_nodes_idxs = p.crack_front_node_id
	crack_kernel = CrackKernel(elastic_kernel, weight_function)

	# crack_kernel = (source, target) -> begin
	# x, y = Inti.coords(source), Inti.coords(target)
	# return w(source, target) * elastic_kernel(target, source)::SMatrix{3, 3, Float64, 9}
	# end

	Hop = Inti.IntegralOperator(crack_kernel, Q_crack)

	# Assemble the matrix
	H₀ = Inti.assemble_hmatrix(Hop)
	@info "𝓗-Matrix compression ratio = $(HMatrices.compression_ratio(H₀))"
	δH = Inti.adaptive_correction(Hop)

	t = node_vals_to_quadrature(Q_crack, collect(values(p.traction_crack_condition)))
	φ = solve(H₀, δH, t)

	@info "Solving static problem... done"

	@info "Mean COD = $(mean(φ))"

	nvals = Inti.quadrature_to_node_vals(Q_crack, φ)

	@info "Computing SIFs..."

	SIFs = compute_sifs(p, nvals; weighted = true)

	@info "Computing SIFs... done"

	p.SIFs = SIFs
	p.status = Solved
end

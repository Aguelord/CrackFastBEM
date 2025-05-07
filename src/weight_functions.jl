"""
function compute_distance_to_nodes(y::SVector{3, Float64}, nodes_list::AbstractVector)

	Given a point `y` and a list of nodes, compute the distance from `y` to the nearest node in the list.
"""
function compute_distance_to_nodes(y::SVector{3, Float64}, nodes_list::AbstractVector)
	d_list = [norm(y - node) for node in nodes_list]
	return minimum(d_list)
end

"""
function distance_and_nearest_node_and_node_id(y::SVector{3, Float64}, nodes_list::AbstractVector, nodes_idxs_list::AbstractVector)

	Given a point `y`, a list of nodes, and their corresponding indices, compute the distance from `y` to the nearest node in the list.
	Returns the distance, the nearest node, and its index.
"""
function distance_and_nearest_node_and_node_id(y::SVector{3, Float64}, nodes_list::AbstractVector, nodes_idxs_list::AbstractVector)
	dists = map(p -> norm(p .- y), nodes_list)
	min_dist, local_idx = findmin(dists)
	return min_dist, nodes_list[local_idx], nodes_idxs_list[local_idx]
end

function smoothing_function(x, x_inf, x_sup)
	if x < x_inf
		return 1.0
	elseif x > x_sup
		return 0.0
	else
		u = (x - x_inf) / (x_sup - x_inf)
		return exp(2 * exp(-1 / u) / (u - 1))
	end
end

function old_create_weight_function(p::CrackProblem, front_nodes, front_nodes_idxs; treshold_inf, treshold_sup)
	@info "Create (old) weight function..."
	frames = frenet_frames(p)
	return (source, target) -> begin
		y = target.coords
		D, nearest_node, nearest_node_id = distance_and_nearest_node_and_node_id(y, front_nodes, front_nodes_idxs)
		τ = frames[nearest_node_id].τ
		d = sqrt(D^2 - dot(y - nearest_node, τ)^2)
		f = smoothing_function(d, treshold_inf, treshold_sup)
		return (sqrt(d) * f + (1 - f))::Float64
	end
end

function create_weight_function(p::CrackProblem; rtol, threshold_inf = 0.1, threshold_sup = 0.2)
	front_elements_dict = build_front_mesh_from_edges(p.crack_mesh)
	distance_func = dist_to_msh_func(front_elements_dict; kneighbors = 1, rel_tol = rtol)
	return (source, target) -> begin
		d = distance_func(target.coords)
		f = smoothing_function(d, threshold_inf, threshold_sup)
		return (sqrt(d) * f + (1 - f))::Float64
	end
end

function create_custom_weight_function(; d_inf, d_sup)
	return (source, target) -> begin
		a = 2.0
		b = 1.0
		x, y = Inti.coords(source), Inti.coords(target)
		d = ((a * b) / (sqrt(b^2 * y[1]^2 + a^2 * y[2]^2)) - 1) * (b^2 * y[1]^2 + a^2 * y[2]^2) / (sqrt(b^4 * y[1]^2 + a^4 * y[2]^2))
		if d < 0
			d = 0.0
		end
		f = smoothing_function(d, d_inf, d_sup)
		return (sqrt(d) * f + (1 - f))::Float64
	end
end

function create_elliptic_weight_function(a, b; d_inf, d_sup)
	return (source, target) -> begin
		x, y = Inti.coords(source), Inti.coords(target)
		d = ((a * b) / (sqrt(b^2 * y[1]^2 + a^2 * y[2]^2)) - 1) * (b^2 * y[1]^2 + a^2 * y[2]^2) / (sqrt(b^4 * y[1]^2 + a^4 * y[2]^2))
		if d < 0
			d = 0.0
		end
		f = smoothing_function(d, d_inf, d_sup)
		return (sqrt(d) * f + (1 - f))::Float64
	end
end

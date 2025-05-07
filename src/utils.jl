#=

Utility functions that have nowhere obvious to go.

=#

"""
	hermitian_distance(y::SVector{3, Float64}, first_node::SVector{3, Float64},
					   second_node::SVector{3, Float64},
					   first_tangent::SVector{3, Float64},
					   second_tangent::SVector{3, Float64})

Given two nodes and their tangents, compute the distance between any point `y` and the Hermitian curve 
defined between the two nodes. This uses the Newton algorithm to find the parameter `u ∈ [0,1]` that 
minimizes the squared distance `‖r(u) - y‖²`, where `r(u)` is the Hermite cubic interpolant.

Returns a tuple `(distance::Float64, u::Float64)`.
"""
function hermitian_distance(y::SVector{3, Float64},
	first_node::SVector{3, Float64},
	second_node::SVector{3, Float64},
	first_tangent::SVector{3, Float64},
	second_tangent::SVector{3, Float64},
)::Tuple{Float64, Float64}

	# Hermite basis functions and their derivatives
	h00(u) = 2u^3 - 3u^2 + 1
	h10(u) = u^3 - 2u^2 + u
	h01(u) = -2u^3 + 3u^2
	h11(u) = u^3 - u^2

	dh00(u) = 6u^2 - 6u
	dh10(u) = 3u^2 - 4u + 1
	dh01(u) = -6u^2 + 6u
	dh11(u) = 3u^2 - 2u

	ddh00(u) = 12u - 6
	ddh10(u) = 6u - 4
	ddh01(u) = -12u + 6
	ddh11(u) = 6u - 2

	# Newton initialization
	u = 0.5
	max_iter = 10
	tol = 1e-10

	for _ in 1:max_iter
		r = h00(u) * first_node + h10(u) * first_tangent +
			h01(u) * second_node + h11(u) * second_tangent

		r1 = dh00(u) * first_node + dh10(u) * first_tangent +
			 dh01(u) * second_node + dh11(u) * second_tangent

		r2 = ddh00(u) * first_node + ddh10(u) * first_tangent +
			 ddh01(u) * second_node + ddh11(u) * second_tangent

		diff = r - y
		f1 = 2 * dot(diff, r1)
		f2 = 2 * (dot(r1, r1) + dot(diff, r2))

		if abs(f2) < 1e-12
			break  # avoid division by near-zero
		end

		delta = f1 / f2
		u = clamp(u - delta, 0.0, 1.0)

		if abs(delta) < tol
			break
		end
	end

	# Final projection point
	r_final = h00(u) * first_node + h10(u) * first_tangent +
			  h01(u) * second_node + h11(u) * second_tangent

	return r_final, norm(r_final - y), u
end

"""
function dist_to_msh_func(msh::Inti.AbstractMesh; kneighbors = 1)

	Build a function that computes the distance from a point to the mesh `msh`.
	The function uses a KDTree to find the nearest elements and then computes the 
	distance to the closest point on each element.
	
	Arguments:
	- `msh`: The mesh to compute distances to.
	- `kneighbors`: The number of nearest neighbors to consider (default is 1).
	
	Returns:
	A function that takes a point `x` and returns the distance from `x` to the mesh.
"""
function dist_to_msh_func(els_dict::Dict; kneighbors = 1, rel_tol = sqrt(eps()))
	els = els_dict |> values |> collect
	centers = map(x -> Inti.center(x), els)
	kdtree = KDTree(centers)
	idxs = Vector{Int}(undef, kneighbors)
	dists = Vector{Float64}(undef, kneighbors)
	method = Fminbox()
	options = Optim.Options(iterations = 1)
	lb = [0.0]
	ub = [1.0]
	initial_y = [0.5]
	return function dist_to_msh(x)
		knn!(idxs, dists, kdtree, x, kneighbors)
		d2min = Inf
		for i in idxs
			el = els[i]
			f(y) = transpose(el(y) - x) * (el(y) - x)
			function g!(storage, y)
				storage = 2 * transpose((el(y) - x)) * Inti.jacobian(el, y)
			end
			od = OnceDifferentiable(f, g!, initial_y)
			result = optimize(od, lb, ub, initial_y, method, options)
			d2min = min(d2min, result.minimum)
		end
		return sqrt(d2min)
	end
end

"""
function suppress_output(f::Function)

	Redirects the standard output and error streams to `/dev/null` while executing the function `f`.
	This is useful for suppressing output during function execution.
	
	Arguments:
	- `f`: The function to execute with suppressed output.
	
	Returns:
	The result of the function `f`.
"""
function suppress_output(f::Function)
	redirect_stdout(devnull) do
		redirect_stderr(devnull) do
			return f()
		end
	end
end

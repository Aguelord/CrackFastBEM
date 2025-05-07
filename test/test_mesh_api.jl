using Test
using Inti, Gmsh
using CrackFastBEM

@testset "Mesh nodes to quadrature nodes extrapolation" begin
	meshsize = 0.5
	radius_x = 2.0
	radius_y = 1.0

	qorder = 2

	p = CrackFastBEM.suppress_output() do
		CrackFastBEM.elliptic_crack_in_infinite_domain_under_uniform_normal_loading_p1_square(radius_x, radius_y, meshsize)
	end
	Q = Inti.Quadrature(p.crack_mesh; qorder = 2)
	n_tests = 100
	for itest in 1:n_tests
		coeffs = randn(4)
		function random_polynomial_function(x::SVector{3, Float64})
			return coeffs[1] + coeffs[2] * x[1] + coeffs[3] * x[2] + coeffs[4] * x[3]
		end
		nvals = [random_polynomial_function(node) for node in p.crack_mesh.nodes]
		qvals = CrackFastBEM.node_vals_to_quadrature(Q, nvals)
		nvals_bis = Inti.quadrature_to_node_vals(Q, qvals)
		error = norm((nvals_bis .- nvals) ./ nvals)
		@test error < 1e-10
	end
end

@testset "Border mesh nodes" begin
	@testset "Ellipse" begin
		meshsize = 0.2
		radius_x = 1.0
		radius_y = 1.0
		p = CrackFastBEM.suppress_output() do
			CrackFastBEM.elliptic_crack_in_infinite_domain_under_uniform_normal_loading_p1_square(radius_x, radius_y, meshsize)
		end
		border_nodes = CrackFastBEM.get_border_nodes(p.crack_mesh)
		for node in border_nodes
			x, y = node[1], node[2]
			@test (x^2 / radius_x^2 + y^2 / radius_y^2) ≈ 1.0
		end
	end
end

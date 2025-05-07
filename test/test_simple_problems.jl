using Test
using CrackFastBEM
using Inti
using Statistics

@testset "Unit disk" begin
	radius = 1.0
	meshsize = 0.2
	SIFs = CrackFastBEM.suppress_output() do
		p = CrackFastBEM.circular_crack_in_infinite_domain_under_uniform_normal_loading_p1_square(radius, meshsize)
		Q = Inti.Quadrature(p.crack_mesh; qorder = 2)
		weight_func = CrackFastBEM.create_elliptic_weight_function(radius, radius; d_inf = 2 * meshsize, d_sup = 4 * meshsize)

		CrackFastBEM.solve_static!(p; qorder = 2, weight_function = weight_func)
		SIFs = p.SIFs
	end
	mean_sif = mean([SIFs[i][1] for i in keys(SIFs)])
	theorical_mean_sif = 2 / sqrt(π)
	@info "theorical SIF = $(theorical_mean_sif)"
	@info "numerical SIF = $(mean_sif)"
	@info "Mean error on the SIF's mode I = $(abs((mean_sif - theorical_mean_sif) / theorical_mean_sif))"
	@test abs((mean_sif - theorical_mean_sif) / theorical_mean_sif) < 1e-2
end

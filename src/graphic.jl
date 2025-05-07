function plot_sifs(p::CrackProblem)
	SIFs = p.SIFs
	ordered_nodes = order_nodes_idxs(p.crack_mesh)
	SIFs_vals = [SIFs[node_id] for node_id in ordered_nodes]
	if isempty(SIFs_vals)
		error("SIFs are not computed yet. Please compute them first before plotting them.")
	end
	Y = collect(values(SIFs_vals))
	fig = Figure()
	ax = Axis(fig[1, 1], title = "SIFs")
	x = LinRange(0, length(SIFs_vals) - 1, length(SIFs_vals))
	y1 = [Y[i][1] for i in 1:length(SIFs_vals)]
	y2 = [Y[i][2] for i in 1:length(SIFs_vals)]
	y3 = [Y[i][3] for i in 1:length(SIFs_vals)]
	lines!(ax, x, y1, label = "K_I", color = :red)
	lines!(ax, x, y2, label = "K_II", color = :green)
	lines!(ax, x, y3, label = "K_III", color = :blue)
	axislegend(ax, position = (:right, :top))
	return fig
end

function viz_on_mesh(p::CrackProblem, data::Vector{Float64}, title::String)
	msh = p.crack_mesh
	fig = Figure()
	ax = Axis3(fig[1, 1])
	viz!(msh; color = data, interpolate = false, showsegments = true)
	colorrange = extrema(data)
	cb = Colorbar(fig[1, 2]; label = title, colorrange)
	return fig
end

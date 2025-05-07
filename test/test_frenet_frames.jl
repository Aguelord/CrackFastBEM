using Test
using Revise
using CrackFastBEM
using Inti, Gmsh
using StaticArrays
using LinearAlgebra

meshsize = 0.1

@testset "Frenet frames and border nodes" begin
	for elementOrder in 1:2
		@testset "On a disk, element order = $(elementOrder)" begin
			r = 1.0
			p = CrackFastBEM.suppress_output() do
				CrackFastBEM.circular_crack_in_infinite_domain_under_uniform_normal_loading_p1_square(r, meshsize)
			end
			frames = CrackFastBEM.frenet_frames(p)
			for (node_id, frame) in frames
				node = p.crack_mesh.nodes[node_id]
				@test frame.O ≈ node
				@test frame.n ≈ SVector(0.0, 0.0, 1.0)
				@test frame.τ ≈ frame.n × frame.ν
				@test abs.(r .* frame.ν - node) < meshsize * SVector(1.0, 1.0, 1.0)
			end
		end
		@testset "On an ellipse, element order = $(elementOrder)" begin
			radius_x = 1.5
			radius_y = 1.0
			p = CrackFastBEM.suppress_output() do
				CrackFastBEM.elliptic_crack_in_infinite_domain_under_uniform_normal_loading_p1_square(radius_x, radius_y, meshsize)
			end
			frames = CrackFastBEM.frenet_frames(p)
			@testset "Border nodes are correct ?" begin
				for (node_id, frame) in frames
					node = p.crack_mesh.nodes[node_id]
					@test radius_x * radius_y / sqrt(radius_x^2 * node[2]^2 + radius_y^2 * node[1]^2) ≈ 1
				end
			end
			@testset "Frames are correct ?" begin
				for (node_id, frame) in frames
					node = p.crack_mesh.nodes[node_id]
					@test frame.O ≈ node
					@test frame.n ≈ SVector(0.0, 0.0, 1.0)
					N = sqrt(radius_x^4 * node[2]^2 + radius_y^4 * node[1]^2)
					T = SVector(
						-radius_x^2 * node[2],
						radius_y^2 * node[1],
						0.0) / N
					@test ((frame.τ - T) < meshsize * SVector(1.0, 1.0, 1.0) || (frame.τ + T) < meshsize * SVector(1.0, 1.0, 1.0))
					ν = SVector(
						-T[2],
						T[1],
						0.0)
					@test ((frame.ν - ν) < meshsize * SVector(1.0, 1.0, 1.0) || (frame.ν + ν) < meshsize * SVector(1.0, 1.0, 1.0))
				end
			end
		end
	end
end

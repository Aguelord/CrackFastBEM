using CrackFastBEM
using SafeTestsets
using Test

@testset "Intitialization tests" begin
	include("test_basic.jl")
end

@testset "Mesh tests" begin
	include("test_mesh_api.jl")
end

@testset "Pre-processor tests" begin
	include("test_frenet_frames.jl")
end

@testset "General tests" begin
	include("test_inti_for_cracks.jl")
end

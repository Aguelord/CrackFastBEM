"""
Extension of the `reference_interpolation.jl` file in the `Inti` package.
"""
# P2 for ReferenceSquare
function Inti.reference_nodes(::Type{<:Inti.LagrangeSquare{9}})
	return SVector(
		SVector(0, 0),
		SVector(1, 0),
		SVector(1, 1),
		SVector(0, 1),
		SVector(0.5, 0),
		SVector(1, 0.5),
		SVector(0.5, 1),
		SVector(0, 0.5),
		SVector(0.5, 0.5),
	)
end

function (el::Inti.LagrangeElement{Inti.ReferenceSquare, 9})(u)
	v = Inti.vals(el)
	L1 = (1 - u[1]) * (1 - u[2]) * (1 - 2 * u[1]) * (1 - 2 * u[2])
	L2 = u[1] * (1 - u[2]) * (2 * u[1] - 1) * (1 - 2 * u[2])
	L3 = u[1] * u[2] * (2 * u[1] - 1) * (2 * u[2] - 1)
	L4 = (1 - u[1]) * u[2] * (1 - 2 * u[1]) * (2 * u[2] - 1)
	L5 = 4 * u[1] * (1 - u[2]) * (1 - u[1]) * (1 - 2 * u[2])
	L6 = 4 * u[1] * u[2] * (2 * u[1] - 1) * (1 - u[2])
	L7 = 4 * u[1] * u[2] * (1 - u[1]) * (2 * u[2] - 1)
	L8 = 4 * (1 - u[1]) * u[2] * (1 - 2 * u[1]) * (1 - u[2])
	L9 = 16 * u[1] * u[2] * (1 - u[1]) * (1 - u[2])
	return v[1] * L1 +
		   v[2] * L2 +
		   v[3] * L3 +
		   v[4] * L4 +
		   v[5] * L5 +
		   v[6] * L6 +
		   v[7] * L7 +
		   v[8] * L8 +
		   v[9] * L9
end

abstract type AbstractEdge end

struct Edge{T} <: AbstractEdge
	nodes::Vector{T}
end

function Base.:(==)(e1::Edge, e2::Edge)
	return e1.nodes == e2.nodes || e1.nodes == reverse(e2.nodes)
end

function Base.hash(e::Edge, h::UInt)
	n = e.nodes
	return min(hash(n, h), hash(reverse(n), h))
end

"""
function nods_idx(connectivity_matrix)

	Returns the indices of the nodes from a connectivity matrix.
"""
function nodes_idx(connectivity_matrix::Matrix{Int})
	_nodes_idx = Int[]
	for col in 1:size(connectivity_matrix, 2)
		element_node_indices = connectivity_matrix[:, col]
		_nodes_idx = _nodes_idx ∪ element_node_indices
	end
	return unique(_nodes_idx)
end

"""
function Inti.Mesh(etype::DataType, connectivity::Matrix{Int}, nodes::Vector{SVector{3, Float64}})

	Create a Mesh object in a very specific situation where there is only one element type `etype` throughout the whole mesh.
	WARNING : When the connectivity matrix contains nodes idx that are not numbering as 1, 2, 3 etc... (for example 1, 2, 4, 6, ...) DON'T USE IT, it causes a problem when creating the element iterator.
"""
function Inti.Mesh(;
	etype::DataType,
	nodes,
	connectivity::Matrix{Int},
)
	n = length(Inti.reference_nodes(etype))
	@assert size(connectivity, 1) == n "connectivity matrix must have $n rows"
	@assert length(nodes_idx(connectivity)) == length(nodes) "connectivity matrix must be consistent with the nodes"
	etype2mat = Dict(etype => connectivity)
	etype2els = Dict{DataType, AbstractVector}()
	tag = Inti.new_tag(2)
	ent = Inti.GeometricEntity(2, tag, Inti.EntityKey[], String[], nothing)
	etags = collect(1:size(connectivity, 2))
	ent2etag = Dict(ent => Dict(etype => etags))
	msh = Inti.Mesh{3, Float64}(nodes, etype2mat, etype2els, ent2etag)
	msh.etype2els[etype] = Inti.ElementIterator(msh, etype)
	return msh
end

empty_mesh() = Inti.Mesh{3, Float64}()
isempty_msh(msh::Inti.Mesh{3, Float64}) = isempty(msh.nodes)

"""
function invert_connectivity(connectivity::Matrix{Int})::Dict{Int, Vector{Int}}

	Inverts the connectivity matrix to create a dictionary that maps each node to the elements it belongs to.
	Returns a dictionary where the keys are node indices and the values are vectors of element indices.
"""
function invert_connectivity(connectivity::Matrix{Int})::Dict{Int, Vector{Int}}
	inverted_connectivity = Dict{Int, Vector{Int}}()
	for i in 1:size(connectivity, 2)
		nodes = connectivity[:, i]
		for node in nodes
			if !haskey(inverted_connectivity, node)
				inverted_connectivity[node] = Int[]
			end
			push!(inverted_connectivity[node], i)
		end
	end
	return inverted_connectivity
end

"""
function get_edge_nodes_idxs(E::DataType)

	Returns the local indices of the nodes that are on an edge of an element for a given element type `E`.
"""
get_edge_nodes_idxs(::Type{<:Inti.LagrangeSquare{4}}) = 1:4
get_edge_nodes_idxs(::Type{<:Inti.LagrangeSquare{8}}) = 1:8
get_edge_nodes_idxs(::Type{<:Inti.LagrangeSquare{9}}) = 1:8
get_edge_nodes_idxs(::Type{<:Inti.LagrangeTriangle{3}}) = 1:3
get_edge_nodes_idxs(::Type{<:Inti.LagrangeTriangle{4}}) = 1:3
get_edge_nodes_idxs(::Type{<:Inti.LagrangeTriangle{6}}) = 1:6
get_edge_nodes_idxs(::Type{<:Inti.LagrangeTriangle{7}}) = 1:6

"""
function get_edges_local_node_indices(E::DataType)
	Renvoie la liste des arêtes locales sous forme de tuples d'indices de nœuds locaux.
	Les indices respectent la numérotation classique éléments finis comme définis dans gmsh.
"""
get_edges_local_node_indices(::Type{<:Inti.LagrangeSquare{4}}) = [[1, 2], [2, 3], [3, 4], [4, 1]]
get_edges_local_node_indices(::Type{<:Inti.LagrangeSquare{8}}) = [[1, 5, 2], [2, 6, 3], [3, 7, 4], [4, 8, 1]]
get_edges_local_node_indices(::Type{<:Inti.LagrangeSquare{9}}) = get_edges_local_node_indices(Inti.LagrangeSquare{8})
get_edges_local_node_indices(::Type{<:Inti.LagrangeTriangle{3}}) = [[1, 2], [2, 3], [3, 1]]
get_edges_local_node_indices(::Type{<:Inti.LagrangeTriangle{6}}) = [[1, 4, 2], [2, 5, 3], [3, 6, 1]]
get_edges_local_node_indices(::Type{<:Inti.LagrangeTriangle{7}}) = get_edges_local_node_indices(Inti.LagrangeTriangle{6})

"""
function border_connectivity_matrix(E::DataType, connectivity_matrix::Matrix{Int})

	Return the submatrix of the connectivity matrix that contains only the nodes that are on the edges of the elements of type `E`.
"""

function border_connectivity_matrix(E::DataType, connectivity_matrix::Matrix{Int})
	I = get_edge_nodes_idxs(E)
	return connectivity_matrix[I, :]
end

"""
function get_edges(connectivity_matrix::Matrix{Int}, E, el_id::Int)::Vector{Edge{Int}}

	For a given element id `el_id`, return the edges of the element of type `E`.
"""
function get_edges(connectivity_matrix::Matrix{Int}, E, el_id::Int)::Vector{Edge{Int}}
	edges_local_node_indices = get_edges_local_node_indices(E)
	edges = Vector{Edge{Int}}()
	for edge in edges_local_node_indices
		edge = Edge(connectivity_matrix[edge, el_id])
		push!(edges, edge)
	end
	return edges
end

"""
function edged_connectivity(connectivity_matrix::Matrix{Int}, E::DataType)

	For each element of type E, return all the edges of the elements, mapping the element id to the edges.
"""
function edged_connectivity(connectivity_matrix::Matrix{Int}, E)
	# use the previous get_edges function
	conn = Dict{Int, Vector{Edge{Int}}}()
	for el_id in 1:size(connectivity_matrix, 2)
		edges = get_edges(connectivity_matrix, E, el_id)
		conn[el_id] = edges
	end
	return conn
end

"""
function edges_connectivity(connectivity_matrix::Matrix{Int}, E::DataType)

	For each edge of an element of type E, return the element (of type E) ids that contains the edge.
"""
function edges_connectivity(connectivity_matrix::Matrix{Int}, E)::Dict{Edge{Int}, Vector{Int}}
	# use the previous edged_connectivity function
	edges_conn = Dict{Edge{Int}, Vector{Int}}()
	for el_id in 1:size(connectivity_matrix, 2)
		edges = get_edges(connectivity_matrix, E, el_id)
		for edge in edges
			if !haskey(edges_conn, edge)
				edges_conn[edge] = Int[]
			end
			push!(edges_conn[edge], el_id)
		end
	end
	for (key, value) in edges_conn
		edges_conn[key] = unique(value)  # Supprime les doublons de chaque liste d'arêtes
	end
	return edges_conn
end

function count_adjacent_elements(msh::Inti.Mesh{3, Float64})::Dict{Edge{Int}, Int}
	# Count the number of elements adjacent to the edge
	count = Dict{Edge{Int}, Int}()
	for (E, connectivity_matrix) in msh.etype2mat
		edges_conn = edges_connectivity(connectivity_matrix, E)
		for (edge, el_ids) in edges_conn
			if !haskey(count, edge)
				count[edge] = 0
			end
			count[edge] += length(el_ids)
		end
	end
	return count
end

function edges_on_the_border(msh::Inti.Mesh{3, Float64})::Dict{Edge{Int}, Int}
	count = count_adjacent_elements(msh)
	filtered_dict = Dict(k => v for (k, v) in count if v == 1)
	return filtered_dict
end

"""
function get_border_edges(msh::Inti.Mesh{3, Float64})

	Return a dictionnary that maps every element type of the mesh to the edges connected to only one element of that type.
	For each element type, the edges are represented as a dictionnary where the keys are the edges and the values are the element ids.
"""
function get_border_edges(msh::Inti.Mesh{3, Float64})
	connectivity_matrices = msh.etype2mat
	border_edges = Dict{DataType, Dict{Edge{Int}, Int}}()
	count = count_adjacent_elements(msh)
	for (E, connectivity_matrix) in connectivity_matrices
		conn_edges = edges_connectivity(connectivity_matrix, E)
		border_edges[E] = Dict{Edge{Int}, Int}()
		for (edge, el_ids) in conn_edges
			if count[edge] == 1
				border_edges[E][edge] = el_ids[1]
			end
		end
	end
	return border_edges
end

"""
function ordered_border_edges(raw_edges::Dict{DataType, Dict{Edge{Int}, Int}},)::Dict{DataType, Dict{Vector{Int}, Int}}

	Return the ordered edges dictionnay, meaning that every edge is well oriented with all the others.
"""
function ordered_border_edges(
	raw_edges::Dict{DataType, Dict{Edge{Int}, Int}},
)::Dict{DataType, Dict{Vector{Int}, Int}}

	# Étape 1 — Fusionner les arêtes non orientées avec leur type et id
	merged_edges = Dict{CrackBEM.Edge{Int}, Tuple{DataType, Int}}()
	for (etype, edge_dict) in raw_edges
		for (edge, el_id) in edge_dict
			merged_edges[edge] = (etype, el_id)
		end
	end

	# Étape 2 — Construire la connectivité noeud → arêtes
	conn = Dict{Int, Vector{CrackBEM.Edge{Int}}}()
	for edge in keys(merged_edges)
		for node in edge.nodes
			push!(get!(conn, node, Vector{CrackBEM.Edge{Int}}()), edge)
		end
	end

	function is_true_endpoint(node::Int, edges::Vector{CrackBEM.Edge{Int}})
		length(edges) == 1 || return false  # relié à une seule arête
		edge = edges[1]
		nodes = edge.nodes
		return node == nodes[1] || node == nodes[end]
	end

	# Calcul des "vrais" points d'extrémité
	endpoints = Dict{Int, Vector{CrackBEM.Edge{Int}}}()
	for (node, edges) in conn
		if is_true_endpoint(node, edges)
			endpoints[node] = edges
		end
	end

	# Étape 3 — Trouver un point de départ
	start_edge = begin
		if !isempty(endpoints)
			node = first(keys(endpoints))      # un des noeuds d’extrémité
			only(endpoints[node])              # l’unique arête connectée à ce noeud
		else
			first(keys(merged_edges))          # chemin fermé → arête arbitraire
		end
	end

	# On choisit arbitrairement l’orientation, mais on garde une continuité
	current_nodes = start_edge.nodes
	start_node = current_nodes[1]
	path = [(start_node, current_nodes[2])]
	visited = Set{CrackBEM.Edge{Int}}([start_edge])
	current_node = current_nodes[2]

	if !isempty(endpoints)
		start_node = first(keys(endpoints))
		nodes = start_edge.nodes
		if nodes[1] == start_node
			current_node = nodes[end]
			oriented_nodes = nodes
		elseif nodes[end] == start_node
			current_node = nodes[1]
			oriented_nodes = reverse(nodes)
		else
			error("start_node $start_node n'est pas une extrémité de start_edge")
		end
	else
		nodes = start_edge.nodes
		start_node = nodes[1]
		current_node = nodes[end]
		oriented_nodes = nodes
	end

	visited = Set{CrackBEM.Edge{Int}}([start_edge])
	path = [Tuple(oriented_nodes)]  # stocké sous forme de tuple orienté

	# Étape 4 — Suivre la chaîne
	while true
		candidates = filter(e -> e ∉ visited, conn[current_node])
		isempty(candidates) && break
		next_edge = only(candidates)  # au plus une arête non visitée
		nodes = next_edge.nodes
		if nodes[1] == current_node
			oriented_nodes = nodes
			current_node = nodes[end]
		elseif nodes[end] == current_node
			oriented_nodes = reverse(nodes)
			current_node = nodes[1]
		else
			error("Incohérence : current_node $current_node n’est pas une extrémité de l’arête $(nodes)")
		end
		push!(path, Tuple(oriented_nodes))
		push!(visited, next_edge)
	end

	# Étape 5 — Recréer le dictionnaire trié
	sorted_edges = Dict{DataType, Dict{Vector{Int}, Int}}()
	for oriented_nodes in path
		# Convertir le tuple en vecteur
		v = collect(oriented_nodes)
		# Trouver l’arête correspondante dans merged_edges, quelle que soit l’orientation
		all_edges = collect(keys(merged_edges))
		i = findfirst(e -> e.nodes == v || e.nodes == reverse(v), all_edges)
		i === nothing && error("Impossible de retrouver l’arête pour $v")
		raw_edge = all_edges[i]
		(etype, el_id) = merged_edges[raw_edge]
		if !haskey(sorted_edges, etype)
			sorted_edges[etype] = Dict{Vector{Int}, Int}()
		end
		sorted_edges[etype][v] = el_id  # v est déjà bien orienté
	end
	return sorted_edges
end

"""

function ordered_border_edges(msh::Inti.Mesh{3, Float64})::Dict{DataType, Dict{Vector{Int}, Int}}

	Return the ordered edges dictionnay, meaning that every edge is well oriented with all the others.
"""
function ordered_border_edges(msh::Inti.Mesh{3, Float64})::Dict{DataType, Dict{Vector{Int}, Int}}
	border_edges = get_border_edges(msh)
	return ordered_border_edges(border_edges)
end

"""
function edge2el(etype2edge_dict::Dict{DataType, Dict{Vector{Int}, Int}})::Dict{Vector{Int}, Inti.LagrangeElement}

	Return a dictionnary that maps every edge (in the form of a vector of node ids) to the actual Lagrange element that contains the edge. It is useful because it frees us to think about the element type.
"""
function edge2el(msh::Inti.Mesh, etype2edge_dict::Dict{DataType, Dict{Vector{Int}, Int}})::Dict{Vector{Int}, Inti.LagrangeElement}
	edge2el_dict = Dict{Vector{Int}, Inti.LagrangeElement}()
	for (etype, edge_dict) in etype2edge_dict
		els = Inti.elements(msh, etype)
		for (edge, el_id) in edge_dict
			edge2el_dict[edge] = els[el_id]
		end
	end
	return edge2el_dict
end

"""
function ordered_edges_path(msh::Inti.Mesh{3, Float64})::Vector{Vector{Int}}
	
	Return the ordered edges path, meaning that every edge is well oriented with all the others, and the vector of edges is ordered. You can travel through the edges in the order of the vector, and each edge is well oriented with the next one, in the form : 

	N-element Vector{Vector{Int64}}:
	[4, 5]
	[5, 8]
	[8, 7]
	[7, 6]
	[6, 1]
	[1, 2]
	[2, 3]
	[3, 4]
	if it is a closed path, and :

	[4, 5]
	[5, 8]
	[8, 7]
	[7, 6]
	[6, 1]
	
	for an open path (for example).
"""
function ordered_edges_path(msh::Inti.Mesh{3, Float64})::Vector{Vector{Int}}
	ordered_edges_dict = ordered_border_edges(msh)
	edge_2_el_dict = edge2el(msh, ordered_edges_dict)
	ordered_edges_path = Vector{Vector{Int}}()
	edges = collect(keys(edge_2_el_dict))
	remaining = deepcopy(edges)
	ordered = [popfirst!(remaining)]
	while !isempty(remaining)
		last_node = last(ordered[end])
		found = false
		for (i, e) in pairs(remaining)
			if first(e) == last_node
				push!(ordered, e)
				deleteat!(remaining, i)
				found = true
				break
			elseif last(e) == last_node
				push!(ordered, reverse(e))
				deleteat!(remaining, i)
				found = true
				break
			end
		end
		if !found
			error("Impossible de poursuivre le chemin : aucune correspondance trouvée pour le nœud $last_node")
		end
	end
	return ordered
end

"""
function ordered_nodes_idxs(msh::Inti.Mesh{3, Float64})::Vector{Int}

	Return the global indices of the nodes that are ordered in the mesh.
"""
function order_nodes_idxs(msh::Inti.Mesh{3, Float64})::Vector{Int}
	ordered_edges = ordered_edges_path(msh)
	ordered_nodes = Vector{Int}()
	for edge in ordered_edges
		ordered_nodes = ordered_nodes ∪ edge
	end
	return unique(ordered_nodes)
end

"""
function get_border_nodes_idxs(msh::Inti.Mesh{3, Float64})

	Return the global indices of the nodes that are on the border of the mesh.
"""
function get_border_nodes_idxs(msh::Inti.Mesh{3, Float64})
	ordered_edges_dict = ordered_border_edges(msh)
	border_nodes_idxs = Int[]
	for (E, edges) in ordered_edges_dict
		tuple_node_ids = collect(keys(edges))
		node_ids = collect(Iterators.flatten(tuple_node_ids))
		border_nodes_idxs = border_nodes_idxs ∪ node_ids
	end
	return sort(unique(border_nodes_idxs))
end

"""
function get_border_nodes(msh::Inti.Mesh{3, Float64})

	Return the coordinates of the nodes that are on the border of the mesh.
"""
function get_border_nodes(msh::Inti.Mesh{3, Float64})
	border_nodes_idxs = get_border_nodes_idxs(msh)
	return Inti.nodes(msh)[border_nodes_idxs]
end

function find_elements_from_node_id(crack_front_edges::Dict{DataType, Dict{Vector{Int}, Int}}, node_id::Int)
	elements_id = Dict{DataType, Vector{Int}}()
	for (E, d) in crack_front_edges
		elements_id[E] = Vector{Int}()
		for (edge, el_id) in d
			if node_id in edge
				elements_id[E] = push!(elements_id[E], el_id)
			end
		end
	end
	return elements_id
end

function find_elements_from_node_id_all(crack_front_edges::Dict{DataType, Dict{Vector{Int}, Int}}, node_ids::Vector{Int})
	elements_id = Dict{Int, Dict{DataType, Vector{Int}}}()
	for node_id in node_ids
		elements_id[node_id] = find_elements_from_node_id(crack_front_edges, node_id)
	end
	return elements_id
end

function adjacent_elements_from_node_id(crack_front_edges::Dict{DataType, Dict{Vector{Int}, Int}}, node_id::Int)
	elements_id = find_elements_from_node_id(crack_front_edges, node_id)
	adjacent_elements = Dict{DataType, Vector{Int}}()
	for (E, el_ids) in elements_id
		if length(el_ids) == 1
			adjacent_elements[E] = el_ids
		elseif length(el_ids) == 2
			adjacent_elements[E] = el_ids[1:2]
		else
			throw("Node ID $node_id has more than two adjacent elements.") # normally impossible but just in case
		end
	end
	return adjacent_elements
end

function edges_and_adjacent_elements_from_node_id(Γ_msh::Inti.Mesh, crack_front_edges::Dict{DataType, Dict{Vector{Int}, Int}}, node_id::Int)
	edges = Dict{Vector{Int}, Inti.LagrangeElement}()
	for (E, d) in crack_front_edges
		for (edge, el_id) in d
			els = Inti.elements(Γ_msh, E)
			if node_id in edge
				el = els[el_id]
				edges[edge] = el
			end
		end
	end
	if isempty(edges)
		throw("No edges found for node ID $node_id.")
	end
	return edges
end

function lagrange_lines2local_position_from_node_id(Γ_msh::Inti.Mesh, crack_front_edges::Dict{DataType, Dict{Vector{Int}, Int}}, node_id::Int)
	edges_and_adj_elements = edges_and_adjacent_elements_from_node_id(Γ_msh, crack_front_edges, node_id)
	adj_lines = Dict{Inti.LagrangeElement{Inti.ReferenceLine}, Int}()
	for (edge, el) in edges_and_adj_elements
		edge_node_idxs = edge
		nodes = Inti.nodes(Γ_msh)[edge_node_idxs]
		nodes = SVector{length(nodes), SVector{3, Float64}}(nodes)
		line_el = Inti.LagrangeElement{Inti.ReferenceLine}(nodes)
		loc_position_1D = findfirst(x -> x == node_id, edge)
		adj_lines[line_el] = loc_position_1D
	end
	return adj_lines
end

"""
function build_front_mesh_from_edges(Γ_msh::Inti.Mesh, crack_front_edges::Dict{DataType, Dict{Vector{Int}, Int}})

	Return a 1D mesh that contains the edges of the crack front. `crack_front_edges` has to be ordered.
"""
function build_front_mesh_from_edges(Γ_msh::Inti.Mesh)::Dict
	ordered_edges = ordered_edges_path(Γ_msh)
	elements_dict = Dict{Vector{Int}, Inti.LagrangeLine}()
	for edge in ordered_edges
		nodes = Inti.nodes(Γ_msh)[edge]
		n = length(edge)
		elem = Inti.LagrangeElement{Inti.ReferenceHyperCube{1}, n, SVector{3, Float64}}(nodes)
		elements_dict[edge] = elem
	end
	return elements_dict
end

function frenet_frame_from_node_id(Γ_msh::Inti.Mesh, crack_front_edges::Dict{DataType, Dict{Vector{Int}, Int}}, node_id::Int)
	node = Inti.nodes(Γ_msh)[node_id]
	adj_lines = Dict{Inti.LagrangeElement{Inti.ReferenceLine}, Int}()
	edges_and_adj_elements = edges_and_adjacent_elements_from_node_id(Γ_msh, crack_front_edges, node_id)
	τs = Dict{Inti.LagrangeElement{Inti.ReferenceLine}, SVector{3, Float64}}()
	ns = Dict{Inti.LagrangeElement{Inti.ReferenceLine}, SVector{3, Float64}}()
	νs = Dict{Inti.LagrangeElement{Inti.ReferenceLine}, SVector{3, Float64}}()
	for (edge, el) in edges_and_adj_elements
		edge_node_idxs = edge
		el_id = crack_front_edges[typeof(el)][edge]
		node_idx = Inti.connectivity(Γ_msh, typeof(el))[:, el_id]
		nodes = Inti.nodes(Γ_msh)[edge_node_idxs]
		nodes = SVector{length(nodes), SVector{3, Float64}}(nodes)
		line_el = Inti.LagrangeElement{Inti.ReferenceLine}(nodes)
		loc_position_1D = findfirst(x -> x == node_id, edge)
		loc_position_2D = findfirst(x -> x == node_id, node_idx)
		adj_lines[line_el] = loc_position_1D
		param_nodes_1D = Inti.reference_nodes(typeof(line_el))[loc_position_1D]
		param_nodes_2D = Inti.reference_nodes(typeof(el))[loc_position_2D]
		τ = Inti.jacobian(line_el, param_nodes_1D) |> SVector{3, Float64}
		τ = τ / norm(τ)
		n = Inti.normal(el, param_nodes_2D) |> SVector{3, Float64}
		ν = τ × n
		ν = ν / norm(ν)
		τs[line_el] = τ
		ns[line_el] = n
		νs[line_el] = ν
	end
	τ = mean(collect(values(τs)))
	n = mean(collect(values(ns)))
	ν = mean(collect(values(νs)))
	return CrackBEM.Frame(node, n, ν, τ)
end

"""
function nodes(Q::Inti.TensorProductQuadrature)

	Returns the nodes of a tensor product quadrature rule as a vector of 3D vectors.
"""
function nodes(TPQ::Inti.TensorProductQuadrature{2, Tuple{Inti.GaussLegendre{N, T}, Inti.GaussLegendre{N, T}}}) where {N, T}
	qnodes = [T[getindex.(qrule.nodes, 1)...] for qrule in TPQ.quads1d]
	q1 = qnodes[1]
	q2 = qnodes[2]
	return [SVector(x, y) for x in q1, y in q2] |> vec
end

"""
	node_vals_to_quadrature(Q::Inti.Quadrature, nvals::AbstractVector)

Convert the nodal values `nvals` to quadrature node values for a given quadrature.
This function is homologous to `quadrature_to_node_vals`, but in the opposite direction.
"""
function node_vals_to_quadrature(Q::Inti.Quadrature, nvals::AbstractVector)
	msh = Q.mesh isa Inti.SubMesh ? collect(Q.mesh) : Q.mesh
	qvals = zeros(eltype(nvals), length(Q.qnodes))
	for (E, mat) in Inti.etype2mat(msh)
		L = Inti.lagrange_basis(E)
		qrule = Q.etype2qrule[E]
		nel = size(mat, 2)
		for n in 1:nel
			node_tags = mat[:, n]
			qtags = Q.etype2qtags[E][:, n]
			local_nvals = nvals[node_tags]
			q̂ = nodes(qrule)
			qvals[qtags] .= map(q -> reduce(+, L(q) .* local_nvals), q̂)
		end
	end
	return qvals
end

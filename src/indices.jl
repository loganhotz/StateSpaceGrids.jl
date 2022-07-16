# file for indexing into AbstractMesh objects



Base.getindex(G::AbstractMesh, I...) = G[I]


Base.getindex(G::AbstractMesh{T, 1}, I::Integer) where {T} = G.dims[1][I]
Base.getindex(G::AbstractMesh{T, 1}, I...) where {T}       = throw(BoundsError(G, I))


Base.getindex(G::AbstractMesh, I::Integer) = (cartesian = Tuple(G.R[I]); G[cartesian])
Base.getindex(G::AbstractMesh, I::Tuple)   = throw(BoundsError(G, I))
function Base.getindex(G::AbstractMesh{T, N}, I::NTuple{N, <:Integer}) where {T, N}
	g = Vector{T}(undef, N)
	for (i, (d, c)) in enumerate(zip(G.dims, I))
		g[i] = d[c]
	end
	return g
end

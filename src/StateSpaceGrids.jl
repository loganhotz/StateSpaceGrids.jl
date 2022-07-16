module StateSpaceGrids

using Interpolations

export AbstractMesh
export AbstractGrid
export Grid

export AbstractStateSpaceDimension
export AbstractMeshDimension
export AbstractGridDimension

export MeshDimension
export GridDimension

export allocate_solution



abstract type AbstractMesh{T, N} end
abstract type AbstractGrid{T, N} <: AbstractMesh{T, N} end

abstract type AbstractStateSpaceDimension{T<:Number} end
abstract type AbstractMeshDimension{T} <: AbstractStateSpaceDimension{T} end
abstract type AbstractGridDimension{T} <: AbstractStateSpaceDimension{T} end



include("indices.jl")
include("solutions.jl")



struct GridDimension{T} <: AbstractGridDimension{T}
	a::T
	b::T
	n::Int
	points::AbstractVector{T}

	function GridDimension(a, b, n::Int, p::AbstractVector)
		if b <= a
			error("lower bound must be below upper bound")
		end

		uniformΔ = allapprox(diff(p))
		!uniformΔ && error("GridDimension requires uniform step size")

		a, b, _ = promote(a, b, first(p))
		return new{eltype(p)}(a, b, n, p)
	end
end
function GridDimension(a::Number, b::Number, n::Int)
	return GridDimension(a, b, n, collect(range(a, b, n)))
end
GridDimension(a::Number, b::Number) = GridDimension(a, b, 11)
GridDimension() = GridDimension(0.0, 1.0)

function GridDimension(A::AbstractVector{<:Number})
	a, b = minimum(A), maximum(A)
	return GridDimension(a, b, length(A), A)
end
GridDimension{T}(A::AbstractVector{<:Number}) where T = GridDimension(Vector{T}(A))



function Base.getindex(D::AbstractStateSpaceDimension, i::Int)
	1 <= i <= D.n || throw(BoundsError(D, i))
	return D.points[i]
end
Base.eltype(D::AbstractStateSpaceDimension{T}) where {T} = T
Base.length(D::AbstractStateSpaceDimension)              = D.n
Base.iterate(D::AbstractStateSpaceDimension, s=1)        = s > D.n ? nothing : (D[s], s+1)

function Base.show(io::IO, D::AbstractGridDimension)
	T = typeof(D)
	a, b, n = D.a, D.b, D.n
	print(io, "$T($a, $b, $n)")
end



# conversion rules
function Base.convert(
    ::Type{GridDimension{T}},
    D::GridDimension{S}
) where {T, S}
    Q = promote_type(T, S)
    return GridDimension{Q}(D.points)
end
Base.convert(::Type{GridDimension}, A::AbstractVector{<:Number}) = GridDimension(A)



struct Grid{T, N} <: AbstractGrid{T, N}
	dims::NTuple{N, GridDimension{T}}
	R::CartesianIndices{N}

    function Grid(dims::NTuple{N, AbstractGridDimension}) where {N}
		lengths = Tuple(1:length(d) for d in dims)
        T       = promote_eltype(dims...)

		return new{T, N}(dims, CartesianIndices(lengths))
	end

	function Grid(dims::NTuple{N, GridDimension{T}}) where {T, N}
		lengths = Tuple(1:length(d) for d in dims)

		return new{T, N}(dims, CartesianIndices(lengths))
	end
end
Grid(D::AbstractGridDimension) = Grid(tuple(D))
Grid(args...)                  = Grid(args)
Grid()                         = Grid(tuple(GridDimension()))



# mostly replicating the basic functions in the manual:
#	https://docs.julialang.org/en/v1/manual/arrays/#man-array-indexing
Base.size(G::AbstractMesh)                      = Tuple(length(d) for d in G.dims)
Base.size(G::AbstractMesh, n::Integer)          = length(G.dims[n])
Base.axes(G::AbstractMesh)                      = axes(G.R)
Base.axes(G::AbstractMesh, n::Integer)          = axes(G.R, n)
Base.ndims(G::AbstractMesh{T, N}) where {T, N}  = N
Base.eltype(G::AbstractMesh{T, N}) where {T, N} = T
Base.length(G::AbstractMesh)                    = prod(size(G))
Base.iterate(G::AbstractMesh, s=1)              = s > length(G) ? nothing : (G[s], s+1)
Base.eachindex(G::AbstractMesh)                 = eachindex(G.R) # Arrays use linear...

# implements joint iteration over a Mesh and its indices, a la `enumerate`, via the `pairs`
# function. using `enumerate(G::AbstractMesh)` returns linear indices, while iterating over
# the `pairs(G::AbstractMesh)` iterator returns cartesian indices
Base.keys(G::AbstractMesh) = G.R


function Base.show(io::IO, G::AbstractGrid)
	ndims(G) > 1 ? dims = join(size(G), "x") : dims = "1-dim"

	print(io, "$dims $(typeof(G)):")
	for d in G.dims
		print(io, "\n    $d")
	end
end



# used when checking grid steps are equal
@inline function allapprox(x)
	length(x) < 2 && return true
	x1 = x[1]
	i  = 2

	@inbounds for i = 2:length(x)
		isapprox(x[i], x1) || return false
	end
	return true
end

# promotion & conversion of dimensions and meshes
eltypeof(x) = typeof(x)
eltypeof(x::AbstractStateSpaceDimension) = eltype(x)
eltypeof(x::AbstractMesh)                = eltype(x)

Bottom = Union{}

promote_eltypeof()          = Bottom
promote_eltypeof(x1, xs...) = promote_type(eltypeof(x1), promote_eltypeof(xs...))

promote_eltype()          = Bottom
promote_eltype(x1, xs...) = promote_type(eltype(x1), promote_eltype(xs...))



end # module

# concrete implementations of meshes and grids



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
#    https://docs.julialang.org/en/v1/manual/arrays/#man-array-indexing
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



# promotion & conversion of dimensions and meshes. see `src/grid.jl` for other `eltypeof`
#   methods
eltypeof(x::AbstractMesh) = eltype(x)

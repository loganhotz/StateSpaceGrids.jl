# concrete implementations of meshes and grids



# interfaces for all Meshes. required for subtypes of AbstractArray
Base.size(M::AbstractMesh)                 = Tuple(length(d) for d in M.D)
Base.getindex(M::AbstractMesh, i::Integer) = getindex(M, M.R[i])
function Base.getindex(M::AbstractMesh{T, N}, I::Vararg{Integer, N}) where {T, N}
    g = Vector{T}(undef, N)
    for (i, (d, cartesian)) in enumerate(zip(M.D, I))
        g[i] = d[cartesian]
    end
    return g
end

# Meshes and Grids, in contrast to their supertype AbstractArray, produce `Vector`s of
#   numbers when indexed into
Base.eltype(M::AbstractMesh{T, N}) where {T, N}  = Vector{T}
Base.valtype(M::AbstractMesh{T, N}) where {T, N} = Vector{T}
Base.keytype(M::AbstractMesh{T, N}) where {T, N} = CartesianIndex{N}
Base.IndexStyle(M::AbstractMesh)                 = IndexCartesian()


# not entirely sure why I have to do this... without this method, and specifically without
#   the `keys(M)` in the second argument, calls to `pairs(M::AbstractMesh)` generate error
# ```julia-repl
# ERROR: MethodError: Cannot `convert` an object of type Vector[T} to an object of type T
# ```
#   i would have thought the redefined `eltype` and `valtype` methods above would have
#   obviated this method
function Base.Pairs(M::AbstractMesh, itr::I) where {I}
    return Base.Pairs{eltype(I), eltype(M), I, typeof(M)}(M, keys(M))
end



# printing to REPL
function Base.show(io::IO, G::AbstractGrid)
    ndimens(G) > 1 ? dimens = join(size(G), "x") : dimens = "1-dim"

    print(io, "$dimens $(typeof(G)):")
    for d in dims(G)
        print(io, "\n    $d")
    end
end



"""
    Grid{T<:Real, N}

an `N`-dimensional grid whose constituent elements are `GridDimension`s with element
type `T`. each dimension must have a constant step size between elements, but the step
size of different dimensions has no constraints
"""
struct Grid{T, N} <: AbstractGrid{T, N}
    D::NTuple{N, GridDimension{T}}
    R::CartesianIndices{N}
end



"""
    Grid(D::NTuple{N, AbstractVector})
    Grid{T}(D::NTuple{N, AbstractVector})

creates an `N`-dimensional `Grid` of type `T`. if `T` is not provided, it is inferred
by the element types of the vectors in `D`
"""
function Grid(D::NTuple{N, AbstractVector}) where {N}
    lengths = Tuple(1:length(d) for d in D)
    T       = Base.promote_eltype(D...)

    return Grid{T, N}(D, CartesianIndices(lengths))
end
function Grid{T}(D::NTuple{N, AbstractVector}) where {T, N}
    lengths = Tuple(1:length(d) for d in D)
    return Grid{T, N}(D, CartesianIndices(lengths))
end

Grid(D::AbstractVector)              = Grid((D, ))
Grid{T}(D::AbstractVector) where {T} = Grid((D, ))


"""
    grid(D...)
    grid(D::NTuple{N, AbstractVector})
    grid(::Type{T<:Real}, D...)
    grid(::Type{T<:Real}, D::NTuple{N, AbstractVector})

creates an `N`-dimensional `Grid` of type `T`. if `T` is not provided, it is inferred
by the element types of the vectors in `D`
"""
grid(D...)                                                 = Grid(D)
grid(::Type{T}, D...) where {T}                            = Grid{T}(D)

grid(D::NTuple{N, AbstractVector}) where {N}               = Grid(D)
grid(::Type{T}, D::NTuple{N, AbstractVector}) where {T, N} = Grid{T}(D)

"""
    grid(D::AbstractVector, n::Integer)
    grid(::Type{T<:Real}, D::AbstractVector, n::Integer)

creates a `Grid` whose `n` constituent dimensions are replicas of `D`
"""
grid(D::AbstractVector, n::Integer)                      = Grid(Tuple(D for _ in 1:n))
grid(::Type{T}, D::AbstractVector, n::Integer) where {T} = Grid{T}(Tuple(D for _ in 1:n))



"""
    Mesh{T<:Real, N}

an `N`-dimensional mesh whose constituent elements are `MeshDimension`s with element
type `T`.
"""
struct Mesh{T, N} <: AbstractMesh{T, N}
    D::NTuple{N, MeshDimension{T}}
    R::CartesianIndices{N}
end

"""
    Mesh(D::NTuple{N, AbstractVector})
    Mesh{T}(D::NTuple{N, AbstractVector})

creates an `N`-dimensional `Mesh` of type `T`. if `T` is not provided, it is inferred by
the element types of the vectors in `D`
"""
function Mesh(D::NTuple{N, AbstractVector}) where {N}
    lengths = Tuple(1:length(d) for d in D)
    T       = Base.promote_eltype(D...)

    return Mesh{T, N}(D, CartesianIndices(lengths))
end
function Mesh{T}(D::NTuple{N, AbstractVector}) where {T, N}
    lengths = Tuple(1:length(d) for d in D)
    return Mesh{T, N}(D, CartesianIndices(lengths))
end

Mesh(D::AbstractVector)              = Mesh((D, ))
Mesh{T}(D::AbstractVector) where {T} = Mesh((D, ))


"""
    mesh(D...)
    mesh(D::NTuple{N, AbstractVector})
    mesh(::Type{T<:Real}, D...)
    mesh(::Type{T<:Real}, D::NTuple{N, AbstractVector})

creates an `N`-dimensional `Mesh` of type `T`. if `T` is not provided, it is inferred
by the element types of the vectors in `D`
"""
mesh(D...)                                                 = Mesh(D)
mesh(::Type{T}, D...) where {T}                            = Mesh{T}(D)

mesh(D::NTuple{N, AbstractVector}) where {N}               = Mesh(D)
mesh(::Type{T}, D::NTuple{N, AbstractVector}) where {T, N} = Mesh{T}(D)

"""
    mesh(D::AbstractVector, n::Integer)
    mesh(::Type{T<:Real}, D::AbstractVector, n::Integer)

creates a `Mesh` whose `n` constituent dimensions are replicas of `D`
"""
mesh(D::AbstractVector, n::Integer)                      = Mesh(Tuple(D for _ in 1:n))
mesh(::Type{T}, D::AbstractVector, n::Integer) where {T} = Mesh{T}(Tuple(D for _ in 1:n))



# interface for both grids and meshes
dims(M::AbstractMesh) = M.D



# interface for AbstractGrid objects
steps(G::AbstractGrid) = Tuple(step(d) for d in dims(G))

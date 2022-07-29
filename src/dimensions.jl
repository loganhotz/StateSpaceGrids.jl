# concrete implementations of mesh and grid dimensions



# interfaces shared between grid and mesh dimensions. required for AbstractVector subtypes
Base.size(D::AbstractStateSpaceDimension) = (D.n, )
Base.@propagate_inbounds Base.getindex(D::AbstractStateSpaceDimension, i::Int) = D.points[i]

# interfaces; not required
Base.eltype(D::AbstractStateSpaceDimension{T}) where {T} = T
Base.IndexStyle(D::AbstractStateSpaceDimension)          = IndexLinear()

Base.similar(D::AbstractStateSpaceDimension) = typeof(D)(D.points)

Base.step(D::AbstractGridDimension) = D.step


# printing to REPL
function Base.show(io::IO, D::AbstractStateSpaceDimension)
    T       = typeof(D)
    a, b, n = minimum(D.points), maximum(D.points), D.n
    print(io, "$T($a, $b, points=$n)")
end
function Base.show(io::IO, D::AbstractGridDimension)
    T = typeof(D)
    a, b, n, Δ = minimum(D.points), maximum(D.points), D.n, D.step
    print(io, "$T($a, $b, points=$n, step = $Δ)")
end



"""
    GridDimension{T<:Real}

a monotonically increasing vector of points of type `T`. the defining aspect of a
`GridDimension` is that the difference between consecutive points must be constant along
ts entire length
"""
struct GridDimension{T} <: AbstractGridDimension{T}
    points::Vector{T}
    n::Int64
    step::T

    function GridDimension(points::AbstractVector{<:Real})

        validate_grid_points(points)

        step = points[2] - points[1]
        n    = length(points)

        return new{eltype(points)}(points, n, step)

    end
end

"""
    GridDimension(points::AbstractVector{<:Real})

create a `GridDimension` of the elements of `points`. checks are made as to whether those
elements are equidistant from their neighbors
"""
GridDimension{T}(A::AbstractVector{<:Real}) where {T<:Real} = GridDimension(Vector{T}(A))
# the definition above allows for e.g. promotions of A::Vector{Int64} to a Dimension whose
#   desired eltype is Float64

"""
    GridDimension(start::Real, stop::Real, n::Int)
    GridDimension(start::Real, stop::Real; n::Int, step::Real)
    GridDimension(start::Real; stop::Real, n::Int, step::Real)
    GridDimension(; start::Real, stop::Real, n::Int, step::Real)

create a `GridDimension` uniformly spaced points between `a` and `b`, inclusive. the points
in the Dimension are uniquely determined by any three of `start`, `stop`, `step`, and `n`.
valid invocations of `GridDimension` are:
* call `GridDimension` with any three of `start`, `stop`, `step`, or `n`
* call `GridDimension` with any two of `start`, `stop`, or `n`. in this case, `step` will
    be assumed
"""
GridDimension(A::AbstractRange) = GridDimension(collect(A))
GridDimension(start; stop=nothing, n::Union{Integer, Nothing}=nothing, step=nothing) =
    GridDimension(range(start=start, stop=stop, length=n, step=step))
GridDimension(start, stop; n::Union{Integer, Nothing}=nothing, step=nothing) =
    GridDimension(range(start=start, stop=stop, length=n, step=step))
GridDimension(start, stop, n::Integer) =
    GridDimension(range(start=start, stop=stop, length=n))
GridDimension(;start=nothing, stop=nothing, n::Union{Integer, Nothing}=nothing, step=nothing) =
    GridDimension(range(start=start, stop=stop, length=n, step=step))



"""
    MeshDimension{T<:Real}

a monotonically increasing vector of state space points of type `T`. as opposed to the
[`GridDimension`](@ref), the difference between consecutive elements need not be equal
"""
struct MeshDimension{T} <: AbstractMeshDimension{T}
    points::Vector{T}
    n::Int64

    function MeshDimension(points::AbstractVector{<:Real})

        validate_mesh_points(points)
        n = length(points)

        return new{eltype(points)}(points, n)
    end
end

"""
    MeshDimension(points::AbstractVector{<:Real})

create a `MeshDimension` of the elements of `points`. checks are made as to whether those
elements are strictly increasing
"""
MeshDimension{T}(A::AbstractVector{<:Real}) where {T<:Real} = MeshDimension(Vector{T}(A))
# the definition above allows for e.g. promotions of A::Vector{Int64} to a Dimension whose
#   desired eltype is Float64



# conversion rules
for D = (:GridDimension, :MeshDimension)
    @eval Base.convert(::Type{$D{T}}, A::AbstractVector{<:Real}) where {T} = $D{T}(A)
    @eval Base.convert(::Type{$D}, A::AbstractVector{<:Real})              = $D(A)
end



# constructor validations
function validate_grid_points(x::AbstractVector)
    uniformΔ = allapprox(diff(x))
    !uniformΔ && error("subtypes of AbstractGridDimension require a uniform step size")

    a, b = x[1], x[end]
    (a >= b) && error("state space dimensions must be increasing")

    return true
end

@inline function allapprox(x::AbstractVector)
    length(x) < 2 && return true
    x1 = x[1]
    i  = 2

    @inbounds for i = 2:length(x)
        isapprox(x[i], x1) || return false
    end
    return true
end


function validate_mesh_points(x::AbstractVector)
    dx = diff(x)
    any(dx .<= 0) && error("state space dimensions must be increasing")

    return true
end

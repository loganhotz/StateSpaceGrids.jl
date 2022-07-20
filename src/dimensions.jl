# concrete implementations of mesh and grid dimensions



# interfaces shared between grid and mesh dimensions
function Base.getindex(D::AbstractStateSpaceDimension, i::Int)
    1 <= i <= D.n || throw(BoundsError(D, i))
    return D.points[i]
end
Base.eltype(D::AbstractStateSpaceDimension{T}) where {T} = T
Base.length(D::AbstractStateSpaceDimension)              = D.n
Base.iterate(D::AbstractStateSpaceDimension, s=1)        = s > D.n ? nothing : (D[s], s+1)



"""
    GridDimension{T}

a monotonically increasing/decreasing vector of points of type `T`. the defining aspect
of a `GridDimension` is that the difference between consecutive points must be constant
along its entire length
"""
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

"""
    GridDimension(a::Real, b::Real, n::Int)

create a `GridDimension` of `n` uniformly spaced points between `a` and `b`, inclusive
"""
function GridDimension(a::Real, b::Real, n::Int)
    return GridDimension(a, b, n, collect(range(a, b, n)))
end

"""
    GridDimension(A::AbstractVector{<:Real})

create a `GridDimension` of the elements of `A`. checks are made as to whether those
elements are equidistant from their neighbors
"""
function GridDimension(A::AbstractVector{<:Real})
    a, b = minimum(A), maximum(A)
    return GridDimension(a, b, length(A), A)
end
GridDimension{T}(A::AbstractVector{<:Real}) where {T<:Real} = GridDimension(Vector{T}(A))
# second case above is to allow for, e.g. promotions of A::Vector{Int64} to a dimension
#   whose desired eltyep is Float64



# printing to REPL
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

# promotion & conversion of dimensions
eltypeof(x) = typeof(x)
eltypeof(x::AbstractStateSpaceDimension) = eltype(x)

Bottom = Union{}

promote_eltypeof()          = Bottom
promote_eltypeof(x1, xs...) = promote_type(eltypeof(x1), promote_eltypeof(xs...))

promote_eltype()          = Bottom
promote_eltype(x1, xs...) = promote_type(eltype(x1), promote_eltype(xs...))

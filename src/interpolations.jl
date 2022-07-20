# interpolating over state space meshes
#   this file makes use of EllipsisNotation's `..` type for gready gobbling of axes



"""
    interpolate(
        M::AbstractMesh,
        A::AbstractArray,
        I::Interpolations.InterpolationType
        E::Interpolations.BoundaryCondition
    )

constructs an interpolating function over the domain of the given Mesh object. the
interpolation mode must be one of `NoInterp`, `Gridded(Constant())`, or
`Gridded(Linear())`, as per the documentation in the `Interpolations` package. similarly,
the extrapolation mode can be set as `Throw`, `Flat`, `Line`, `Periodic`, or `Reflect`
"""
function Interpolations.interpolate(
    M::AbstractMesh,
    A::AbstractArray,
    I::IT,
    E::ET
) where {IT <: Interpolations.InterpolationType, ET <: Interpolations.BoundaryCondition}

    dims = Tuple(D.points for D in M.dims)
    n    = size(A, 1)

    V = Vector{Any}(nothing, n)
    for dec_var = 1:n
        itp = interpolate(dims, A[dec_var, ..], I)
        ext = extrapolate(itp, E)
        V[dec_var] = ext
    end

    vitp(args...; kwargs...) = [itp(args...; kwargs...) for itp in V]
    return vitp
end



"""
    interpolate(
        M::AbstractMesh,
        A::AbstractArray,
        I::Interpolations.InterpolationType
    )

constructs an interpolating function over the domain of the given Mesh. the interpolation
method is set via the `I` parameter, and should be `NoInterp`, `Gridded(Constant())`, or
`Gridded(Linear())`, just as in the `Interpolations` package. extrapolation is assumed to
be linear outside the interpolating region
"""
function Interpolations.interpolate(
    M::AbstractMesh,
    A::AbstractArray,
    I::IT
) where {IT <: Interpolations.InterpolationType}
    return interpolate(M, A, I, Line())
end



"""
    interpolate(
        M::AbstractMesh,
        A::AbstractArray,
        E::Interpolations.BoundaryCondition
    )

constructs an interpolation function over the domain of the given Mesh. the extrapolation
method is determined by the `E` parameter, and must be one of `Throw`, `Flat`, `Line`,
`Periodic`, or `Reflect`. the interpolation is uses `Gridded(Linear())`
"""
function Interpolations.interpolate(
    M::AbstractMesh,
    A::AbstractArray,
    E::ET
) where {ET <: Interpolations.BoundaryCondition}
    return interpolate(M, A, Gridded(Linear()), E)
end



"""
    interpolate(
        M::AbstractMesh,
        A::AbstractArray
    )

constructs an interpolation function over the domain of the given Mesh. the method of
interpolation is `Gridded(Linear())`, and that of extrapolation is just `Line()`
"""
function Interpolations.interpolate(
    M::AbstractMesh,
    A::AbstractArray
)
    return interpolate(M, A, Gridded(Linear()), Line())
end

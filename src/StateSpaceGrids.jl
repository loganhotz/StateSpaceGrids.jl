module StateSpaceGrids



using EllipsisNotation
using Interpolations
using Reexport

# exposes all the inter- extrapolation structures
@reexport using Interpolations



# abstract grids
export AbstractMesh
export AbstractGrid

# concrete grids
export Grid

# abstract dimensions
export AbstractStateSpaceDimension
export AbstractMeshDimension
export AbstractGridDimension

# concrete dimensions
export GridDimension

# creating state-space solutions, and an `Interpolations.jl` function with methods that
#   interact with AbstractMesh structs
export allocate_solution
export interpolate



abstract type AbstractMesh{T, N} end
abstract type AbstractGrid{T, N} <: AbstractMesh{T, N} end

abstract type AbstractStateSpaceDimension{T<:Number} end
abstract type AbstractMeshDimension{T} <: AbstractStateSpaceDimension{T} end
abstract type AbstractGridDimension{T} <: AbstractStateSpaceDimension{T} end


include("dimensions.jl")
include("grids.jl")
include("indices.jl")
include("interpolations.jl")
include("solutions.jl")



end # module

module StateSpaceGrids



using EllipsisNotation
using Interpolations
using Reexport

# exposes all the inter- extrapolation structures
@reexport using Interpolations



# abstract grids and meshes
export AbstractMesh
export AbstractGrid

# concrete grids and meshes
export Grid
export Mesh

# abstract dimensions
export AbstractStateSpaceDimension
export AbstractMeshDimension
export AbstractGridDimension

# concrete dimensions
export GridDimension
export MeshDimension

# instantiation functions
export grid
export mesh

# creating state-space solutions, and an `Interpolations.jl` function with methods that
#   interact with AbstractMesh structs
export allocate_solution
export interpolate

# utils for meshes, grids, and dimensions
export dims
export steps



abstract type AbstractMesh{T<:Real, N} <: AbstractArray{T, N} end
abstract type AbstractGrid{T, N}       <: AbstractMesh{T, N} end

abstract type AbstractStateSpaceDimension{T<:Real} <: AbstractVector{T} end
abstract type AbstractMeshDimension{T} <: AbstractStateSpaceDimension{T} end
abstract type AbstractGridDimension{T} <: AbstractStateSpaceDimension{T} end



include("dimensions.jl")
include("grids.jl")
include("interpolations.jl")
include("solutions.jl")



end # module

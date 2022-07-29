# file for generating solution arrays


"""
    allocate_solution(G::AbstractMesh, P::Integer, [init=1])

utility function for initializing an array at which to write the policy rule solutions

# Arguments
`G::AbstractMesh`: a state space mesh
`P::Integer`: the number of policy rules in the solution
`init::Real`: a scalar for every entry of the array
"""
function allocate_solution(
    G::AbstractMesh{T, N},
    P::Integer,
    init::Real
) where {T, N}
    A = Array{T, N+1}(undef, P, size(G)...)
    fill!(A, init)
    return A
end
allocate_solution(G::AbstractMesh, P::Integer) = allocate_solution(G, P, 1.0)

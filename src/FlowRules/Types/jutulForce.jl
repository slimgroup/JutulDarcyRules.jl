export jutulForce

struct jutulForce{D, T}
    irate::T
    name::Vector{Symbol}
    loc::Vector{NTuple{D, T}}
end

jutulForce(irate::T, loc::NTuple{D, T}) where {D, T} = jutulForce(irate, [loc])
jutulForce(irate::T, loc::Vector{NTuple{D, T}}) where {D, T}= jutulForce(irate, vcat(:Injector, [:Producer for i = 1:length(loc)]), loc)

# display(f::jutulForce{D, T}) where {D, T} =
#     println("$(D)D jutulForce structure with $(length(f.loc)) injection/production wells and rate $(f.irate) m^3/s")

==(A::jutulForce{D, T}, B::jutulForce{D, T}) where {D,T} = (A.irate == B.irate && A.name == B.name && A.loc == B.loc)
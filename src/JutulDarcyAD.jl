__precompile__()
module JutulDarcyAD

    using LinearAlgebra
    using Jutul
    using ChainRulesCore
    import Jutul: JutulGeometry, get_facepos, compute_face_trans, compute_half_face_trans, expand_perm

    include("PropertyConversion/PropertyConversion.jl")
    include("FlowAD/FlowAD.jl")

end # module JutulDarcyAD

__precompile__()
module JutulDarcyAD

    using LinearAlgebra
    using Statistics
    using Jutul
    using Optim
    using Flux
    using ChainRulesCore
    import Jutul: JutulGeometry, get_facepos, compute_face_trans, compute_half_face_trans, expand_perm
    import Jutul: SimulationModel

    const Darcy = 9.869232667160130e-13
    const md = Darcy * 1e-3
    
    include("PropertyConversion/PropertyConversion.jl")
    include("FlowAD/FlowAD.jl")

end # module JutulDarcyAD

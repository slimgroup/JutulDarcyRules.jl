__precompile__()
module JutulDarcyAD

    export JutulDarcyADPATH, Darcy, md

    using LinearAlgebra
    using Statistics
    using JutulDarcy
    using Jutul
    using Optim
    using Flux
    using ChainRulesCore
    import Jutul: JutulGeometry, get_facepos, compute_face_trans, compute_half_face_trans, expand_perm
    import Jutul: SimulationModel
    import Base: display

    JutulDarcyADPATH = dirname(pathof(JutulDarcyAD))

    const Darcy = 9.869232667160130e-13
    const md = Darcy * 1e-3
    
    const sys = ImmiscibleSystem((LiquidPhase(), VaporPhase()))
    const day = 24*3600.0

    include("PropertyConversion/PropertyConversion.jl")
    include("FlowAD/FlowAD.jl")

end # module JutulDarcyAD

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
    import Base: +, -, *, /, ==
    import Base: display, length, size, getindex, setindex!, IndexStyle, vec, firstindex, lastindex
    import LinearAlgebra: norm, dot
    import ChainRulesCore: rrule

    JutulDarcyADPATH = dirname(pathof(JutulDarcyAD))

    visCO2 = 1e-4
    visH2O = 1e-3
    ρCO2 = 501.9
    ρH2O = 1053.0

    const Darcy = 9.869232667160130e-13
    const md = Darcy * 1e-3
    
    const sys = ImmiscibleSystem((VaporPhase(), AqueousPhase()))
    const day = 24*3600.0

    include("PropertyConversion/PropertyConversion.jl")
    include("FlowAD/FlowAD.jl")

    ### hack
    import Jutul: solve!, linear_solve_return
    function solve!(sys; dx = sys.dx, r = sys.r, J = sys.jac)
        dx .= -(J\r)
        @assert all(isfinite, dx) "Linear solve resulted in non-finite values."
        return linear_solve_return()
    end

end # module JutulDarcyAD

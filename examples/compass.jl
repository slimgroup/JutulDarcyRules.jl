## A 64×64 2D example for permeability inversion of a tortuous channel

using DrWatson
@quickactivate "JutulDarcyAD-example"

using JutulDarcyAD
using LinearAlgebra
using PyPlot
using Flux
using LineSearches
using JLD2
using JUDI

sim_name = "2D-K-inv"
exp_name = "compass"

mkpath(datadir())
mkpath(plotsdir())

## grid size
JLD2.@load datadir("BGCompass_tti_625m.jld2") m d;
v = Float64.(sqrt.(1f0./m));
d = Float64.(d);
function VtoK(v::Matrix{T}, d::Tuple{T, T}; α=T(1.03f0)) where T

    n = size(v)
    idx_wb = find_water_bottom(v.-minimum(v))
    idx_ucfmt = find_water_bottom((v.-T(3.5)).*(v.>T(3.5)))
    Kh = zeros(T, n)
    capgrid = Int(round(T(50)/d[2]))
    for i = 1:n[1]
        Kh[i,1:idx_wb[i]-1] .= T(1e-10)  # water layer
        Kh[i,idx_wb[i]:idx_ucfmt[i]-capgrid-1] .= α*v[i,idx_wb[i]:idx_ucfmt[i]-capgrid-1].+T(15)
        Kh[i,idx_ucfmt[i]-capgrid:idx_ucfmt[i]-1] .= T(1e-3)
        Kh[i,idx_ucfmt[i]:end] .= α*v[i,idx_ucfmt[i]:end].+T(200)
    end
    return Kh
end
Kh = VtoK(v, d);
K = Float64.(Kh * md);
n = (size(K,1), 1, size(K,2))
d = (d[1], 10.0, d[2])

ϕ = 0.25
model = jutulModel(n, d, ϕ, K1to3(K))

## simulation time steppings
tstep = 40 * ones(50)
tot_time = sum(tstep)

## injection & production
inj_loc = (3, 1, n[end]) .* d
prod_loc = (n[1]-3, 1, n[end]) .* d
irate = 5e-3
q = jutulForce(irate, [inj_loc, prod_loc])

## set up modeling operator
S = jutulModeling(model, tstep)

## simulation
mesh = CartesianMesh(model)
T(x) = log.(KtoTrans(mesh, K1to3(exp.(x))))

logK = log.(K)

@time state = S(T(logK), q)

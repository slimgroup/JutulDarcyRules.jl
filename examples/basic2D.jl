## A simple 2D example for fluid-flow simulation

using JutulDarcyAD
using LinearAlgebra
using PyPlot
using SlimPlotting

## grid size
n = (30, 1, 15)
d = (30.0, 30.0, 30.0)

## permeability
K = 20 * md * ones(n)
ϕ = 0.25
model = jutulModel(n, d, ϕ, K1to3(K))

## simulation time steppings
tstep = 20 * ones(50)
tot_time = sum(tstep)

## injection & production
inj_loc = (3, 1, 9) .* d
prod_loc = (28, 1, 9) .* d
irate = 5e-3
q = jutulForce(irate, [inj_loc, prod_loc])

## set up modeling operator
S = jutulModeling(model, tstep)

## simulation
Trans = KtoTrans(CartesianMesh(model), K1to3(K))
@time result = S(log.(Trans), q)
@time result = S(log.(Trans), q)

## plotting
fig=figure(figsize=(20,12));
subplot(1,2,1);
imshow(reshape(result.states[end].Saturations, n[1], n[end])', vmin=0, vmax=1); colorbar(); title("saturation")
subplot(1,2,2);
imshow(reshape(result.states[end].Pressure, n[1], n[end])', vmin=0, vmax=maximum(result.states[end].Pressure)); colorbar(); title("pressure")

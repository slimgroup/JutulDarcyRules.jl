function test_config()
    n = (30, 1, 15)
    d = (30.0, 30.0, 30.0)

    ## permeability
    K0 = 200 * md * ones(n)
    ϕ = 0.25
    K = deepcopy(K0)
    K[:,:,1:2:end] .*= 100

    model0 = jutulModel(n, d, ϕ, K1to3(K0))
    model = jutulModel(n, d, ϕ, K1to3(K))

    ## simulation time steppings
    tstep = 50 * ones(10)

    ## injection & production
    inj_loc = (15, 1, 10) .* d
    prod_loc = (30, 1, 10) .* d
    irate = 5e-3
    q = jutulForce(irate, inj_loc)
    q1 = jutulSource(irate, [inj_loc, prod_loc])
    q2 = jutulVWell(irate, inj_loc[1:2]; startz = 9 * d[3], endz = 11 * d[3])
    state0 = jutulState(JutulDarcyAD.setup_well_model(model, q, tstep)[3])
    state1 = jutulSimpleState(model)
    return model, model0, q, q1, q2, state0, state1, tstep
end

mean(x) = sum(x)/length(x)

## adapted from https://github.com/slimgroup/JUDI.jl/blob/master/test/testing_utils.jl
function grad_test(misfit, x0, dx, g; maxiter=6, h0=5f-2, stol=1f-1)
    # init
    err1 = zeros(Float32, maxiter)
    err2 = zeros(Float32, maxiter)
    
    gdx = dot(g, dx)
    f0 = misfit(x0)
    h = h0

    @printf("%11.5s, %11.5s, %11.5s, %11.5s, %11.5s, %11.5s \n", "h", "gdx", "e1", "e2", "rate1", "rate2")
    for j=1:maxiter
        f = misfit(x0 + h*dx)
        err1[j] = norm(f - f0, 1)
        err2[j] = norm(f - f0 - h*gdx, 1)
        j == 1 ? prev = 1 : prev = j - 1
        @printf("%5.5e, %5.5e, %5.5e, %5.5e, %5.5e, %5.5e \n", h, h*norm(gdx, 1), err1[j], err2[j], err1[prev]/err1[j], err2[prev]/err2[j])
        h = h * .8f0
    end

    rate1 = err1[1:end-1]./err1[2:end]
    rate2 = err2[1:end-1]./err2[2:end]
    @test isapprox(mean(rate1), 1.25f0; atol=stol)
    @test isapprox(mean(rate2), 1.5625f0; atol=stol)
end

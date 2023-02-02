function setup_model_state()

    ## grid size
    n = (rand(3:5), rand(3:5), rand(3:5))
    d = 10 .* (rand(), rand(), rand())

    ## permeability
    K = 200 * rand() * md * ones(n)
    ϕ = rand()
    model = jutulModel(n, d, ϕ, K1to3(K))
    init_state = jutulState(model)
    return model, init_state

end

function test_config()
    n = (10, 1, 10)
    d = (100.0, 100.0, 100.0)

    K0 = 20 * md * ones(n)
    K = deepcopy(K0)
    K[:,:,6:end] .*= 5
    ϕ = 0.25

    model0 = jutulModel(n, d, ϕ, K1to3(K0))
    model = jutulModel(n, d, ϕ, K1to3(K))

    ## simulation time steppings
    tstep = 20 * ones(50)

    ## injection & production
    inj_loc = (2, 1, 8) .* d
    prod_loc = (8, 1, 8) .* d
    irate = 5e-3
    q = jutulForce(irate, [inj_loc, prod_loc])

    return model, model0, q, tstep
end

mean(x) = sum(x)/length(x)

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

model, model0, q, q1, q2, init_state, init_state_simple, tstep = test_config();

@testset "Test jutulState" begin
    @info "length, size"
    @test size(init_state) == (length(init_state),)
    @test length(init_state) == length(vec(init_state))
    
    @info "getindex, setindex!"
    rdm_idx = randperm(length(init_state))[1:5]
    rdm_num = rand(5)
    init_state[rdm_idx] = rdm_num
    @test init_state[rdm_idx] == rdm_num
    @test vec(init_state) == init_state[1:end]
    
    @info "test =="
    init_state1 = deepcopy(init_state)
    init_state1[1] = 0.2
    @test init_state1 != init_state
    init_state[1] = 0.2
    @test init_state1 == init_state

    @info "test dict"
    dict(init_state) == init_state.state

    @info "test translation"
    half_init_state = init_state(vec(init_state) ./ 2)
    JutulDarcyRules.check_valid_state(half_init_state)
    @test vec(half_init_state) == 0.5 * vec(init_state)
end
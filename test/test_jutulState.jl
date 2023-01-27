model, init_state = setup_model_state()

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
end
@testset "TLSOptions defaults" begin
    o = TLSOptions()
    @test o.R_star == 1.0
    @test o.M_star == 1.0
    @test o.u == [0.4804, 0.1867]
    @test o.oversampling_factor == 3
    @test o.duration_grid_step == 1.1
    @test o.n_transits_min == 2
end

@testset "TLSOptions overrides" begin
    o = TLSOptions(R_star = 1.2, M_star = 0.9, oversampling_factor = 5)
    @test o.R_star == 1.2
    @test o.M_star == 0.9
    @test o.oversampling_factor == 5
    @test o.u == [0.4804, 0.1867]          # unchanged
end

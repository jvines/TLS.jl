@testset "period_grid" begin
    pg = TLS.period_grid(1.0, 1.0, 30.0;
                         period_min = 0.5, period_max = 10.0,
                         oversampling_factor = 3, n_transits_min = 2)
    @test !isempty(pg)
    @test issorted(pg)
    @test first(pg) >= 0.5
    @test last(pg) <= 10.0

    # finer oversampling yields at least as many periods
    pg2 = TLS.period_grid(1.0, 1.0, 30.0;
                          period_min = 0.5, period_max = 10.0,
                          oversampling_factor = 6, n_transits_min = 2)
    @test length(pg2) >= length(pg)

    # argument validation
    @test_throws ArgumentError TLS.period_grid(-1.0, 1.0, 30.0)
    @test_throws ArgumentError TLS.period_grid(1.0, -1.0, 30.0)
    @test_throws ArgumentError TLS.period_grid(1.0, 1.0, 0.0)
end

@testset "duration_grid" begin
    periods = TLS.period_grid(1.0, 1.0, 30.0;
                              period_min = 1.0, period_max = 10.0,
                              oversampling_factor = 2)
    dg = TLS.duration_grid(periods, 1.0, 1.0; step = 1.1)
    @test !isempty(dg)
    @test issorted(dg)
    @test all(0 .< dg .< 1)
end

@testset "T14" begin
    # Central transit duration for b=0, small planet: ~13 h for Earth.
    T_earth = TLS.T14(; R_s = 1.0, M_s = 1.0, P = 365.25, upper_limit = false)
    @test 10 / 24 <= T_earth <= 16 / 24
    # upper_limit doubles it as a safety margin for duration-grid bounds.
    T_earth_upper = TLS.T14(; R_s = 1.0, M_s = 1.0, P = 365.25,
                              upper_limit = true)
    @test T_earth_upper ≈ 2 * T_earth
    # Shorter period → shorter duration (scales as P^(1/3)).
    T_hot = TLS.T14(; R_s = 1.0, M_s = 1.0, P = 1.0, upper_limit = false)
    @test T_hot < T_earth
end

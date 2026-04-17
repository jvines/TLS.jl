@testset "ldc_from_params: solar" begin
    u = TLS.ldc_from_params(5777.0, 4.44, 0.0; passband = :Kepler)
    @test u ≈ [0.4804, 0.1867] atol = 1e-4
end

@testset "ldc_from_params: interpolation monotonicity" begin
    # u1 decreases with Teff across the FGK range in the Kepler passband.
    u_cool = TLS.ldc_from_params(4000.0, 4.5, 0.0; passband = :Kepler)
    u_warm = TLS.ldc_from_params(7000.0, 4.5, 0.0; passband = :Kepler)
    @test u_cool[1] > u_warm[1]
end

@testset "ldc_from_params: out of range clamp" begin
    u_low  = TLS.ldc_from_params(2000.0, 4.5, 0.0; passband = :Kepler)
    u_high = TLS.ldc_from_params(20000.0, 4.5, 0.0; passband = :Kepler)
    # Clamped to table endpoints.
    @test u_low == TLS.ldc_from_params(3500.0, 4.5, 0.0; passband = :Kepler)
    @test u_high == TLS.ldc_from_params(8000.0, 4.5, 0.0; passband = :Kepler)
end

@testset "ldc_from_params: unknown passband" begin
    @test_throws ArgumentError TLS.ldc_from_params(5777.0, 4.5, 0.0;
                                                   passband = :JWST)
end

@testset "ldc_from_params: passband differs" begin
    @test TLS.ldc_from_params(6500.0, 4.5, 0.0; passband = :Kepler) !=
          TLS.ldc_from_params(6500.0, 4.5, 0.0; passband = :TESS)
end

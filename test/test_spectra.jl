@testset "spectra: SR definition matches Python TLS" begin
    # SR = min(chi2) / chi2, peak = 1 at argmin(chi2)
    chi2 = [10.0, 8.0, 4.0, 2.0, 6.0, 12.0]
    SR, _, _, _, _ = TLS.spectra(chi2, 3)
    @test SR[argmin(chi2)] ≈ 1.0
    @test all(0 .< SR .<= 1)
    @test SR ≈ minimum(chi2) ./ chi2
end

@testset "spectra: SDE_raw formula" begin
    # Short grid falls back to SDE = SDE_raw = (1 - mean(SR)) / std(SR)
    chi2 = [10.0, 8.0, 4.0, 2.0, 6.0, 12.0, 9.0, 7.0]
    SR, _, _, SDE_raw, SDE = TLS.spectra(chi2, 3)
    @test SDE_raw ≈ (1 - Statistics.mean(SR)) / Statistics.std(SR)
    @test SDE == SDE_raw                     # below detrend threshold
end

@testset "spectra: power peak equals SDE on long grid" begin
    # Grid long enough to trigger the detrend branch. Inject a single deep chi2
    # minimum in flat noise and check that max(power) ≈ SDE.
    import Random
    Random.seed!(99)
    n = 2500
    chi2 = 100 .+ 0.5 .* Random.randn(n)
    chi2[1234] = 50.0                          # strong single peak
    _, _, power, _, SDE = TLS.spectra(chi2, 3)
    @test isapprox(maximum(power), SDE; atol = 1e-8)
    @test SDE > 10                             # strong injection → strong SDE
end

@testset "spectra: flat chi2 yields tiny SDE" begin
    chi2 = fill(100.0, 2000) .+ 0.001 .* (1:2000)
    _, _, _, SDE_raw, SDE = TLS.spectra(chi2, 3)
    @test SDE_raw < 5
    @test SDE < 20
end

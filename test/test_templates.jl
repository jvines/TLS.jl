@testset "reference_transit normalization" begin
    ref = TLS.reference_transit()
    @test length(ref) > 1000
    @test all(0.0 .<= ref .<= 1.0)
    @test minimum(ref) == 0.0
    @test ref[length(ref) ÷ 2] < 0.01      # deepest near center
    @test ref[1] > 0.3
    @test ref[end] > 0.3
end

@testset "build_templates shape" begin
    ref = TLS.reference_transit()
    cache = TLS.build_templates([0.05, 0.1, 0.2], ref, 1024)
    @test length(cache.templates) == 3
    @test length(cache.signals) == 3
    @test length(cache.signal_sq) == 3

    for k in 1:3
        t = cache.templates[k]
        @test length(t) == 1024
        @test t[1] == 1.0
        @test t[end] == 1.0
        @test minimum(t) < 1.0

        s = cache.signals[k]
        s2 = cache.signal_sq[k]
        @test length(s) == cache.intransit_counts[k]
        @test length(s2) == cache.intransit_counts[k]
        @test all(s .>= 0)
        @test all(isapprox.(s2, s .^ 2))
    end
end

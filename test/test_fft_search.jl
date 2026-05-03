@testset "fft_search" begin
    Random.seed!(42)

    durations = [0.005, 0.01, 0.02, 0.05, 0.1, 0.2]
    ref = TLS.reference_transit()
    Nphase = 256
    templates = TLS.build_templates(durations, ref, Nphase)

    # Synthetic transit
    N = 4000
    t = sort(60 * rand(N))
    y = zeros(N)
    period = 3.5
    T0 = 1.2
    duration_frac = 0.04
    for i in 1:N
        phase = mod(t[i] - T0, period) / period
        if phase > 0.5
            phase -= 1.0
        end
        if abs(phase) < duration_frac / 2
            y[i] = -0.005 * (1 - (2 * phase / duration_frac)^2)
        end
    end
    y .+= 1e-3 .* randn(N)
    w = ones(N)
    sum_wy2 = sum(@. y * y * w)

    fy_d = Vector{Float64}(undef, Nphase); fw_d = Vector{Float64}(undef, Nphase)
    fy_h = Vector{Float64}(undef, Nphase); fw_h = Vector{Float64}(undef, Nphase)
    scratch = TLS.fft_scratch(Nphase)

    @testset "FFT off matches direct exactly" begin
        cache = TLS.build_template_fft(templates; threshold = typemax(Int))
        @test isempty(cache.fft_indices)
        for P in (2.5, 3.5, 5.0)
            pb_d = TLS.fold_and_score!(fy_d, fw_d, t, y, w, sum_wy2, P, templates)
            pb_h = TLS.fold_and_score_hybrid!(fy_h, fw_h, scratch, t, y, w,
                                              sum_wy2, P, templates, cache)
            @test isapprox(pb_d.chi2, pb_h.chi2; rtol = 1e-12, atol = 1e-14)
            @test pb_d.t0 == pb_h.t0
            @test pb_d.duration_idx == pb_h.duration_idx
            @test pb_d.depth == pb_h.depth
        end
    end

    @testset "Full FFT path matches direct to fp tol" begin
        cache = TLS.build_template_fft(templates; threshold = 0)
        @test length(cache.fft_indices) == length(durations)
        for P in (2.0, 3.5, 5.5, 7.0)
            pb_d = TLS.fold_and_score!(fy_d, fw_d, t, y, w, sum_wy2, P, templates)
            pb_h = TLS.fold_and_score_hybrid!(fy_h, fw_h, scratch, t, y, w,
                                              sum_wy2, P, templates, cache)
            @test isapprox(pb_d.chi2, pb_h.chi2; rtol = 1e-10, atol = 1e-12)
            @test pb_d.duration_idx == pb_h.duration_idx
            @test isapprox(pb_d.depth, pb_h.depth; rtol = 1e-10, atol = 1e-12)
            @test pb_d.t0 == pb_h.t0
        end
    end

    @testset "Mixed threshold matches direct" begin
        cache = TLS.build_template_fft(templates; threshold = 16)
        for P in (2.0, 3.5, 5.5)
            pb_d = TLS.fold_and_score!(fy_d, fw_d, t, y, w, sum_wy2, P, templates)
            pb_h = TLS.fold_and_score_hybrid!(fy_h, fw_h, scratch, t, y, w,
                                              sum_wy2, P, templates, cache)
            @test isapprox(pb_d.chi2, pb_h.chi2; rtol = 1e-10, atol = 1e-12)
            @test pb_d.duration_idx == pb_h.duration_idx
        end
    end

    @testset "TLSOptions.fft_threshold round-trips through tls()" begin
        # End-to-end recovery should be unchanged whether FFT is on or off.
        r_off = TLS.tls(t, y .+ 1; period_min = 1.0, period_max = 8.0,
                        oversampling_factor = 2, fft_threshold = typemax(Int))
        r_auto = TLS.tls(t, y .+ 1; period_min = 1.0, period_max = 8.0,
                         oversampling_factor = 2)
        @test isapprox(r_off.period, r_auto.period; rtol = 1e-6)
        @test isapprox(r_off.depth, r_auto.depth; rtol = 1e-3)
    end
end

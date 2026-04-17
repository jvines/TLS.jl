"""
Inject a box transit at known (P, T0, dur, depth), return (time, flux).
"""
function inject_box(; P, T0, dur, depth, N = 3000, span = 40.0, noise = 5e-4,
                    seed = 0)
    Random.seed!(seed)
    t = collect(range(0.0, span, length = N))
    f = ones(N) .+ noise .* randn(N)
    for i in 1:N
        ph = (t[i] - T0) / P
        ph -= floor(ph + 0.5)
        if abs(ph) < (dur / 2) / P
            f[i] -= depth
        end
    end
    return t, f
end

@testset "recovery: canonical case (P=3.5, 1% depth)" begin
    t, f = inject_box(; P = 3.5, T0 = 1.2, dur = 0.15, depth = 0.01)
    r = tls(t, f; period_min = 1.0, period_max = 8.0, oversampling_factor = 2)
    @test abs(r.period - 3.5) / 3.5 < 0.01
    @test r.SDE > 5
    @test r.depth > 0.005
    @test r.transit_count >= 9
end

@testset "recovery: short period (P=1.8)" begin
    t, f = inject_box(; P = 1.8, T0 = 0.4, dur = 0.08, depth = 0.015,
                        N = 4000, span = 40.0, seed = 1)
    r = tls(t, f; period_min = 1.0, period_max = 5.0, oversampling_factor = 2)
    @test abs(r.period - 1.8) / 1.8 < 0.02
    @test r.SDE > 5
end

@testset "recovery: long period (P=12)" begin
    t, f = inject_box(; P = 12.0, T0 = 4.0, dur = 0.3, depth = 0.02,
                        N = 3000, span = 80.0, seed = 2)
    r = tls(t, f; period_min = 5.0, period_max = 20.0, oversampling_factor = 2)
    @test abs(r.period - 12.0) / 12.0 < 0.02
    @test r.SDE > 5
end

@testset "no-transit baseline: low SDE" begin
    Random.seed!(42)
    N = 3000
    t = collect(range(0.0, 40.0, length = N))
    f = ones(N) .+ 1e-3 .* randn(N)
    r = tls(t, f; period_min = 1.0, period_max = 8.0, oversampling_factor = 2)
    @test r.SDE < 8                        # nothing real; spurious peak small
end

@testset "odd/even depth symmetry on symmetric transit" begin
    t, f = inject_box(; P = 2.7, T0 = 0.3, dur = 0.1, depth = 0.012,
                        N = 4000, span = 50.0, seed = 3, noise = 3e-4)
    r = tls(t, f; period_min = 1.0, period_max = 5.0, oversampling_factor = 2)
    @test abs(r.period - 2.7) / 2.7 < 0.01
    if isfinite(r.depth_mean_odd[1]) && isfinite(r.depth_mean_even[1])
        # symmetric injection → odd/even mismatch should be modest
        @test r.odd_even_mismatch < 5
    end
end

@testset "TLSResult show methods" begin
    t, f = inject_box(; P = 3.5, T0 = 1.2, dur = 0.15, depth = 0.01)
    r = tls(t, f; period_min = 1.0, period_max = 8.0, oversampling_factor = 2)
    # both MIME paths should run without error
    buf = IOBuffer()
    show(IOContext(buf, :compact => true), r)
    @test occursin("TLSResult", String(take!(buf)))
    buf = IOBuffer()
    show(buf, MIME"text/plain"(), r)
    s = String(take!(buf))
    @test occursin("period", s)
    @test occursin("SDE", s)
end

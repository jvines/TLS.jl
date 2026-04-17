@testset "hot loop produces finite outputs and is fast enough" begin
    Random.seed!(7)
    N = 3000
    t = collect(range(0.0, 40.0, length = N))
    f = ones(N) .+ 5e-4 .* randn(N)
    for i in 1:N
        ph = (t[i] - 1.2) / 3.5
        ph -= floor(ph + 0.5)
        if abs(ph) < 0.05 / 3.5
            f[i] -= 0.01
        end
    end

    # warmup
    r0 = tls(t, f; period_min = 1.0, period_max = 8.0, oversampling_factor = 2)
    @test isfinite(r0.period)

    elapsed = @elapsed tls(t, f; period_min = 1.0, period_max = 8.0,
                           oversampling_factor = 2)
    # Generous wall-clock ceiling so CI noise doesn't flake. On a
    # modern CPU single-thread we run this in ~0.5s; CI gets 15s.
    @test elapsed < 15.0
end

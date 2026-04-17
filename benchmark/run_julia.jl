using TLS
using DelimitedFiles
using Printf

data = readdlm(joinpath(@__DIR__, "lc.csv"), ',', Float64)
t, f = data[:, 1], data[:, 2]
@printf "loaded %d points; span=%.2fd; median flux=%.6f\n" length(t) (t[end]-t[1]) (sum(f)/length(f))

# warmup
TLS.tls(t, f; period_min = 1.0, period_max = 8.0, oversampling_factor = 2)

t0 = time_ns()
r = TLS.tls(t, f; period_min = 1.0, period_max = 8.0, oversampling_factor = 2)
elapsed = (time_ns() - t0) / 1e9

println("\n=== TLS.jl ===")
@printf "elapsed       : %.3f s\n" elapsed
@printf "period        : %.6f d  (true 3.5)\n" r.period
@printf "T0            : %.6f d  (true 1.2 mod P)\n" r.T0
@printf "duration      : %.6f d  (true 0.15)\n" r.duration
@printf "depth         : %.6f\n" r.depth
@printf "SDE           : %.3f\n" r.SDE
@printf "SR            : %.6f\n" r.SR
@printf "transit_count : %d\n" r.transit_count
@printf "periods_grid  : %d points\n" length(r.periods)

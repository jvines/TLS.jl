using Random, DelimitedFiles

# Reproducible synthetic light curve with an injected Mandel-Agol-style
# box transit. Written as a two-column CSV (time,flux) so Python and
# Julia see bit-identical input.
function inject(; P = 3.5, T0 = 1.2, dur = 0.15, depth = 0.01,
                N = 5000, span = 60.0, noise = 5e-4, seed = 0)
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

t, f = inject()
writedlm(joinpath(@__DIR__, "lc.csv"), hcat(t, f), ',')
println("wrote $(length(t)) rows to lc.csv")
println("true: P=3.5, T0=1.2, dur=0.15, depth=0.01")

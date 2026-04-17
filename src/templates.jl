"""
    reference_transit(; u=[0.4804, 0.1867], supersample_factor=10000,
                        rp_over_rs=0.1) -> Vector{Float64}

Compute a canonical transit flux profile with quadratic limb darkening
and small planet-to-star radius ratio. Sampled on a symmetric grid over
the in-transit region (phase ∈ [-0.5, 0.5]) and normalized so the deepest
point has depth 1 (value 0) and the out-of-transit level is 1.

This is the Mandel-Agol small-planet analytic form used by the reference
Python `transitleastsquares` package.
"""
function reference_transit(; u::AbstractVector{<:Real} = [0.4804, 0.1867],
                           supersample_factor::Integer = 10_000,
                           rp_over_rs::Real = 0.1)
    length(u) == 2 || throw(ArgumentError("quadratic LD requires 2 coefficients"))
    u1, u2 = float(u[1]), float(u[2])

    # Mandel & Agol (2002), Seager & Mallen-Ornelas (2003): for a small
    # planet the fractional flux drop at normalized impact parameter
    # z = d / R_star (center of planet to stellar center, normalized) is
    #    dF(z) = (Rp/R*)^2 * I(r) / <I>
    # where r = z for the planet-over-star-center limit, I(mu) is the LD
    # law, and <I> is the disk-averaged intensity.

    # Disk-averaged intensity for quadratic LD:
    # <I> = 1 - u1/3 - u2/6
    mean_I = 1 - u1 / 3 - u2 / 6

    # Sample phase from -0.5 to 0.5 (fraction of transit duration); the
    # planet crosses the star center at phase = 0. z = 2 * |phase| spans
    # 0 (center) to 1 (limb) linearly with phase for the canonical
    # equatorial chord.
    N = supersample_factor
    flux = Vector{Float64}(undef, N)
    rp2 = rp_over_rs^2
    @inbounds for i in 1:N
        phase = (i - 1) / (N - 1) - 0.5   # -0.5 .. 0.5
        z = min(2 * abs(phase), 1.0)      # chord radius at this phase
        μ = sqrt(max(0.0, 1 - z * z))
        I = 1 - u1 * (1 - μ) - u2 * (1 - μ)^2
        flux[i] = 1 - rp2 * I / mean_I
    end

    # Normalize to unit depth at the deepest point.
    d = 1 - minimum(flux)
    if d > 0
        @inbounds for i in eachindex(flux)
            flux[i] = 1 - (1 - flux[i]) / d
        end
    end
    return flux
end

"""
    TemplateCache

A cache of canonical transit profiles resampled at varying fractional
durations. Each row `templates[k]` is the template resampled onto a
regular phase grid of `Nphase` points covering one full period, with
the in-transit portion occupying `fractional_duration[k] * Nphase` bins
centered on phase 0.
"""
struct TemplateCache
    reference::Vector{Float64}              # canonical shape (supersample_factor pts)
    durations::Vector{Float64}              # fractional durations (per-period)
    Nphase::Int                             # phase bins per period
    templates::Vector{Vector{Float64}}      # templates[k] ∈ R^{Nphase} (full grid, OOT=1)
    signals::Vector{Vector{Float64}}        # signal s = 1 - template, compact (length nin[k])
    signal_sq::Vector{Vector{Float64}}      # s .* s for each template (reduction)
    intransit_counts::Vector{Int}           # in-transit bin counts per template
    lo_offsets::Vector{Int}                 # 0-based index of first in-transit bin in full grid
end

"""
    build_templates(durations, reference, Nphase) -> TemplateCache

Resample the canonical `reference` onto a length-`Nphase` phase grid for
each fractional duration in `durations`. Out-of-transit bins hold 1.0.
"""
function build_templates(durations::AbstractVector{<:Real},
                         reference::AbstractVector{<:Real},
                         Nphase::Integer)
    Nphase > 0 || throw(ArgumentError("Nphase must be positive"))
    ref = collect(Float64, reference)
    durs = collect(Float64, durations)
    templates = Vector{Vector{Float64}}(undef, length(durs))
    signals = Vector{Vector{Float64}}(undef, length(durs))
    signal_sq = Vector{Vector{Float64}}(undef, length(durs))
    intransit = Vector{Int}(undef, length(durs))
    lo_offsets = Vector{Int}(undef, length(durs))
    for k in eachindex(durs)
        nin = max(3, round(Int, durs[k] * Nphase))
        nin = min(nin, Nphase)
        t = ones(Float64, Nphase)
        s = Vector{Float64}(undef, nin)
        s2 = Vector{Float64}(undef, nin)
        lo = (Nphase - nin) ÷ 2                # 0-based offset
        @inbounds for j in 1:nin
            phase_ref = (j - 1) / (nin - 1)    # 0..1
            idx = clamp(round(Int, phase_ref * (length(ref) - 1)) + 1, 1, length(ref))
            val = ref[idx]
            t[lo + j] = val
            sig = 1.0 - val
            s[j] = sig
            s2[j] = sig * sig
        end
        templates[k] = t
        signals[k] = s
        signal_sq[k] = s2
        intransit[k] = nin
        lo_offsets[k] = lo
    end
    return TemplateCache(ref, durs, Nphase, templates, signals, signal_sq,
                         intransit, lo_offsets)
end

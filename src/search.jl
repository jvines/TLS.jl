"""
    PeriodBest

Per-period best-fit summary produced by the hot loop.
"""
struct PeriodBest
    chi2::Float64
    t0::Float64
    duration_idx::Int
    depth::Float64
end

"""
    fold_and_score!(scratch, time, y, w, period, templates) -> PeriodBest

Core per-period operation. For the given period, phase-fold `(time, y, w)`
into `scratch.folded_y` and `scratch.folded_w`, then for each template
slide it across all phase offsets to find the minimum χ². Returns the
best (χ², T0, duration_idx, depth).

`w` is inverse-variance weight (1/σ²). `y` should already be
mean-centered (flux - 1).
"""
function fold_and_score!(folded_y::Vector{Float64},
                         folded_w::Vector{Float64},
                         time::AbstractVector{<:Real},
                         y::AbstractVector{<:Real},
                         w::AbstractVector{<:Real},
                         sum_wy2::Real,
                         period::Real,
                         templates::TemplateCache)
    Nphase = templates.Nphase
    fill!(folded_y, 0.0)
    fill!(folded_w, 0.0)

    invP = 1.0 / period
    @inbounds for i in eachindex(time)
        φ = time[i] * invP
        φ -= floor(φ)
        bin = unsafe_trunc(Int, φ * Nphase) + 1
        if bin > Nphase
            bin = Nphase
        end
        folded_w[bin] += w[i]
        folded_y[bin] += y[i] * w[i]
    end

    best_chi2 = Inf
    best_offset = 0
    best_k = 1
    best_depth = 0.0

    @inbounds for k in eachindex(templates.signals)
        s  = templates.signals[k]
        s2 = templates.signal_sq[k]
        nin = templates.intransit_counts[k]
        # For each offset o ∈ [0, Nphase), the transit spans bins
        #   ((lo + j - 1 + o) mod Nphase) + 1  for j = 1..nin
        # We only iterate over nin bins per offset — the hot-loop cost is
        # O(Nphase * nin * Ndur) per period rather than O(Nphase^2 * Ndur).
        for o in 0:(Nphase - 1)
            # unrolled to avoid `%` cost when wrap is not needed
            start = o + 1       # first 1-based bin index of the transit window
            # two contiguous segments: [start, min(start+nin-1, Nphase)] and
            # [1, nin - (Nphase - start + 1)] (the wrap-around).
            last_no_wrap = min(start + nin - 1, Nphase)
            len1 = last_no_wrap - start + 1
            len2 = nin - len1

            sum_wys = 0.0
            sum_ws2 = 0.0
            @simd for j in 1:len1
                bin = start + j - 1
                sum_wys += folded_y[bin] * s[j]
                sum_ws2 += folded_w[bin] * s2[j]
            end
            if len2 > 0
                @simd for j in 1:len2
                    bin = j
                    sj = s[len1 + j]
                    sum_wys += folded_y[bin] * sj
                    sum_ws2 += folded_w[bin] * s2[len1 + j]
                end
            end
            sum_ws2 > 0 || continue
            depth = sum_wys / sum_ws2
            chi2 = sum_wy2 - 2 * depth * sum_wys + depth * depth * sum_ws2
            if chi2 < best_chi2
                best_chi2 = chi2
                best_offset = o
                best_k = k
                best_depth = depth
            end
        end
    end

    # Convert offset to T0. With offset `o`, the compact signal s[j] lines
    # up with data bin `o + j` (1-based), so the transit midpoint in data
    # coordinates sits at bin (o + (nin+1)/2). Convert to phase fraction.
    nin = templates.intransit_counts[best_k]
    mid_bin = mod(best_offset + (nin + 1) / 2 - 0.5, Nphase)
    best_t0 = (mid_bin / Nphase) * period

    return PeriodBest(best_chi2, best_t0, best_k, best_depth)
end

"""
    _tls(time, flux, flux_err, opts) -> TLSResult

Top-level pipeline. Computes period grid, duration grid, template cache,
runs the threaded period search, extracts statistics, and assembles the
result struct.
"""
function _tls(time::Vector{Float64},
              flux::Vector{Float64},
              flux_err::Union{Nothing,Vector{Float64}},
              opts::TLSOptions)
    length(time) == length(flux) ||
        throw(ArgumentError("time and flux must be same length"))
    length(time) >= 4 ||
        throw(ArgumentError("need at least 4 samples"))

    # normalize flux
    normalize_flux!(flux)

    # weights
    w = if flux_err === nothing
        # initial uniform weights; updated below using out-of-transit std
        ones(Float64, length(flux))
    else
        length(flux_err) == length(flux) ||
            throw(ArgumentError("flux_err length mismatch"))
        [fe > 0 ? inv(fe * fe) : 0.0 for fe in flux_err]
    end

    y = similar(flux)
    @inbounds for i in eachindex(flux)
        y[i] = flux[i] - 1.0
    end

    baseline = maximum(time) - minimum(time)
    period_min = opts.period_min === nothing ? 0.0 : opts.period_min
    period_max = opts.period_max === nothing ? baseline / max(1, opts.n_transits_min) :
                 opts.period_max
    periods = period_grid(opts.R_star, opts.M_star, baseline;
                          period_min = period_min,
                          period_max = period_max,
                          oversampling_factor = opts.oversampling_factor,
                          n_transits_min = opts.n_transits_min)

    isempty(periods) && throw(ArgumentError(
        "period grid is empty; check period_min/period_max and baseline"))

    durations = duration_grid(periods, opts.R_star, opts.M_star;
                              step = opts.duration_grid_step)

    ref = reference_transit(u = opts.u)
    # heuristic: bin count per period ~ samples-per-period at shortest period,
    # bounded to a reasonable range.
    Nphase = clamp(round(Int, length(time) / opts.n_transits_min / 1), 256, 4096)
    templates = build_templates(durations, ref, Nphase)

    nperiods = length(periods)
    chi2 = Vector{Float64}(undef, nperiods)
    t0_best = Vector{Float64}(undef, nperiods)
    k_best = Vector{Int}(undef, nperiods)
    depth_best = Vector{Float64}(undef, nperiods)

    # per-thread scratch (allocate for the max addressable thread id)
    nt = Threads.maxthreadid()
    folded_ys = [Vector{Float64}(undef, Nphase) for _ in 1:nt]
    folded_ws = [Vector{Float64}(undef, Nphase) for _ in 1:nt]

    # sum_wy^2 is period-independent; precompute once.
    sum_wy2 = 0.0
    @inbounds for i in eachindex(y)
        sum_wy2 += y[i] * y[i] * w[i]
    end

    Threads.@threads :static for ip in 1:nperiods
        tid = Threads.threadid()
        fy = folded_ys[tid]
        fw = folded_ws[tid]
        pb = fold_and_score!(fy, fw, time, y, w, sum_wy2, periods[ip], templates)
        chi2[ip] = pb.chi2
        t0_best[ip] = pb.t0
        k_best[ip] = pb.duration_idx
        depth_best[ip] = pb.depth
    end

    # Build the power spectrum and pick best period using the Python-TLS
    # normalization (see `spectra`).
    idx_best = argmin(chi2)
    SR, power_raw, power, SDE_raw, SDE = spectra(chi2, opts.oversampling_factor)

    best_period = periods[idx_best]
    best_t0 = t0_best[idx_best]
    best_k = k_best[idx_best]
    best_duration = durations[best_k] * best_period
    best_depth = depth_best[idx_best]

    res = TLSResult(
        periods = periods,
        power = power,
        power_raw = power_raw,
        SDE = SDE,
        SDE_raw = SDE_raw,
        SR = maximum(SR),
        chi2 = chi2,
        chi2red = chi2 ./ max(1, length(time) - 1),
        period = best_period,
        T0 = best_t0,
        duration = best_duration,
        depth = -best_depth,
        rs_of_periods = [durations[k] for k in k_best],
    )

    postprocess!(res, time, flux, w, templates, opts)
    return res
end

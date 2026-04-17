"""
    postprocess!(res, time, flux, w, templates, opts)

Fill in the remaining fields of `res` that were not computed in the
hot loop: per-transit depths, odd/even depths, transit count, SNR, and
the model light-curve arrays.

This is a v0.1 implementation that covers the common fields; deeper
parity with Python TLS (e.g. `snr_pink_per_transit`) will land in v0.2.
"""
function postprocess!(res::TLSResult,
                      time::AbstractVector{<:Real},
                      flux::AbstractVector{<:Real},
                      w::AbstractVector{<:Real},
                      templates::TemplateCache,
                      opts::TLSOptions)
    P = res.period
    T0 = res.T0
    dur = res.duration

    # in-transit mask and per-transit window
    phase = similar(time, Float64)
    @inbounds for i in eachindex(time)
        ph = (time[i] - T0) / P
        ph -= floor(ph + 0.5)       # (-0.5, 0.5]
        phase[i] = ph
    end
    half_dur_phase = (dur / P) / 2
    in_transit = [abs(p) <= half_dur_phase for p in phase]
    res.in_transit_count = count(in_transit)

    # count distinct transit epochs
    t_min = minimum(time)
    t_max = maximum(time)
    epoch_lo = floor((t_min - T0) / P)
    epoch_hi = ceil((t_max - T0) / P)
    transit_epochs = Int[]
    for e in epoch_lo:epoch_hi
        tc = T0 + e * P
        if tc >= t_min - dur && tc <= t_max + dur
            push!(transit_epochs, Int(e))
        end
    end
    res.transit_count = length(transit_epochs)

    depths = Float64[]
    depth_errs = Float64[]
    per_counts = Int[]
    distinct = 0
    empty = 0
    for e in transit_epochs
        tc = T0 + e * P
        sum_w = 0.0
        sum_wy = 0.0
        cnt = 0
        @inbounds for i in eachindex(time)
            if abs(time[i] - tc) <= dur / 2
                sum_w += w[i]
                sum_wy += w[i] * (1 - flux[i])
                cnt += 1
            end
        end
        push!(per_counts, cnt)
        if cnt > 0 && sum_w > 0
            d = sum_wy / sum_w
            push!(depths, d)
            push!(depth_errs, sqrt(1 / sum_w))
            distinct += 1
        else
            push!(depths, NaN)
            push!(depth_errs, NaN)
            empty += 1
        end
    end
    res.transit_depths = depths
    res.transit_depths_uncertainties = depth_errs
    res.per_transit_count = per_counts
    res.distinct_transit_count = distinct
    res.empty_transit_count = empty

    finite_depths = filter(isfinite, depths)
    if !isempty(finite_depths)
        μ = mean(finite_depths)
        σ = length(finite_depths) > 1 ? std(finite_depths) / sqrt(length(finite_depths)) : NaN
        res.depth_mean = (μ, σ)

        odd  = [d for (e, d) in zip(transit_epochs, depths) if isfinite(d) && isodd(e)]
        even = [d for (e, d) in zip(transit_epochs, depths) if isfinite(d) && iseven(e)]
        if !isempty(odd)
            res.depth_mean_odd = (mean(odd),
                                  length(odd) > 1 ? std(odd) / sqrt(length(odd)) : NaN)
        end
        if !isempty(even)
            res.depth_mean_even = (mean(even),
                                   length(even) > 1 ? std(even) / sqrt(length(even)) : NaN)
        end
        if !isempty(odd) && !isempty(even)
            δ = res.depth_mean_odd[1] - res.depth_mean_even[1]
            σ_combo = sqrt(res.depth_mean_odd[2]^2 + res.depth_mean_even[2]^2)
            res.odd_even_mismatch = σ_combo > 0 ? abs(δ) / σ_combo : 0.0
        end
    end

    # folded series
    order = sortperm(phase)
    res.folded_phase = phase[order]
    res.folded_y = flux[order]
    # dy from weights (σ = sqrt(1/w))
    res.folded_dy = [w_i > 0 ? sqrt(1 / w_i) : NaN for w_i in w[order]]

    # SNR from overall depth and out-of-transit scatter
    oot = .!in_transit
    σ_oot = std(flux[oot])
    n_in = res.in_transit_count
    if σ_oot > 0 && n_in > 0
        res.snr = res.depth * sqrt(n_in) / σ_oot
    end

    return res
end

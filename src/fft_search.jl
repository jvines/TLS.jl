using FFTW

"""
    TemplateFFTCache

FFT-domain cache for templates whose in-transit length `nin` is large
enough that an FFT-based circular cross-correlation beats the direct
O(Nphase * nin) inner loop.

For each cached template `k` we store
- `F_signals[k]    = conj(rfft(pad(s_k,  Nphase)))`
- `F_signal_sq[k]  = conj(rfft(pad(s2_k, Nphase)))`

`fft_indices[i]` gives the index into the parent `TemplateCache` for
the i-th cached template; templates with `nin < threshold` are not
cached and continue to use the direct path.

We deliberately store kernels as separate `Vector{ComplexF64}` rather
than columns of a `Matrix`: in microbenchmarks the matrix-column form
inhibits @simd vectorization of the per-template product loop and
costs ~25% throughput. We also avoid FFTW's batched IFFT plan; for
non-power-of-two `Nphase` (e.g. 2500) the batched codepath is up to
7× slower than a loop of independent 1D IFFTs.
"""
struct TemplateFFTCache
    Nphase::Int
    threshold::Int
    fft_indices::Vector{Int}
    F_signals::Vector{Vector{ComplexF64}}
    F_signal_sq::Vector{Vector{ComplexF64}}
    rfft_plan::FFTW.rFFTWPlan{Float64,-1,false,1}
    irfft_plan::AbstractFFTs.ScaledPlan
end

"""
    build_template_fft(templates; threshold) -> TemplateFFTCache

Precompute FFT-domain copies of every template in `templates` whose
`intransit_counts[k] >= threshold`. Per period the two forward FFTs
of the folded data are amortized across all FFT-cached templates, so
the marginal FFT cost per template is one IFFT plus a length-Nphase
scan: ~`Nphase * log2(Nphase)` ops, vs `Nphase * nin` for the direct
SIMD loop. The asymptotic crossover sits near `nin ≈ log2(Nphase)`,
but FFTW's per-IFFT constant is larger than the @simd direct-loop
constant, so the empirical optimum on a power-of-2 `Nphase` is closer
to `1.5 · log2(Nphase)`.

Pass `threshold = 0` to force every template through the FFT path,
or `typemax(Int)` to disable FFT entirely.
"""
function build_template_fft(templates::TemplateCache;
                            threshold::Integer = ceil(Int, 1.5 * log2(max(2, templates.Nphase))))
    Nphase = templates.Nphase
    nf = Nphase ÷ 2 + 1
    rplan = plan_rfft(Vector{Float64}(undef, Nphase); flags = FFTW.MEASURE)
    iplan = plan_irfft(Vector{ComplexF64}(undef, nf), Nphase; flags = FFTW.MEASURE)

    indices = Int[]
    F_sigs = Vector{Vector{ComplexF64}}()
    F_sig_sq = Vector{Vector{ComplexF64}}()

    pad_buf = Vector{Float64}(undef, Nphase)
    for k in eachindex(templates.signals)
        nin = templates.intransit_counts[k]
        nin >= threshold || continue
        s  = templates.signals[k]
        s2 = templates.signal_sq[k]

        fill!(pad_buf, 0.0)
        @inbounds for j in 1:nin
            pad_buf[j] = s[j]
        end
        Fs = rplan * pad_buf
        @inbounds for i in eachindex(Fs)
            Fs[i] = conj(Fs[i])
        end

        fill!(pad_buf, 0.0)
        @inbounds for j in 1:nin
            pad_buf[j] = s2[j]
        end
        Fs2 = rplan * pad_buf
        @inbounds for i in eachindex(Fs2)
            Fs2[i] = conj(Fs2[i])
        end

        push!(indices, k)
        push!(F_sigs, Fs)
        push!(F_sig_sq, Fs2)
    end

    return TemplateFFTCache(Nphase, Int(threshold), indices, F_sigs, F_sig_sq,
                            rplan, iplan)
end

"""
    fft_scratch(Nphase) -> NamedTuple

Per-thread scratch buffers for the FFT inner-loop path. The single
`F_tmp` complex buffer is reused across templates inside one call to
`fold_and_score_hybrid!`; `wys` / `ws2` hold the IFFT outputs for the
current template before the chi² scan consumes them.
"""
function fft_scratch(Nphase::Integer)
    nf = Nphase ÷ 2 + 1
    (F_y    = Vector{ComplexF64}(undef, nf),
     F_w    = Vector{ComplexF64}(undef, nf),
     F_tmp  = Vector{ComplexF64}(undef, nf),
     wys    = Vector{Float64}(undef, Nphase),
     ws2    = Vector{Float64}(undef, Nphase))
end

"""
    fold_and_score_hybrid!(folded_y, folded_w, scratch, time, y, w,
                           sum_wy2, period, templates, fft_cache) -> PeriodBest

Phase-fold the data and search every template-offset pair for the
minimum χ². Templates with `nin >= fft_cache.threshold` are evaluated
through a single batched FFT cross-correlation; the rest use the
direct SIMD path identical to `fold_and_score!`.
"""
function fold_and_score_hybrid!(folded_y::Vector{Float64},
                                folded_w::Vector{Float64},
                                scratch,
                                time::AbstractVector{<:Real},
                                y::AbstractVector{<:Real},
                                w::AbstractVector{<:Real},
                                sum_wy2::Real,
                                period::Real,
                                templates::TemplateCache,
                                fft_cache::TemplateFFTCache)
    Nphase = templates.Nphase
    @assert fft_cache.Nphase == Nphase
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

    # Track the best in terms of the partial form `-sum_wys² / sum_ws2`,
    # which is monotone with χ² (full χ² = sum_wy2 + this value). Skipping
    # the constant `sum_wy2` and the redundant `2·depth·sum_wys` term tightens
    # the inner argmin.
    best_partial = Inf
    best_offset = 0
    best_k = 1
    best_depth = 0.0

    K = length(fft_cache.fft_indices)
    threshold = fft_cache.threshold

    if K > 0
        F_y, F_w = scratch.F_y, scratch.F_w
        F_tmp = scratch.F_tmp
        wys_buf, ws2_buf = scratch.wys, scratch.ws2
        nf = length(F_y)

        mul!(F_y, fft_cache.rfft_plan, folded_y)
        mul!(F_w, fft_cache.rfft_plan, folded_w)

        @inbounds for ki in 1:K
            k = fft_cache.fft_indices[ki]
            Fs  = fft_cache.F_signals[ki]
            Fs2 = fft_cache.F_signal_sq[ki]

            @simd for i in 1:nf
                F_tmp[i] = F_y[i] * Fs[i]
            end
            mul!(wys_buf, fft_cache.irfft_plan, F_tmp)

            @simd for i in 1:nf
                F_tmp[i] = F_w[i] * Fs2[i]
            end
            mul!(ws2_buf, fft_cache.irfft_plan, F_tmp)

            @simd for o in 1:Nphase
                sum_wys = wys_buf[o]
                sum_ws2 = ws2_buf[o]
                if sum_ws2 > 0
                    depth = sum_wys / sum_ws2
                    partial = -sum_wys * depth
                    if partial < best_partial
                        best_partial = partial
                        best_offset = o - 1
                        best_k = k
                        best_depth = depth
                    end
                end
            end
        end
    end

    @inbounds for k in eachindex(templates.signals)
        nin = templates.intransit_counts[k]
        nin >= threshold && continue
        s  = templates.signals[k]
        s2 = templates.signal_sq[k]
        for o in 0:(Nphase - 1)
            start = o + 1
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
                    sum_wys += folded_y[bin] * s[len1 + j]
                    sum_ws2 += folded_w[bin] * s2[len1 + j]
                end
            end
            sum_ws2 > 0 || continue
            depth = sum_wys / sum_ws2
            partial = -sum_wys * depth
            if partial < best_partial
                best_partial = partial
                best_offset = o
                best_k = k
                best_depth = depth
            end
        end
    end

    nin = templates.intransit_counts[best_k]
    mid_bin = mod(best_offset + (nin + 1) / 2 - 0.5, Nphase)
    best_t0 = (mid_bin / Nphase) * period
    best_chi2 = sum_wy2 + best_partial

    return PeriodBest(best_chi2, best_t0, best_k, best_depth)
end

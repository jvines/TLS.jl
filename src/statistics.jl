"""
    spectra(chi2, oversampling_factor) -> (SR, power_raw, power, SDE_raw, SDE)

Reference implementation matching
`transitleastsquares.stats.spectra` (Python TLS v1.32, Apr 2024):

    SR       = min(chi2) / chi2                           (∈ (0, 1], peak = 1)
    SDE_raw  = (1 - mean(SR)) / std(SR)
    power_raw = (SR - mean(SR)) * SDE_raw / max(SR - mean(SR))
    kernel   = oversampling_factor * 30 (forced odd)
    power    = power_raw - running_median(power_raw, kernel)
    power    = power - mean(power)
    SDE      = max(power) / std(power)
    power    = power * (SDE / max(power))                 (so max(power) == SDE)

Falls back to `power = power_raw`, `SDE = SDE_raw` when the grid is too
short for the detrend pass (< 2*kernel points).
"""
function spectra(chi2::AbstractVector{<:Real}, oversampling_factor::Integer)
    chi2_min = minimum(chi2)
    n = length(chi2)

    SR = Vector{Float64}(undef, n)
    @inbounds for i in 1:n
        SR[i] = chi2[i] > 0 ? chi2_min / chi2[i] : 0.0
    end

    meanSR = mean(SR)
    stdSR = std(SR)
    SDE_raw = stdSR > 0 ? (1 - meanSR) / stdSR : 0.0

    power_raw = Vector{Float64}(undef, n)
    @inbounds for i in 1:n
        power_raw[i] = SR[i] - meanSR
    end
    maxpr = maximum(power_raw)
    if maxpr > 0
        s = SDE_raw / maxpr
        @inbounds for i in 1:n
            power_raw[i] *= s
        end
    end

    kernel = oversampling_factor * SDE_MEDIAN_KERNEL_SIZE
    if iseven(kernel)
        kernel += 1
    end

    if n > 2 * kernel
        trend = running_median(power_raw, kernel)
        power = Vector{Float64}(undef, n)
        @inbounds for i in 1:n
            power[i] = power_raw[i] - trend[i]
        end
        meanp = mean(power)
        @inbounds for i in 1:n
            power[i] -= meanp
        end
        stdp = std(power)
        maxp = maximum(power)
        if stdp > 0 && maxp > 0
            SDE = maxp / stdp
            s = SDE / maxp
            @inbounds for i in 1:n
                power[i] *= s
            end
        else
            SDE = 0.0
        end
        return SR, power_raw, power, SDE_raw, SDE
    else
        return SR, power_raw, copy(power_raw), SDE_raw, SDE_raw
    end
end

const SDE_MEDIAN_KERNEL_SIZE = 30

"""
    running_median(x, window) -> Vector{Float64}

Running median with a centered window of `window` samples (forced odd,
clipped at array boundaries).
"""
function running_median(x::AbstractVector{<:Real}, window::Integer)
    n = length(x)
    window = max(3, window | 1)
    half = window ÷ 2
    out = Vector{Float64}(undef, n)
    buf = Vector{Float64}(undef, window)
    @inbounds for i in 1:n
        lo = max(1, i - half)
        hi = min(n, i + half)
        k = hi - lo + 1
        for j in 1:k
            buf[j] = x[lo + j - 1]
        end
        out[i] = median(view(buf, 1:k))
    end
    return out
end

# Kept for backward compatibility with earlier callers; the new pipeline
# calls `spectra(chi2, oversampling_factor)` instead.
function chi2_to_power(chi2::AbstractVector{<:Real})
    _, power_raw, _, _, _ = spectra(chi2, 3)
    return power_raw
end

function detrend_power(power_raw::AbstractVector{<:Real}; window_frac::Real = 0.05)
    n = length(power_raw)
    w = max(3, round(Int, window_frac * n))
    trend = running_median(power_raw, w)
    resid = Vector{Float64}(undef, n)
    @inbounds for i in 1:n
        resid[i] = power_raw[i] - trend[i]
    end
    σ = std(resid)
    σ = σ > 0 ? σ : 1.0
    @inbounds for i in 1:n
        resid[i] /= σ
    end
    return resid
end

"""
    normalize_flux!(flux) -> median

Divide `flux` by its median in-place and return the median used.
Safe for flux already normalized to unity (median ~= 1).
"""
function normalize_flux!(flux::AbstractVector{<:AbstractFloat})
    m = median(flux)
    m > 0 || throw(ArgumentError("flux median must be positive"))
    @inbounds for i in eachindex(flux)
        flux[i] /= m
    end
    return m
end

"""
    out_of_transit_std(flux, in_transit_mask) -> σ

Standard deviation of the out-of-transit samples. Used as a noise
estimator when no `flux_err` is supplied.
"""
function out_of_transit_std(flux::AbstractVector{<:Real},
                            in_transit::AbstractVector{Bool})
    n = count(!, in_transit)
    n > 1 || return NaN
    μ = 0.0
    @inbounds for i in eachindex(flux)
        in_transit[i] || (μ += flux[i])
    end
    μ /= n
    s = 0.0
    @inbounds for i in eachindex(flux)
        if !in_transit[i]
            d = flux[i] - μ
            s += d * d
        end
    end
    return sqrt(s / (n - 1))
end

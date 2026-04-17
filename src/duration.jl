"""
    T14(; R_s=1.0, M_s=1.0, P=1.0, small=false) -> duration (days)

Maximum transit duration for a planet grazing the stellar limb (b=1 for
total; b=0 for full) assuming a small planet.
"""
function T14(; R_s::Real = 1.0, M_s::Real = 1.0, P::Real = 1.0,
             upper_limit::Bool = true)
    G     = 6.67384e-11
    R_sun = 696_342e3
    M_sun = 1.98892e30
    SECDAY = 86_400.0
    P_s = P * SECDAY
    M = M_s * M_sun
    R = R_s * R_sun
    # central-transit duration (b=0, small planet, circular orbit):
    #   T = R_star * P / (π * a),    a = (G M P^2 / 4π²)^(1/3)
    # = R_star / π * (2π)^(2/3) * (P / (G M))^(1/3).
    # upper_limit=true doubles this as a safe grazing-configuration bound
    # for the duration-grid construction.
    T = (R / π) * (2π)^(2/3) * (P_s / (G * M))^(1/3)
    return upper_limit ? 2 * T / SECDAY : T / SECDAY
end

"""
    duration_grid(periods, R_s, M_s; step=1.1, log_step=true) -> Vector{Float64}

Fractional duration grid (duration / period). Geometric spacing with the
given step factor. Clamped so that physical T14 at the shortest period
is the upper bound and the shortest resolvable duration at the longest
period is the lower bound.
"""
function duration_grid(periods::AbstractVector{<:Real},
                       R_s::Real = 1.0, M_s::Real = 1.0;
                       step::Real = 1.1,
                       N_min::Integer = 2)
    isempty(periods) && return Float64[]
    p_min, p_max = extrema(periods)
    # longest physical transit at shortest period, expressed as fraction
    T_upper = T14(; R_s = R_s, M_s = M_s, P = p_min, upper_limit = true) / p_min
    T_upper = clamp(T_upper, 1e-4, 0.5)
    # shortest resolvable duration: floor at ~2 cadence points equivalent.
    # Without knowing cadence here we pick a conservative 1e-4 * p_max / p_min
    T_lower = 1e-4
    T_lower < T_upper || return [T_upper]

    n = max(N_min, ceil(Int, log(T_upper / T_lower) / log(step)) + 1)
    out = Vector{Float64}(undef, n)
    f = (T_upper / T_lower)^(1 / (n - 1))
    @inbounds for i in 1:n
        out[i] = T_lower * f^(i - 1)
    end
    return out
end

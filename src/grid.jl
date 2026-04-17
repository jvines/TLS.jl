"""
    period_grid(R_star, M_star, time_span; period_min=0.0,
                period_max=Inf, oversampling_factor=3, n_transits_min=2)

Optimal period grid of Ofir (2014), uniform in `f^(1/3)`. All units:
solar radii, solar masses, days. Returns ascending periods in days.
"""
function period_grid(R_star::Real, M_star::Real, time_span::Real;
                     period_min::Real = 0.0,
                     period_max::Real = Inf,
                     oversampling_factor::Integer = 3,
                     n_transits_min::Integer = 2)
    R_star > 0 || throw(ArgumentError("R_star must be positive"))
    M_star > 0 || throw(ArgumentError("M_star must be positive"))
    time_span > 0 || throw(ArgumentError("time_span must be positive"))
    n_transits_min >= 1 || throw(ArgumentError("n_transits_min must be >= 1"))

    R_sun = 696_342e3          # m
    M_sun = 1.98892e30         # kg
    G     = 6.67384e-11        # m^3 kg^-1 s^-2
    SECDAY = 86_400.0

    R_star_m = R_star * R_sun
    M_star_kg = M_star * M_sun
    span_s = time_span * SECDAY

    f_min = n_transits_min / span_s
    f_max = (1 / (2π)) * sqrt(G * M_star_kg / (3 * R_star_m)^3)

    A = (2π)^(2/3) / π * R_star_m / (G * M_star_kg)^(1/3) /
        (span_s * oversampling_factor)
    C = f_min^(1/3) - A / 3
    N_opt_f = (f_max^(1/3) - f_min^(1/3) + A / 3) * 3 / A
    N_opt = max(1, floor(Int, N_opt_f))

    freqs_cbrt = Vector{Float64}(undef, N_opt)
    @inbounds for i in 1:N_opt
        freqs_cbrt[i] = A / 3 * i + C
    end

    periods_s = Vector{Float64}(undef, N_opt)
    @inbounds for i in 1:N_opt
        periods_s[i] = 1 / freqs_cbrt[i]^3
    end

    periods_days = periods_s ./ SECDAY
    sort!(periods_days)

    lo = isfinite(period_min) ? Float64(period_min) : 0.0
    hi = isfinite(period_max) ? Float64(period_max) : Inf
    lo_idx = searchsortedfirst(periods_days, lo)
    hi_idx = searchsortedlast(periods_days, hi)
    return periods_days[lo_idx:hi_idx]
end

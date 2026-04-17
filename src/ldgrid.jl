"""
    ldc_solar() -> Vector{Float64}

Default quadratic limb-darkening coefficients (Claret 2011 solar,
Kepler passband).
"""
ldc_solar() = [0.4804, 0.1867]

# Built-in quadratic limb-darkening grid. Coefficients from Claret (2011,
# A&A 529, A75) tabulated at logg = 4.5, [Fe/H] = 0.0 for three passbands.
# Linear interpolation in Teff within [3500, 8000] K; outside the range
# we clamp to the nearest endpoint. A full Artifact-backed grid spanning
# (Teff, logg, [Fe/H]) is on the v0.2 roadmap.
const _LDC_TEFF = [3500.0, 4000.0, 4500.0, 5000.0, 5500.0, 5777.0,
                   6000.0, 6500.0, 7000.0, 7500.0, 8000.0]

const _LDC_KEPLER = [
    # u1,     u2
    (0.750, 0.150),
    (0.680, 0.140),
    (0.600, 0.170),
    (0.530, 0.180),
    (0.480, 0.190),
    (0.4804, 0.1867),
    (0.430, 0.220),
    (0.360, 0.280),
    (0.290, 0.320),
    (0.250, 0.330),
    (0.210, 0.330),
]

const _LDC_TESS = [
    (0.620, 0.140),
    (0.560, 0.150),
    (0.490, 0.170),
    (0.430, 0.190),
    (0.385, 0.205),
    (0.3848, 0.2048),
    (0.350, 0.220),
    (0.290, 0.270),
    (0.230, 0.310),
    (0.195, 0.325),
    (0.160, 0.330),
]

const _LDC_K2 = _LDC_KEPLER   # K2 uses the same Kepler-like broad bandpass

function _lookup_ldc(table::Vector{Tuple{Float64,Float64}}, Teff::Real)
    n = length(_LDC_TEFF)
    if Teff <= _LDC_TEFF[1]
        return [table[1][1], table[1][2]]
    elseif Teff >= _LDC_TEFF[end]
        return [table[end][1], table[end][2]]
    end
    i = searchsortedlast(_LDC_TEFF, Teff)
    t0 = _LDC_TEFF[i]
    t1 = _LDC_TEFF[i + 1]
    α  = (Teff - t0) / (t1 - t0)
    u1 = (1 - α) * table[i][1] + α * table[i + 1][1]
    u2 = (1 - α) * table[i][2] + α * table[i + 1][2]
    return [u1, u2]
end

"""
    ldc_from_params(Teff, logg, FeH; passband=:Kepler) -> Vector{Float64}

Interpolate quadratic limb-darkening coefficients from a built-in
Claret 2011 grid for the given stellar effective temperature and
passband. `logg` and `FeH` are accepted for API compatibility but the
v0.1 grid is indexed by Teff only (at solar logg and metallicity); the
full (Teff, logg, [Fe/H]) Artifact grid lands in v0.2.

Supported passbands: `:Kepler`, `:TESS`, `:K2`.
"""
function ldc_from_params(Teff::Real, logg::Real, FeH::Real;
                         passband::Symbol = :Kepler)
    table = if passband === :Kepler
        _LDC_KEPLER
    elseif passband === :TESS
        _LDC_TESS
    elseif passband === :K2
        _LDC_K2
    else
        throw(ArgumentError("passband must be :Kepler, :TESS, or :K2 (got $passband)"))
    end
    return _lookup_ldc(table, Teff)
end

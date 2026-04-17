using HTTP
using JSON3
using Scratch: @get_scratch!

struct CatalogLookupError <: Exception
    target::String
    status::Int
    message::String
end

Base.showerror(io::IO, e::CatalogLookupError) =
    print(io, "CatalogLookupError(", e.target, ", HTTP ", e.status, "): ", e.message)

const MAST_TIC_URL = "https://mast.stsci.edu/api/v0.1/catalogs/tic/"
const VIZIER_TAP_URL = "https://tapvizier.cds.unistra.fr/TAPVizieR/tap/sync"

"""
    catalog_info(; TIC=nothing, KIC=nothing, EPIC=nothing) -> NamedTuple

Look up stellar parameters for a target by TIC (MAST), KIC, or EPIC
(Vizier). Returns a NamedTuple containing `R_star`, `R_star_min`,
`R_star_max`, `M_star`, `M_star_min`, `M_star_max`, `Teff`, `logg`,
`FeH`, and limb-darkening coefficients `u` suitable for the relevant
passband.

Results are cached in a package-local scratch directory keyed by target
so repeated lookups are free.

If the upstream service returns no row or a row with missing fields,
solar defaults are substituted for the missing values.
"""
function catalog_info(; TIC::Union{Nothing,Integer} = nothing,
                      KIC::Union{Nothing,Integer} = nothing,
                      EPIC::Union{Nothing,Integer} = nothing)
    n_spec = count(!isnothing, (TIC, KIC, EPIC))
    n_spec == 1 || throw(ArgumentError(
        "catalog_info: specify exactly one of TIC, KIC, or EPIC"))

    if TIC !== nothing
        return _catalog_tic(Int(TIC))
    elseif KIC !== nothing
        return _catalog_kic(Int(KIC))
    else
        return _catalog_epic(Int(EPIC))
    end
end

function _default_stellar(passband::Symbol)
    return (R_star = 1.0, R_star_min = 0.13, R_star_max = 3.5,
            M_star = 1.0, M_star_min = 0.1,  M_star_max = 1.0,
            Teff = 5778.0, logg = 4.44, FeH = 0.0,
            u = ldc_from_params(5778.0, 4.44, 0.0; passband = passband))
end

function _catalog_tic(tic::Int)
    cache_file = joinpath(@get_scratch!("tic_cache"), string(tic, ".json"))
    payload = if isfile(cache_file)
        read(cache_file, String)
    else
        url = string(MAST_TIC_URL, tic)
        resp = HTTP.get(url; status_exception = false,
                        retry = false, readtimeout = 30)
        if resp.status != 200
            throw(CatalogLookupError("TIC $tic", resp.status,
                                     String(resp.body)))
        end
        body = String(resp.body)
        write(cache_file, body)
        body
    end

    data = JSON3.read(payload)
    row = _first_tic_row(data, tic)
    R = _coalesce_float(get(row, :rad, nothing), 1.0)
    R_lo = _coalesce_float(get(row, :e_rad, nothing), 0.0)
    M = _coalesce_float(get(row, :mass, nothing), 1.0)
    M_lo = _coalesce_float(get(row, :e_mass, nothing), 0.0)
    Teff = _coalesce_float(get(row, :Teff, nothing), 5778.0)
    logg = _coalesce_float(get(row, :logg, nothing), 4.44)
    FeH  = _coalesce_float(get(row, :MH, nothing), 0.0)
    return (R_star = R, R_star_min = max(R - R_lo, 0.0), R_star_max = R + R_lo,
            M_star = M, M_star_min = max(M - M_lo, 0.0), M_star_max = M + M_lo,
            Teff = Teff, logg = logg, FeH = FeH,
            u = ldc_from_params(Teff, logg, FeH; passband = :TESS))
end

function _first_tic_row(data, tic::Int)
    if data isa JSON3.Object && haskey(data, :data)
        arr = data.data
        length(arr) == 0 && throw(CatalogLookupError("TIC $tic", 200, "no rows returned"))
        return arr[1]
    end
    return data
end

"""
Query Vizier's TAP service with an ADQL statement and return the parsed
JSON response. Caches by query string.
"""
function _vizier_tap_query(query::AbstractString; cache_key::AbstractString)
    cache_file = joinpath(@get_scratch!("vizier_cache"), cache_key * ".json")
    if isfile(cache_file)
        return JSON3.read(read(cache_file, String))
    end
    params = Dict("request" => "doQuery", "lang" => "adql",
                  "format" => "json", "query" => query)
    resp = HTTP.post(VIZIER_TAP_URL, [], HTTP.Form(params);
                     status_exception = false, retry = false,
                     readtimeout = 30)
    if resp.status != 200
        throw(CatalogLookupError("Vizier", resp.status, String(resp.body)))
    end
    body = String(resp.body)
    write(cache_file, body)
    return JSON3.read(body)
end

function _catalog_kic(kic::Int)
    # Kepler Stellar Properties Catalog (Mathur et al. 2017, J/ApJS/229/30/catalog)
    q = """SELECT "Teff","logg","[Fe/H]","Rad","Mass","e_Teff","e_Rad","e_Mass"
           FROM "J/ApJS/229/30/catalog" WHERE "KIC"=$kic"""
    row = try
        _first_tap_row(_vizier_tap_query(q; cache_key = "kic_$kic"))
    catch err
        if err isa CatalogLookupError
            @warn "KIC lookup failed; returning solar defaults" kic err
            return _default_stellar(:Kepler)
        else
            rethrow()
        end
    end
    row === nothing && return _default_stellar(:Kepler)
    R = _coalesce_float(get(row, :Rad, nothing), 1.0)
    Re = _coalesce_float(get(row, :e_Rad, nothing), 0.0)
    M = _coalesce_float(get(row, :Mass, nothing), 1.0)
    Me = _coalesce_float(get(row, :e_Mass, nothing), 0.0)
    Teff = _coalesce_float(get(row, :Teff, nothing), 5778.0)
    logg = _coalesce_float(get(row, :logg, nothing), 4.44)
    FeH  = _coalesce_float(get(row, Symbol("[Fe/H]"), nothing), 0.0)
    return (R_star = R, R_star_min = max(R - Re, 0.0), R_star_max = R + Re,
            M_star = M, M_star_min = max(M - Me, 0.0), M_star_max = M + Me,
            Teff = Teff, logg = logg, FeH = FeH,
            u = ldc_from_params(Teff, logg, FeH; passband = :Kepler))
end

function _catalog_epic(epic::Int)
    # EPIC catalog (Huber et al. 2016, IV/34/epic)
    q = """SELECT "Teff","logg","[Fe/H]","Rad","Mass","e_Rad","e_Mass"
           FROM "IV/34/epic" WHERE "ID"=$epic"""
    row = try
        _first_tap_row(_vizier_tap_query(q; cache_key = "epic_$epic"))
    catch err
        if err isa CatalogLookupError
            @warn "EPIC lookup failed; returning solar defaults" epic err
            return _default_stellar(:K2)
        else
            rethrow()
        end
    end
    row === nothing && return _default_stellar(:K2)
    R = _coalesce_float(get(row, :Rad, nothing), 1.0)
    Re = _coalesce_float(get(row, :e_Rad, nothing), 0.0)
    M = _coalesce_float(get(row, :Mass, nothing), 1.0)
    Me = _coalesce_float(get(row, :e_Mass, nothing), 0.0)
    Teff = _coalesce_float(get(row, :Teff, nothing), 5778.0)
    logg = _coalesce_float(get(row, :logg, nothing), 4.44)
    FeH  = _coalesce_float(get(row, Symbol("[Fe/H]"), nothing), 0.0)
    return (R_star = R, R_star_min = max(R - Re, 0.0), R_star_max = R + Re,
            M_star = M, M_star_min = max(M - Me, 0.0), M_star_max = M + Me,
            Teff = Teff, logg = logg, FeH = FeH,
            u = ldc_from_params(Teff, logg, FeH; passband = :K2))
end

function _first_tap_row(resp)
    # Vizier TAP JSON: { "metadata": [...], "data": [[...], [...]] }
    haskey(resp, :data) || return nothing
    data = resp.data
    length(data) == 0 && return nothing
    meta = resp.metadata
    names = [Symbol(m.name) for m in meta]
    row = data[1]
    return (; (names[i] => row[i] for i in eachindex(names))...)
end

_coalesce_float(x::Nothing, default) = Float64(default)
_coalesce_float(x::Missing, default) = Float64(default)
_coalesce_float(x, default) = try
    Float64(x)
catch
    Float64(default)
end

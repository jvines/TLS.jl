"""
    TLSOptions

Configuration for a TLS run. All fields have sensible defaults matching
the reference Python `transitleastsquares` package.
"""
Base.@kwdef struct TLSOptions
    R_star::Float64 = 1.0
    R_star_min::Float64 = 0.13
    R_star_max::Float64 = 3.5
    M_star::Float64 = 1.0
    M_star_min::Float64 = 0.1
    M_star_max::Float64 = 1.0
    u::Vector{Float64} = [0.4804, 0.1867]
    period_min::Union{Nothing,Float64} = nothing
    period_max::Union{Nothing,Float64} = nothing
    n_transits_min::Int = 2
    oversampling_factor::Int = 3
    duration_grid_step::Float64 = 1.1
    transit_template::Symbol = :default
    threads::Int = Threads.nthreads()
    verbose::Bool = false
    T0_fit_margin::Float64 = 0.01
end

"""
    TLSResult

Output of a TLS search. Field set is a superset of Python TLS so the
struct can be used as a drop-in replacement in downstream pipelines.
"""
Base.@kwdef mutable struct TLSResult
    periods::Vector{Float64} = Float64[]
    power::Vector{Float64} = Float64[]
    power_raw::Vector{Float64} = Float64[]
    SDE::Float64 = NaN
    SDE_raw::Float64 = NaN
    SR::Float64 = NaN
    chi2::Vector{Float64} = Float64[]
    chi2red::Vector{Float64} = Float64[]
    period::Float64 = NaN
    period_uncertainty::Float64 = NaN
    T0::Float64 = NaN
    duration::Float64 = NaN
    depth::Float64 = NaN
    depth_mean::Tuple{Float64,Float64} = (NaN, NaN)
    depth_mean_odd::Tuple{Float64,Float64} = (NaN, NaN)
    depth_mean_even::Tuple{Float64,Float64} = (NaN, NaN)
    odd_even_mismatch::Float64 = NaN
    transit_count::Int = 0
    distinct_transit_count::Int = 0
    empty_transit_count::Int = 0
    transit_depths::Vector{Float64} = Float64[]
    transit_depths_uncertainties::Vector{Float64} = Float64[]
    per_transit_count::Vector{Int} = Int[]
    snr::Float64 = NaN
    snr_per_transit::Vector{Float64} = Float64[]
    snr_pink_per_transit::Vector{Float64} = Float64[]
    rs_of_periods::Vector{Float64} = Float64[]
    model_lightcurve_time::Vector{Float64} = Float64[]
    model_lightcurve_model::Vector{Float64} = Float64[]
    model_folded_phase::Vector{Float64} = Float64[]
    model_folded_model::Vector{Float64} = Float64[]
    folded_phase::Vector{Float64} = Float64[]
    folded_y::Vector{Float64} = Float64[]
    folded_dy::Vector{Float64} = Float64[]
    in_transit_count::Int = 0
end

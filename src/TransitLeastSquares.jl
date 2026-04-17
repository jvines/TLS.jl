module TransitLeastSquares

using LinearAlgebra
using Statistics
using Printf
using PrecompileTools: @setup_workload, @compile_workload

export tls, TLSResult, TLSOptions, catalog_info

include("types.jl")
include("grid.jl")
include("duration.jl")
include("templates.jl")
include("detrend.jl")
include("search.jl")
include("statistics.jl")
include("stats_post.jl")
include("ldgrid.jl")       # must precede catalog.jl (provides ldc_from_params)
include("catalog.jl")
include("io.jl")

"""
    tls(time, flux; kwargs...) -> TLSResult

Run the Transit Least Squares period search (Hippke & Heller 2019) on a
light curve.

See [`TLSOptions`](@ref) for the full list of keyword arguments.
"""
function tls(time::AbstractVector{<:Real},
             flux::AbstractVector{<:Real};
             flux_err::Union{Nothing,AbstractVector{<:Real}} = nothing,
             kwargs...)
    opts = TLSOptions(; kwargs...)
    return _tls(collect(Float64, time),
                collect(Float64, flux),
                flux_err === nothing ? nothing : collect(Float64, flux_err),
                opts)
end

end # module TransitLeastSquares

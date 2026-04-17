using TLS
using Test
using Random
using Statistics

@testset "TLS.jl" verbose = true begin
    include("test_grid.jl")
    include("test_templates.jl")
    include("test_options.jl")
    include("test_ldgrid.jl")
    include("test_catalog.jl")
    include("test_spectra.jl")
    include("test_recovery.jl")
    include("test_performance.jl")
end

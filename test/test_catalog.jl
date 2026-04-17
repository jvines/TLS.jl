@testset "catalog_info argument validation" begin
    @test_throws ArgumentError catalog_info()
    @test_throws ArgumentError catalog_info(TIC = 1, KIC = 2)
end

@testset "CatalogLookupError shows cleanly" begin
    err = TLS.CatalogLookupError("TIC 1", 404, "not found")
    buf = IOBuffer()
    showerror(buf, err)
    s = String(take!(buf))
    @test occursin("TIC 1", s)
    @test occursin("404", s)
    @test occursin("not found", s)
end

@testset "_coalesce_float handles missing / nothing / bad values" begin
    @test TLS._coalesce_float(nothing, 1.5) == 1.5
    @test TLS._coalesce_float(missing, 1.5) == 1.5
    @test TLS._coalesce_float(2.25, 1.5) == 2.25
    @test TLS._coalesce_float("not a number", 1.5) == 1.5
end

# The actual HTTP calls to MAST / Vizier are not exercised in CI to keep
# tests deterministic and network-free. A tagged end-to-end check can be
# added under a `NETWORK_TESTS=1` environment flag in a later revision.

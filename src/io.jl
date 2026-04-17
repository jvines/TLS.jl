function Base.show(io::IO, ::MIME"text/plain", r::TLSResult)
    println(io, "TLSResult:")
    @printf(io, "  period     = %.6f d\n", r.period)
    @printf(io, "  T0         = %.6f d\n", r.T0)
    @printf(io, "  duration   = %.6f d\n", r.duration)
    @printf(io, "  depth      = %.6f\n", r.depth)
    @printf(io, "  SDE        = %.3f\n", r.SDE)
    @printf(io, "  SDE_raw    = %.3f\n", r.SDE_raw)
    @printf(io, "  SR         = %.3f\n", r.SR)
    @printf(io, "  snr        = %.3f\n", r.snr)
    @printf(io, "  transits   = %d (%d distinct, %d empty)\n",
            r.transit_count, r.distinct_transit_count, r.empty_transit_count)
    @printf(io, "  grid size  = %d periods\n", length(r.periods))
end

Base.show(io::IO, r::TLSResult) =
    @printf(io, "TLSResult(period=%.6f, SDE=%.2f)", r.period, r.SDE)

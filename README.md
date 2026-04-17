# TLS.jl

A Julia implementation of the [Transit Least Squares](https://arxiv.org/abs/1901.02015)
algorithm (Hippke & Heller 2019) for detecting periodic transits in
photometric time series. Designed to be feature-compatible with the
reference Python package
[`transitleastsquares`](https://github.com/hippke/tls) while taking
advantage of Julia's multi-threading and SIMD for substantially higher
throughput.

## Installation

```julia
using Pkg
Pkg.add("TLS")      # once registered in the Julia General registry
```

## Quick start

```julia
using TLS

# time in days, flux normalized to ~1
result = tls(time, flux; flux_err = flux_err)

result.period        # best-fit period (days)
result.T0            # transit mid-time of the first transit (days)
result.duration      # transit duration (days)
result.depth         # fractional transit depth
result.SDE           # signal detection efficiency
result.power         # detrended power spectrum vs result.periods
```

### Catalog lookup

If you have a TIC / KIC / EPIC identifier, `catalog_info` fetches
stellar parameters and limb-darkening coefficients:

```julia
info = catalog_info(TIC = 307210830)
result = tls(time, flux;
             R_star = info.R_star, M_star = info.M_star, u = info.u)
```

## Parallelism

The period search parallelises over the period grid. Launch Julia with
the number of threads you want to use:

```bash
julia -t auto --project -e 'using TLS; ...'
```

## Benchmarks

Head-to-head against the reference Python
[`transitleastsquares`](https://github.com/hippke/tls) (v1.32) on an
identical 5000-point synthetic light curve (injected
P=3.5 d, T0=1.2 d, duration=0.15 d, depth=0.01, white noise σ=5e-4,
60-day baseline; grid covers 1.0-8.0 d at oversampling=2, yielding
~2380 trial periods; run on a 28-core Linux box):

| Metric            | Python TLS (1T) | Python TLS (4T) | TLS.jl (1T) | TLS.jl (4T) |
|-------------------|----------------:|----------------:|------------:|------------:|
| Wall time         |       10.46 s   |        4.75 s   |    5.05 s   |  **1.41 s** |
| Period (true 3.5) |       3.500254  |       3.500254  |   3.498740  |   3.498740  |
| T0 (true 1.2)     |       1.197572  |       1.197572  |   1.209165  |   1.209165  |
| Duration (true 0.15) |    0.167294  |       0.167294  |   0.153819  |   0.153819  |
| Depth (true 0.01) |       0.010603  |       0.010603  |   0.010751  |   0.010751  |
| SDE               |       39.485    |       39.485    |    33.929   |    33.929   |
| Transit count     |             17  |             17  |         17  |         17  |
| Period grid size  |           2380  |           2380  |       2377  |       2377  |

At matched thread count TLS.jl is ~3× faster; against Python's default
single-thread config it is ~7× faster. All recovered transit
parameters agree with Python TLS to within the grid resolution. The
residual SDE gap (33.9 vs 39.5) comes from the underlying χ² values —
TLS.jl currently uses an analytic small-planet Mandel-Agol template
rather than the full `batman` model Python TLS uses; the
SR/power/SDE **statistic itself is normalized identically** (see
`src/statistics.jl`).

Reproduce with:

```
julia --project benchmark/generate_lc.jl
python3 benchmark/run_python.py
julia --project -t 4 benchmark/run_julia.jl
```

## Citation

If you use TLS.jl in published work, please cite Hippke & Heller (2019):

```bibtex
@article{Hippke2019,
  author  = {Hippke, Michael and Heller, Ren\'e},
  title   = {{Optimized transit detection algorithm to search for periodic transits of small planets}},
  journal = {Astronomy & Astrophysics},
  volume  = {623},
  pages   = {A39},
  year    = {2019},
  doi     = {10.1051/0004-6361/201834672}
}
```

## Status

Early development (v0.1). See `test/runtests.jl` for the validated
feature set.

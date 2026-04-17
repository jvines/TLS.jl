#!/usr/bin/env python3
"""Run Python transitleastsquares on the shared lc.csv and print
recovered parameters + wall-clock time."""
import sys, types, time, os

# distutils shim for Python 3.12 (batman-package imports distutils).
d = types.ModuleType('distutils'); c = types.ModuleType('distutils.ccompiler')
class _C:
    def has_function(self, *a, **k): return False
c.new_compiler = lambda *a, **k: _C()
d.ccompiler = c
sys.modules['distutils'] = d
sys.modules['distutils.ccompiler'] = c

import numpy as np
from transitleastsquares import transitleastsquares

here = os.path.dirname(os.path.abspath(__file__))
data = np.loadtxt(os.path.join(here, "lc.csv"), delimiter=",")
t, f = data[:, 0], data[:, 1]
print(f"loaded {len(t)} points; span={t[-1]-t[0]:.2f}d; median flux={np.median(f):.6f}")

n_threads = int(os.environ.get("TLS_THREADS", "1"))
model = transitleastsquares(t, f)
t0 = time.perf_counter()
r = model.power(
    period_min=1.0, period_max=8.0,
    oversampling_factor=2,
    duration_grid_step=1.1,
    n_transits_min=2,
    show_progress_bar=False,
    use_threads=n_threads,
)
elapsed = time.perf_counter() - t0
print(f"(ran with use_threads={n_threads})")

print(f"\n=== Python TLS ===")
print(f"elapsed       : {elapsed:.3f} s")
print(f"period        : {r.period:.6f} d  (true 3.5)")
print(f"T0            : {r.T0:.6f} d  (true 1.2 mod P)")
print(f"duration      : {r.duration:.6f} d  (true 0.15)")
print(f"depth         : {r.depth:.6f}")
print(f"SDE           : {r.SDE:.3f}")
try:
    print(f"SR            : {float(np.max(r.SR)):.6f}")
except Exception:
    pass
print(f"transit_count : {r.transit_count}")
print(f"periods_grid  : {len(r.periods)} points")

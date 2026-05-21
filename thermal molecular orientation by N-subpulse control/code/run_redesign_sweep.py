"""
Production thermally-aware-redesign sweep.

For each (Jmax, T) case:
  - Compute baseline (Hong et al.) thermal peak
  - Optimise the N pulse areas via Nelder-Mead
  - Record both, save to pickle
"""
import sys
sys.path.insert(0, '/home/claude/work/extension')
import numpy as np
import pickle
import time
from scipy.optimize import minimize

from thermal_simulation import (
    design_pulse_M0, T_rot, B,
)
from thermally_aware_redesign import thermal_peak
from schur_horn_bound import schur_horn_bound

# Fast grid
Nt = 30000
t = np.linspace(-15*T_rot, 225*T_rot, Nt)
dt = t[1] - t[0]
J_max_basis = 12
M_max = 6
print(f"Nt={Nt}, dt={dt:.0f} a.u. = {dt/T_rot:.4f} T_rot", file=sys.stderr)
print(f"J_max_basis={J_max_basis}, M_max={M_max}", file=sys.stderr)

# Cases (Jmax, T_K, maxiter) — keep maxiter modest for time budget
cases = [
    (5,  10.0,  60),
    (10, 10.0,  80),
    (10, 20.0,  80),
]

results = {}
out_path = '/home/claude/work/extension/optimised_results.pkl'
for Jmax, T_K, maxit in cases:
    print(f"\n=== Jmax={Jmax}, T={T_K} K, maxiter={maxit} ===", file=sys.stderr, flush=True)
    _, theta0, lam = design_pulse_M0(Jmax, t)

    pk0 = thermal_peak(theta0, T_K, t, dt,
                       J_max_basis=J_max_basis, M_max=M_max)
    bound = schur_horn_bound(T_K, Jmax)
    print(f"  baseline {pk0:.4f} | bound {bound:.4f} | "
          f"baseline/bound {pk0/bound:.3f}", file=sys.stderr, flush=True)

    def neg(x):
        return -thermal_peak(x, T_K, t, dt,
                             J_max_basis=J_max_basis, M_max=M_max)
    t0 = time.time()
    res = minimize(neg, theta0, method='Nelder-Mead',
                   options={'xatol': 5e-3, 'fatol': 1e-4,
                            'maxiter': maxit, 'adaptive': True})
    elapsed = time.time() - t0
    pk_opt = -res.fun
    print(f"  optimised {pk_opt:.4f} | optimised/bound {pk_opt/bound:.3f} | "
          f"Delta = {pk_opt - pk0:+.4f} ({(pk_opt-pk0)/pk0*100:+.1f}%) | "
          f"nfev={res.nfev} | {elapsed:.0f} s", file=sys.stderr, flush=True)

    results[(Jmax, T_K)] = {
        'lambda':    lam,
        'theta0':    theta0,
        'peak0':     pk0,
        'theta_opt': res.x,
        'peak_opt':  pk_opt,
        'bound':     bound,
        'nfev':      int(res.nfev),
        'elapsed':   elapsed,
    }
    # Save after EACH case in case we time out
    with open(out_path, 'wb') as f:
        pickle.dump(results, f)
    print(f"  saved after {len(results)} cases", file=sys.stderr, flush=True)

# Summary
print("\n" + "=" * 84, file=sys.stderr)
print("THERMALLY-AWARE REDESIGN SUMMARY", file=sys.stderr)
print("=" * 84, file=sys.stderr)
print(f"{'Jmax':>4} | {'T(K)':>5} | {'bound':>6} | {'baseline':>9} | "
      f"{'optimised':>9} | {'Delta':>7} | {'bd ratio'}", file=sys.stderr)
for (Jmax, T_K), r in results.items():
    delta = r['peak_opt'] - r['peak0']
    pct = delta / r['peak0'] * 100
    print(f"{Jmax:>4} | {T_K:>5.1f} | {r['bound']:>6.4f} | {r['peak0']:>9.4f} | "
          f"{r['peak_opt']:>9.4f} | {delta:>+7.4f} | "
          f"{r['peak0']/r['bound']:.3f} -> {r['peak_opt']/r['bound']:.3f}",
          file=sys.stderr)

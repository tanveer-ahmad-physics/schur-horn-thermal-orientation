"""
Thermally-aware redesign of the Hong et al. (2025) N-subpulse protocol.

We keep the Hong et al. pulse structure fixed (N subpulses, Gaussian rms width
T_p = 3 T_rot, central times tau_n = 5(n-1) T_p, carrier frequencies
omega_n = 2 n B, carrier phases phi_n = pi/2 + omega_n tau_n) and optimise the
N pulse areas {theta_n} to maximise the thermally averaged peak orientation
max_t <<cos theta>>_T(t).

The T=0 protocol of Hong et al. is recovered as the special case T -> 0; we
seed the optimiser with the Hong et al. areas and run Nelder-Mead.

Output:
  optimised_results.pkl  -- per-(Jmax, T_K), the optimised areas, the
                            optimised peak, the Hong et al. baseline peak,
                            and the Schur-Horn bound.
"""
import numpy as np
import pickle
import time
import sys
from scipy.optimize import minimize

sys.path.insert(0, '/home/claude/work/extension')
from thermal_simulation import (
    B, mu0, T_rot, Msector_ops, boltzmann_weights, design_pulse_M0,
)


# -------------------------------------------------------------------------
# Build E(t) from a list of pulse areas
# -------------------------------------------------------------------------
def build_E_from_areas(areas, t, Tn=None):
    """Build E(t) from {theta_n}, using the Hong et al. layout (rest fixed)."""
    if Tn is None:
        Tn = 3 * T_rot
    N = len(areas)
    if N == 0:
        return np.zeros_like(t)
    tau   = np.arange(N) * 5 * Tn
    omega = 2 * np.arange(1, N + 1) * B
    phi   = np.pi/2 + omega * tau
    E = np.zeros_like(t)
    for n in range(N):
        n_p   = n + 1
        mu_n  = mu0 * n_p / np.sqrt((2*n_p + 1) * (2*n_p - 1))
        F0_n  = np.sqrt(2/np.pi) * areas[n] / (mu_n * Tn)
        env   = np.exp(-((t - tau[n])**2) / (2 * Tn**2))
        car   = np.cos(omega[n] * (t - tau[n]) + phi[n])
        E    += F0_n * env * car
    return E


# -------------------------------------------------------------------------
# Thermal peak: simulate ensemble and return max <<cos theta>>_T(t)
# -------------------------------------------------------------------------
def propagate_M_sector_only_peak(E_field, dt, M, J_max_basis, mu0_val,
                                  t_arr, t_window):
    """Propagate one M-sector and return ONLY the cos-theta trace for
    times inside t_window (faster: skips observable computation when not
    needed)."""
    H0_diag, cos_M, Jlist = Msector_ops(M, J_max_basis)
    n = len(Jlist)
    MU_M = mu0_val * cos_M
    Nt = len(E_field)
    D_mu, V_mu = np.linalg.eigh(MU_M)
    phase_half = np.exp(-1j * H0_diag * dt/2)
    Psi = np.eye(n, dtype=complex)
    E_half = 0.5 * (E_field[:-1] + E_field[1:])
    t_Trot = t_arr / T_rot
    # Slice indices for the window
    win_idx = np.where((t_Trot >= t_window[0]) & (t_Trot <= t_window[1]))[0]
    cos_track = np.zeros((len(win_idx), n))
    out_idx = 0
    for j in range(1, Nt):
        Psi = phase_half[:, None] * Psi
        Psi = V_mu.conj().T @ Psi
        Psi = np.exp(1j * D_mu * E_half[j-1] * dt)[:, None] * Psi
        Psi = V_mu @ Psi
        Psi = phase_half[:, None] * Psi
        if out_idx < len(win_idx) and j == win_idx[out_idx]:
            aux = cos_M @ Psi
            cos_track[out_idx, :] = np.real(np.sum(Psi.conj() * aux, axis=0))
            out_idx += 1
    return cos_track, Jlist, win_idx


def thermal_peak(areas, T_K, t_arr, dt, J_max_basis=15, M_max=10,
                 t_window=(222.4, 225.0)):
    """Max_t <<cos theta>>_T(t) for the given pulse areas."""
    E = build_E_from_areas(areas, t_arr)
    weights = boltzmann_weights(T_K, J_max_basis)
    total = None
    for M in range(0, M_max + 1):
        cos_t_M, Jlist, win_idx = propagate_M_sector_only_peak(
            E, dt, M, J_max_basis, mu0, t_arr, t_window)
        w = weights[M]
        mult = 1 if M == 0 else 2
        contrib = (cos_t_M @ w) * mult
        if total is None:
            total = contrib
        else:
            total += contrib
    return float(total.max())


# -------------------------------------------------------------------------
# Optimisation driver
# -------------------------------------------------------------------------
def optimise(Jmax, T_K, t_arr, dt, J_max_basis=15, M_max=10,
             maxiter=200, verbose=True):
    """Optimise pulse areas to maximise thermal peak orientation.

    Returns dict with keys:
        areas0     -- Hong et al. baseline pulse areas
        peak0      -- baseline thermal peak
        areas_opt  -- optimised pulse areas
        peak_opt   -- optimised thermal peak
        success    -- scipy optimiser flag
        nfev       -- number of function evaluations
    """
    # Hong et al. baseline pulse areas
    _, theta0, lam = design_pulse_M0(Jmax, t_arr)
    peak0 = thermal_peak(theta0, T_K, t_arr, dt,
                         J_max_basis=J_max_basis, M_max=M_max)
    if verbose:
        print(f"  Jmax={Jmax}, T={T_K} K: baseline peak = {peak0:.4f}",
              file=sys.stderr)

    def neg_peak(x):
        return -thermal_peak(x, T_K, t_arr, dt,
                             J_max_basis=J_max_basis, M_max=M_max)

    t0 = time.time()
    res = minimize(neg_peak, theta0, method='Nelder-Mead',
                   options={'xatol': 1e-3, 'fatol': 1e-5,
                            'maxiter': maxiter, 'adaptive': True})
    elapsed = time.time() - t0
    peak_opt = -res.fun
    if verbose:
        print(f"  Jmax={Jmax}, T={T_K} K: optimised peak = {peak_opt:.4f} "
              f"  (Delta = {peak_opt - peak0:+.4f}, nfev={res.nfev}, "
              f"{elapsed:.0f} s)", file=sys.stderr)

    return {
        'areas0':    theta0,
        'peak0':     peak0,
        'areas_opt': res.x,
        'peak_opt':  peak_opt,
        'success':   bool(res.success),
        'nfev':      int(res.nfev),
        'elapsed':   elapsed,
    }


# -------------------------------------------------------------------------
# Main: run for representative (Jmax, T) cases
# -------------------------------------------------------------------------
if __name__ == "__main__":
    # Use a coarser time grid for the optimiser to keep runtime modest
    Nt = 80000
    t  = np.linspace(-15*T_rot, 225*T_rot, Nt)
    dt = t[1] - t[0]
    print(f"Nt={Nt}, dt={dt:.2f} a.u. = {dt/T_rot:.5f} T_rot",
          file=sys.stderr)

    # Cross-check: at this Nt, the baseline at T=0 should still recover
    # lambda within ~ 1e-4
    for Jmax in [5, 10]:
        _, theta0, lam = design_pulse_M0(Jmax, t)
        pk0 = thermal_peak(theta0, 0.0, t, dt, J_max_basis=15, M_max=10)
        print(f"  baseline check: Jmax={Jmax} lam={lam:.4f}, "
              f"sim peak T=0 = {pk0:.4f}, |diff|={abs(pk0 - lam):.2e}",
              file=sys.stderr)

    # Cases of interest:
    #   - Jmax=5,  T=10 K (moderate Jmax, beyond B/kB)
    #   - Jmax=10, T=10 K (high Jmax, moderate T)
    #   - Jmax=10, T=20 K (high Jmax, high T)
    cases = [
        (5,  10.0),
        (10, 10.0),
        (10, 20.0),
    ]
    results = {}
    for Jmax, T_K in cases:
        res = optimise(Jmax, T_K, t, dt,
                       J_max_basis=15, M_max=10, maxiter=300)
        results[(Jmax, T_K)] = res

    out = '/home/claude/work/extension/optimised_results.pkl'
    with open(out, 'wb') as f:
        pickle.dump(results, f)
    print(f"\nSaved {out}", file=sys.stderr)

    # Summary
    print("\n" + "=" * 78, file=sys.stderr)
    print("SUMMARY OF THERMALLY-AWARE REDESIGN",
          file=sys.stderr)
    print("=" * 78, file=sys.stderr)
    print(f"{'Jmax':>4} | {'T(K)':>5} | {'Hong baseline':>14} | "
          f"{'optimised':>10} | {'improvement':>11}", file=sys.stderr)
    for (Jmax, T_K), r in results.items():
        delta = r['peak_opt'] - r['peak0']
        rel   = delta / r['peak0'] * 100
        print(f"{Jmax:>4} | {T_K:>5.1f} | {r['peak0']:>14.4f} | "
              f"{r['peak_opt']:>10.4f} | {delta:+.4f} ({rel:+.1f}%)",
              file=sys.stderr)

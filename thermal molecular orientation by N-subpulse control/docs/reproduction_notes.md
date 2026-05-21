# Reproduction Notes

This file documents exactly what is reproducible from this minimal repository, what is not, and how to fill the gap.

## What is reproducible from this repository alone

- **Reading Table 2 of the manuscript.** The pickle `data/optimised_results.pkl` is the direct output of `code/run_redesign_sweep.py` and contains, for each $(J_{\max}, T)$ case studied at $T = 10$ K, the baseline thermal peak, the optimised pulse areas, the optimised thermal peak, and the analytical Schur–Horn bound. Open it with `pickle.load()` (see Quickstart in the main README).

- **Visual inspection of Figures 1–4.** The PNG files in `figures/` are exactly as published in the manuscript.

## What is not bundled in this repository

The two Python files in `code/` import two helper modules that are part of the author's larger working directory and are not redistributed here:

| Module | Purpose | Manuscript reference |
|---|---|---|
| `thermal_simulation.py` | Full $M$-sector thermal propagator: builds the dipole-coupled Hamiltonian, integrates the time-dependent Schrödinger equation with a second-order Strang split-operator scheme in matrix form over all initial $\|J, M\rangle$ states in parallel, returns the thermally averaged orientation trace $\langle\!\langle\cos\theta\rangle\!\rangle_T(t)$. Also exposes the LiH constants `B`, `mu0`, `T_rot`. | Section 2 (Hamiltonian + pulse-train construction) + Appendix A (numerical implementation). |
| `schur_horn_bound.py` | Analytical Schur–Horn upper bound: builds the $\cos\theta$ matrix per $M$-sector, diagonalises, sorts both eigenvalue lists in decreasing order, contracts with the sorted Boltzmann eigenvalue list, sums across sectors. | Equations (8)–(10) of the manuscript. |

Each module is approximately 100 lines of NumPy. The mathematical recipe is fully specified in the manuscript text, so reconstruction is straightforward.

## Why these were not bundled

Two reasons:

1. **Scope.** The Zenodo / GitHub release is intended as a minimal, archivable companion to the manuscript: the optimisation driver, the Schur–Horn bound output, and the figures. The full working directory contains additional scratch files, experimental variants, and intermediate development artefacts that would inflate the deposit without aiding reproducibility.

2. **Author preference.** The full simulation directory is available on request — see Contact in the main README. The author is happy to share it for verification, extension, or collaboration.

## How to fill the gap

If you would like a self-contained, runnable version of this repository:

**Option A — request the full directory.** Email `tanveer.quantum@gmail.com` or [open an issue](../../issues) and the author will share the complete working directory.

**Option B — reconstruct the two modules from the manuscript.** The minimal API surface that `thermally_aware_redesign.py` and `run_redesign_sweep.py` expect:

```python
# In thermal_simulation.py:
B       # LiH rotational constant in atomic units
mu0     # LiH dipole moment in atomic units
T_rot   # rotational period = pi / B

def design_pulse_M0(Jmax, t):
    """Return (E_field_trace, pulse_areas, eigenvalue)
       for the Hong et al. M=0 analytical protocol."""

def Msector_ops(Jmax, M_max):
    """Return the (cos theta)_M matrices and H_0 per M-sector."""

def boltzmann_weights(T_K, J_max_basis, M_max):
    """Return the Boltzmann weights per (J, M)."""
```

```python
# In schur_horn_bound.py:
def schur_horn_bound(T_K, Jmax, J_max_basis=None, M_max=None):
    """Return the Schur–Horn upper bound on thermal orientation
       for LiH at temperature T_K and pulse-design subspace Jmax."""
```

With these signatures filled in (each function is ~10–30 lines), `python run_redesign_sweep.py` will reproduce `data/optimised_results.pkl`.

## Reproducibility of Table 1 (verification suite)

Table 1 of the manuscript reports eight quantitative cross-checks of the analytical protocol's reproduction. These tests are not bundled here in scripted form, but the eight checks are listed verbatim in Table 1 itself and are straightforward to script once `thermal_simulation.py` is rebuilt:

1. Eigenvalue from characteristic polynomial vs. `numpy.linalg.eigh`.
2. Eigenvector norm and sign convention.
3. Recurrence-based pulse areas vs. numerical-eigenvector pulse areas.
4. Integrated pulse-area vs. designed pulse-area.
5. Peak-orientation timing vs. integer rotational periods.
6. Leakage out of the designed subspace.
7. $\langle\cos^2\theta\rangle$ peak vs. analytical eigenvector formula.
8. $\langle\cos\theta\rangle$ peak vs. eigenvalue.

Each test compares one numerical quantity to one analytical quantity and reports the maximum discrepancy across $J_{\max} \in \{1, \ldots, 15\}$. See Table 1 of the manuscript for the values.

## Provenance of `data/optimised_results.pkl`

| Case | $J_{\max}$ | $T$ (K) | Nelder–Mead iterations |
|---|---:|---:|---:|
| 1 | 5  | 10 | 60 |
| 2 | 10 | 10 | 80 |
| 3 | 10 | 20 | 80 |

Stored fields per case:
- `theta_baseline` — analytical pulse areas (from Eq. (5) of the manuscript)
- `peak_baseline` — thermal peak from the analytical baseline
- `theta_opt` — Nelder–Mead optimised pulse areas
- `peak_opt` — thermal peak from the optimised areas
- `bound` — analytical Schur–Horn upper bound
- `gap_closed` — $(p_\text{opt} - p_\text{base}) / (p_\text{bound} - p_\text{base})$

These map directly to the columns of Table 2 of the manuscript.

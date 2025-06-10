#!/usr/bin/env python3
"""
Standalone λ₁/λ₂ stability scan for the complete 4-variable Adler-2020 model.

Usage example
-------------
python lambda_scan_standalone.py \
       --lambda1_lo .2 --lambda1_hi 2 \
       --lambda2_lo .2 --lambda2_hi 2 \
       --N 4 --nproc 10   > scan.csv
"""

from __future__ import annotations
import argparse, csv, itertools, sys, multiprocessing as mp
from typing import Any, Dict, List, Tuple
import numpy as np
from scipy.optimize  import root
from scipy.optimize import least_squares, brentq
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

def format_row(row):
    return [f"{x:.3}" if isinstance(x, float) and not np.isnan(x) else "nan" if np.isnan(x) else x for x in row]

def phi_F(F, M, p):
    P = P_qss(F, M, p)
    return p["lambda1"] * P / (p["k1"] + P) * (1 - F/p["K"]) - p["mu1"]

def phi_M(F, M, p):
    C = C_qss(F, M, p)
    return p["lambda2"] * C / (p["k2"] + C) - p["mu2"]

def safe_brentq(f, lo, hi, *args):
    """Return brentq root or np.nan if the bracket never changes sign."""
    flo, fhi = f(lo), f(hi)
    if flo == 0.0:                 # exact root at lo
        return lo
    if flo * fhi > 0.0:            # same sign → try to expand hi
        for factor in (1e13, 1e15, 1e18):   # three expansions are enough
            try_hi = hi * factor
            if f(lo) * f(try_hi) < 0.0:
                return brentq(f, lo, try_hi, *args)
        return np.nan              # give up – let caller skip this node
    return brentq(f, lo, hi, *args)


def P_qss(F, M, p):
    f = lambda P: p["beta2"]*M + p["beta3"]*F \
                  - p["alpha2"]*F*P/(p["k1"]+P) - p["gamma"]*P
    return safe_brentq(f, 0.0, 1e12)


def C_qss(F, M, p):
    f = lambda C: p["beta1"]*F \
                  - p["alpha1"]*M*C/(p["k2"]+C) - p["gamma"]*C
    return safe_brentq(f, 0.0, 1e12)

def g_F(F: float, M: float, p: Dict[str, float]) -> float:
    """
    Reduced RHS for fibroblasts.
    Uses the quasi-steady-state ligand level P_qss(F, M).
    """
    P = P_qss(F, M, p)
    growth_term = p["lambda1"] * P / (p["k1"] + P) * (1.0 - F / p["K"])
    return (growth_term - p["mu1"]) * F


def g_M(F: float, M: float, p: Dict[str, float]) -> float:
    """
    Reduced RHS for macrophages.
    Uses the quasi-steady-state ligand level C_qss(F, M).
    """
    C = C_qss(F, M, p)
    growth_term = p["lambda2"] * C / (p["k2"] + C)
    return (growth_term - p["mu2"]) * M

def plot_heat(grid_rows, which: str, title: str, fig_name: str):
    """
    Draw a heat-map of F* for the chosen sink (‘cold’ or ‘hot’) and overlay
    M* contour lines.
    """
    # unpack grid → 2-D arrays
    l1_vals = sorted(set(r[0] for r in grid_rows))
    l2_vals = sorted(set(r[1] for r in grid_rows))
    Fmap = np.full((len(l2_vals), len(l1_vals)), np.nan)
    Mmap = np.full_like(Fmap, np.nan)

    for row in grid_rows:
        i = l2_vals.index(row[1])
        j = l1_vals.index(row[0])
        if which == "cold":
            M = row[4]; F = row[5]
        else:  # hot
            M = row[6]; F = row[7]
        Fmap[i, j] = F
        Mmap[i, j] = M

    positive_vals = Fmap[Fmap > 0]

    if positive_vals.size == 0:
        print("⚠️ No positive values in Fmap — skipping log color scale.")
        return  # or fall back to linear plotting
    else:
        fig, ax = plt.subplots()
        h = ax.imshow(Fmap, origin="lower",
                      extent=[min(l1_vals), max(l1_vals),
                              min(l2_vals), max(l2_vals)],
                      aspect="auto",
                      norm=LogNorm(vmin=np.nanmin(positive_vals),
                                   vmax=np.nanmax(positive_vals)))

        cbar = plt.colorbar(h, ax=ax)
        cbar.set_label("F* (cells)")
        CS = ax.contour(l1_vals, l2_vals, Mmap, colors="k", linewidths=0.6)
        ax.clabel(CS, inline=True, fmt="%.1e", fontsize=7)
        ax.set_xlabel("λ₁  (day⁻¹)")
        ax.set_xscale("log")
        ax.set_ylabel("λ₂  (day⁻¹)")
        ax.set_yscale("log")
        ax.set_title(title)
        plt.tight_layout()
        fig.tight_layout()
        fig.savefig(fig_name, dpi=300, bbox_inches="tight")
        print(f"Saved → {fig_name}")


# ─────────────────────────── parameters & utilities ──────────────────────
DEFAULT_PARS: Dict[str, float] = dict(
    # *** FILL THESE WITH THE NUMBERS FROM YOUR ORIGINAL FILE -------------
    # kinetics
    mu1   = 0.30,     # myofibroblast apoptosis rate  (1/d)
    mu2   = 0.30,     # macrophage apoptosis rate     (1/d)
    beta1 = 470 * 1440,   # CSF-1 secretion  (molecules/(cell·d))
    beta2 = 70 * 1440,   # PDGF secretion  (molecules/(cell·d))
    beta3 = 240 * 1440,
    alpha1= 940 * 1440,   # CSF-1 endocytosis (molecules d⁻¹)
    alpha2= 510 * 1440,   # PDGF endocytosis
    gamma = 2.0,     # ligand clearance (1/d)
    k1    = 1e9,      # PDGF half-sat. (molecules/ml)
    k2    = 1e9,      # CSF-1 half-sat.
    K     = 1e6,    # fibroblast carrying capacity (cells/ml)
    # *** the two parameters we scan:
    lambda1 = 0.8,    # PDGF autocrine strength (d⁻¹)
    lambda2 = 0.8,    # CSF-1 macrophage feedback (d⁻¹)
)

def rhs(t: float, y: np.ndarray, p: Dict[str, float]) -> np.ndarray:
    """Full Adler ODE system:  y = [M,F,P,C]."""
    M, F, P, C = y
    dF = (p["lambda1"]*P/(p["k1"]+P)*(1 - F/p["K"]) - p["mu1"]) * F
    dM = (p["lambda2"]*C/(p["k2"]+C)               - p["mu2"]) * M
    dP =  p["beta2"]*M + p["beta3"]*F \
        - p["alpha2"]*F*P/(p["k1"]+P) - p["gamma"]*P
    dC =  p["beta1"]*F \
        - p["alpha1"]*M*C/(p["k2"]+C) - p["gamma"]*C
    return np.array([dM, dF, dP, dC])

def jacobian(y: np.ndarray, p: Dict[str, float]) -> np.ndarray:
    """Analytic 4×4 Jacobian at state y=[M,F,P,C]."""
    M, F, P, C = y
    k1, k2 = p["k1"], p["k2"]
    mu1, mu2 = p["mu1"], p["mu2"]
    la1, la2 = p["lambda1"], p["lambda2"]
    K        = p["K"]
    # helpers
    s1 = P/(k1+P);  s2 = C/(k2+C)
    ds1 = k1 / (k1+P)**2
    ds2 = k2 / (k2+C)**2
    J = np.zeros((4,4))
    # dM/d• rows
    J[0,0] = la2*s2 - mu2
    J[0,3] = la2*M*ds2
    # dF/d• rows
    J[1,1] = la1*s1*(1 - 2*F/p["K"]) - mu1
    J[1,2] = la1*F*ds1*(1 - F/K)
    # dP/d• rows
    J[2,0] = p["beta2"]
    J[2,1] = p["beta3"] - p["alpha2"]*P/(k1+P)
    J[2,2] = -p["alpha2"]*F*(1 - P/(k1+P)) - p["gamma"]
    # dC/d• rows
    J[3,0] = -p["alpha1"]*C/(k2+C)
    J[3,1] = p["beta1"]
    J[3,3] = -p["alpha1"]*M*(1 - C/(k2+C)) - p["gamma"]
    return J

# ───────────────────── fixed-point search utilities ──────────────────────
def newton_polish(seed: np.ndarray, p: Dict[str, float],
                  tol: float = 1e-10, itmax: int = 1000000) -> np.ndarray|None:
    sol = root(lambda x: rhs(0.0, x, p),
               seed, jac=lambda x: jacobian(x, p),
               method="hybr",
               options={"xtol": tol, "maxfev": itmax})
    return sol.x if sol.success else None

def automatic_seeds(p: Dict[str, float]) -> List[np.ndarray]:
    """
    Return up to three seeds [M,F,P,C]:
        – healing (always)
        – cold fibrosis (if it exists)
        – hot fibrosis  (if it exists)
    Now uses triangle-based sign detection with barycentric seed placement.
    """

    # ── 0. quick test: pure-healing region ──────────────────────────────
    if p["lambda1"] < p["mu1"] and p["lambda2"] < p["mu2"]:
        return [np.array([10., 10.,
                          P_qss(10., 10., p),
                          C_qss(10., 10., p)])]

    # ── 1. sign grid for φ_F and φ_M ────────────────────────────────────
    Ms = np.r_[0.0, np.logspace(1, 7, 160)]  # include 0 explicitly
    Fs = np.logspace(5, 6, 160)
    signs_F = np.empty((len(Ms), len(Fs)), int)
    signs_M = np.empty_like(signs_F)

    for i, M in enumerate(Ms):
        for j, F in enumerate(Fs):
            signs_F[i, j] = np.sign(phi_F(F, M, p))
            signs_M[i, j] = np.sign(phi_M(F, M, p))

    # ── 2. always include a small healing seed ──────────────────────────
    seeds = [np.array([50., 50.,
                       P_qss(50., 50., p),
                       C_qss(50., 50., p)])]

    # ── 3. scan triangles for sign flips ────────────────────────────────
    for i in range(len(Ms) - 1):
        for j in range(len(Fs) - 1):

            triangles = [
                [(i, j), (i+1, j), (i+1, j+1)],  # lower-right triangle
                [(i, j), (i, j+1), (i+1, j+1)]   # upper-left triangle
            ]

            for tri in triangles:
                sF = {signs_F[x][y] for x, y in tri}
                sM = {signs_M[x][y] for x, y in tri}

                if {-1, 1}.issubset(sF) and {-1, 1}.issubset(sM):
                    # Barycentric: (2/3, 1/3) → slightly biased away from vertex
                    M0 = (2/3)*Ms[tri[0][0]] + (1/3)*Ms[tri[2][0]]
                    F0 = (2/3)*Fs[tri[0][1]] + (1/3)*Fs[tri[2][1]]
                    seeds.append(np.array([M0, F0,
                                           P_qss(F0, M0, p),
                                           C_qss(F0, M0, p)]))
                    if len(seeds) >= 3:
                        return seeds

    return seeds


def find_fp_triplet(p: Dict[str, float]) -> Tuple[np.ndarray|None,
                                                  np.ndarray|None,
                                                  np.ndarray|None]:
    """
    Return (heal, cold, hot) fixed points or None where not found.

    Strategy: try a small set of seeds, Newton-polish, deduplicate,
              then sort by fibroblast number F.
    """

    seeds = automatic_seeds(p)
    fps: List[np.ndarray] = []

    def _is_same(a, b, tol=0.05):
        return np.linalg.norm(np.log10(np.maximum(1e-30, a[:2]) /
                                       np.maximum(1e-30, b[:2]))) < tol

    for sd in seeds:
        fp = newton_polish(sd, p, tol=1e-8)
        if fp is None: continue
        if np.linalg.norm(rhs(0, fp, p)) / (p["beta3"] * p["K"]) > 1e-4:
            continue
        print("fp:", fp)
        print("rhs(fp):", rhs(0.0, fp, p))
        print("‖rhs‖:", np.linalg.norm(rhs(0.0, fp, p)))
        if np.linalg.norm(rhs(0.0, fp, p)) > 30000000000:
            print(f"❌ Spurios root found: {sd}")
            continue  # reject spurious root
        print(f"✅ Converged: {sd} → F={fp[1]:.2}, M={fp[0]:.2}")
        # deduplicate (distance in log-space < 5 %)
        if not any(_is_same(fp, q) for q in fps):
            fps.append(fp)

    if not fps:
        return None, None, None
    # sort by F (second component) → heal (lowest F) … hot (highest F)
    fps.sort(key=lambda v: v[1])
    heal = fps[0]
    cold = fps[1] if len(fps) > 2 else None
    hot  = fps[-1] if len(fps) > 1  else None
    return heal, cold, hot

# ─────────────────────────── grid-scan helper ────────────────────────────
def eig_triplet_for(l1: float, l2: float) -> Tuple[Any, ...]:
    """Return a row with λ₁, λ₂, (healing: M, F), (cold: M, F), (hot: M, F)."""
    p = DEFAULT_PARS.copy()
    p["lambda1"] = l1
    p["lambda2"] = l2
    heal, cold, hot = find_fp_triplet(p)

    def pack(fp):
        if fp is None:
            return np.nan, np.nan
        return fp[0], fp[1]      # eig, M, F

    Mh, Fh = pack(heal)
    Mc, Fc = pack(cold)
    Mhhot, Fhhot = pack(hot)

    return (l1, l2,
            Mh,  Fh,
            Mc,  Fc,
            Mhhot, Fhhot)
# ────────────────────────────────── main ─────────────────────────────────
def main() -> None:
    pa = argparse.ArgumentParser("λ₁/λ₂ stability scan – standalone model")
    pa.add_argument("--lambda1_lo", type=float, default=0.05)
    pa.add_argument("--lambda1_hi", type=float, default=2.0)
    pa.add_argument("--lambda2_lo", type=float, default=0.05)
    pa.add_argument("--lambda2_hi", type=float, default=2.0)
    pa.add_argument("--fig_name", type=str, default='hot.png')
    pa.add_argument("--N", type=int, default=20,   help="grid points per axis")
    pa.add_argument("--nproc", type=int, default=0, help="processes (0→serial)")
    args = pa.parse_args()

    l1_vals = np.linspace(args.lambda1_lo, args.lambda1_hi, args.N)
    l2_vals = np.linspace(args.lambda2_lo, args.lambda2_hi, args.N)
    grid = list(itertools.product(l1_vals, l2_vals))

    if args.nproc > 1:
        with mp.Pool(args.nproc) as pool:
            rows = pool.starmap(eig_triplet_for, grid)
    else:
        rows = [eig_triplet_for(l1, l2) for l1, l2 in grid]

    rows.sort(key=lambda r: (r[0], r[1]))
    writer = csv.writer(sys.stdout)
    writer.writerow([
        "lambda1", "lambda2",
        "M_heal", "F_heal",
        "M_cold", "F_cold",
        "M_hot", "F_hot"
    ])
    writer.writerows(map(format_row, rows))
    plot_heat(rows, "hot", "heat map: F*, contour lines: M*", args.fig_name)
if __name__ == "__main__":
    main()


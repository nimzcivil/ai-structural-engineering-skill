# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "numpy>=1.24",
# ]
# ///
"""
Column Buckling Design per Eurocode 3 (EN 1993-1-1).

Computes:
- Euler critical load (N_cr)
- Slenderness ratio (λ)
- Non-dimensional slenderness (λ̄)
- Buckling reduction factor (χ) per clause 6.3.1.2
- Buckling resistance (N_b,Rd)

All inputs in SI: metres, Pascals.
"""

import argparse
import json
import math
import sys


# Imperfection factors per EN 1993-1-1 Table 6.1
ALPHA_IMPERFECTION = {
    "a0": 0.13,
    "a": 0.21,
    "b": 0.34,
    "c": 0.49,
    "d": 0.76,
}


def column_buckling(L: float, A_mm2: float, I_min_mm4: float, fy: float,
                     E: float, curve: str, k: float,
                     gamma_M1: float = 1.0) -> dict:
    """
    EN 1993-1-1 clause 6.3.1 flexural buckling.

    Parameters
    ----------
    L : float       – Column length (m)
    A_mm2 : float   – Cross-section area (mm²)
    I_min_mm4: float– Minimum second moment of area (mm⁴)
    fy : float      – Yield strength (MPa)
    E : float       – Young's modulus (MPa)
    curve : str     – Buckling curve: a0, a, b, c, d
    k : float       – Effective length factor
    """
    if curve not in ALPHA_IMPERFECTION:
        return {"error": f"Unknown buckling curve '{curve}'. Use: a0, a, b, c, d"}

    alpha = ALPHA_IMPERFECTION[curve]

    # Convert everything to base SI for computation
    A = A_mm2 * 1e-6       # m²
    I_min = I_min_mm4 * 1e-12  # m⁴
    fy_Pa = fy * 1e6       # Pa
    E_Pa = E * 1e6         # Pa

    # Effective length
    L_cr = k * L

    # Radius of gyration
    i = math.sqrt(I_min / A)

    # Slenderness ratio
    lam = L_cr / i

    # Euler critical load
    N_cr = math.pi**2 * E_Pa * I_min / L_cr**2

    # Squash load
    N_pl = A * fy_Pa

    # Non-dimensional slenderness
    lam_bar = math.sqrt(N_pl / N_cr)

    # Φ factor
    phi = 0.5 * (1 + alpha * (lam_bar - 0.2) + lam_bar**2)

    # Reduction factor χ
    chi = 1 / (phi + math.sqrt(phi**2 - lam_bar**2))
    chi = min(chi, 1.0)

    # Buckling resistance
    N_b_Rd = chi * A * fy_Pa / gamma_M1

    return {
        "column_design": "EN 1993-1-1 cl.6.3.1",
        "inputs": {
            "L_m": L,
            "k": k,
            "L_cr_m": round(L_cr, 4),
            "A_mm2": A_mm2,
            "I_min_mm4": I_min_mm4,
            "fy_MPa": fy,
            "E_MPa": E,
            "buckling_curve": curve,
            "alpha_imperfection": alpha,
            "gamma_M1": gamma_M1,
        },
        "results": {
            "i_mm": round(i * 1000, 2),
            "slenderness_lambda": round(lam, 2),
            "N_cr_kN": round(N_cr / 1e3, 2),
            "N_pl_kN": round(N_pl / 1e3, 2),
            "lambda_bar": round(lam_bar, 4),
            "Phi": round(phi, 4),
            "chi": round(chi, 4),
            "N_b_Rd_kN": round(N_b_Rd / 1e3, 2),
        },
        "checks": {
            "slender_column": lam_bar > 0.2,
            "stocky_no_buckling_check_needed": lam_bar <= 0.2,
        },
    }

def main():
    parser = argparse.ArgumentParser(description="Column buckling design — EN 1993-1-1")
    parser.add_argument("--length", type=float, required=True, help="Column length (m)")
    parser.add_argument("--A", type=float, required=True, help="Cross-section area (mm²)")
    parser.add_argument("--I_min", type=float, required=True, help="Min moment of inertia (mm⁴)")
    parser.add_argument("--fy", type=float, required=True, help="Yield strength (MPa)")
    parser.add_argument("--E", type=float, default=210000, help="Young's modulus (MPa)")
    parser.add_argument("--buckling_curve", required=True, choices=ALPHA_IMPERFECTION.keys())
    parser.add_argument("--k", type=float, default=1.0, help="Effective length factor")
    parser.add_argument("--gamma_M1", type=float, default=1.0, help="Partial safety factor γ_M1")

    args = parser.parse_args()

    if args.length <= 0 or args.A <= 0 or args.I_min <= 0 or args.fy <= 0:
        print(json.dumps({"error": "All dimensions and material properties must be positive"}))
        sys.exit(1)

    result = column_buckling(
        L=args.length, A_mm2=args.A, I_min_mm4=args.I_min, fy=args.fy,
        E=args.E, curve=args.buckling_curve, k=args.k, gamma_M1=args.gamma_M1,
    )
    print(json.dumps(result, indent=2))

if __name__ == "__main__":
    main()

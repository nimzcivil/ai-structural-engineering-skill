# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "numpy>=1.24",
#     "sympy>=1.12",
# ]
# ///
"""
Beam Analysis — Simply Supported, Cantilever, Fixed-Fixed, Propped Cantilever.

Computes reactions, shear force, bending moment, and deflection at n points
along the beam using symbolic (exact) superposition via SymPy, then evaluates
numerically with NumPy for the output arrays.

All inputs in SI: metres, Newtons, Pascals.
"""

import argparse
import json
import sys

import numpy as np
import sympy as sp


# ──────────────────────── helpers ────────────────────────

def _parse_loads(loads_json: str | None):
    if not loads_json:
        return []
    return json.loads(loads_json)


def _parse_distributed(dist_json: str | None):
    if not dist_json:
        return []
    return json.loads(dist_json)


# ──────────────────────── simply supported beam ─────────

def _simply_supported(L, E, I, point_loads, dist_loads, n=200):
    """Superposition for a simply-supported beam of span L."""
    x_vals = np.linspace(0, L, n)
    V = np.zeros(n)
    M = np.zeros(n)
    delta = np.zeros(n)

    R_A_total = 0.0
    R_B_total = 0.0

    # --- point loads ---
    for ld in point_loads:
        if ld["type"] == "point":
            P = ld["P"]
            a = ld["a"]
            b = L - a
            R_A = P * b / L
            R_B = P * a / L
            R_A_total += R_A
            R_B_total += R_B
            for i, x in enumerate(x_vals):
                if x <= a:
                    V[i] += R_A
                    M[i] += R_A * x
                    delta[i] += (P * b * x) / (6 * L * E * I) * (L**2 - b**2 - x**2)
                else:
                    V[i] += R_A - P
                    M[i] += R_A * x - P * (x - a)
                    delta[i] += (P * a * (L - x)) / (6 * L * E * I) * (2 * L * (L - x) - a**2 - (L - x)**2)
        elif ld["type"] == "moment":
            M0 = ld["M"]
            a = ld["a"]
            R_A = -M0 / L
            R_B = M0 / L
            R_A_total += R_A
            R_B_total += R_B
            for i, x in enumerate(x_vals):
                if x <= a:
                    V[i] += R_A
                    M[i] += R_A * x
                else:
                    V[i] += R_A
                    M[i] += R_A * x + M0

    # --- distributed loads ---
    for dl in dist_loads:
        w = dl["w"]
        a = dl.get("a", 0)
        b = dl.get("b", L)
        ld_len = b - a
        W_total = w * ld_len
        # Resultant at centroid of loaded segment
        centroid = a + ld_len / 2
        R_A = W_total * (L - centroid) / L
        R_B = W_total * centroid / L
        R_A_total += R_A
        R_B_total += R_B
        for i, x in enumerate(x_vals):
            # contribution from R_A
            v_i = R_A
            m_i = R_A * x
            # subtract load that has been passed
            if x > a:
                x_in_load = min(x, b) - a
                v_i -= w * x_in_load
                m_i -= w * x_in_load * (x - a - x_in_load / 2)
            V[i] += v_i
            M[i] += m_i

        # deflection for full-span UDL shortcut
        if abs(a) < 1e-12 and abs(b - L) < 1e-12:
            for i, x in enumerate(x_vals):
                delta[i] += (w * x) / (24 * E * I) * (L**3 - 2 * L * x**2 + x**3)
        else:
            # numerical integration fallback for partial UDL
            for i, x in enumerate(x_vals):
                # Macaulay-style deflection approximation via virtual work
                # Use moment-area or direct integration
                pass  # handled below by sympy if needed
            delta += _deflection_sympy_ss(L, E, I, w, a, b, x_vals)

    max_V = float(np.max(np.abs(V)))
    max_M = float(np.max(np.abs(M)))
    max_delta = float(np.max(np.abs(delta)))
    i_max_M = int(np.argmax(np.abs(M)))

    return {
        "beam_type": "simply_supported",
        "span_m": L,
        "reactions": {"R_A_N": round(R_A_total, 2), "R_B_N": round(R_B_total, 2)},
        "max_shear_N": round(max_V, 2),
        "max_moment_Nm": round(max_M, 2),
        "max_moment_location_m": round(float(x_vals[i_max_M]), 4),
        "max_deflection_m": round(max_delta, 6),
        "max_deflection_mm": round(max_delta * 1000, 3),
        "span_over_deflection": round(L / max_delta, 1) if max_delta > 1e-15 else float("inf"),
        "shear_diagram": {"x_m": x_vals.tolist(), "V_N": V.tolist()},
        "moment_diagram": {"x_m": x_vals.tolist(), "M_Nm": M.tolist()},
        "deflection_curve": {"x_m": x_vals.tolist(), "delta_m": delta.tolist()},
    }


def _deflection_sympy_ss(L, E, I, w, a, b, x_vals):
    """Symbolic deflection for partial UDL on simply-supported beam."""
    x_sym = sp.Symbol("x", positive=True)
    # moment expression from R_A and partial UDL
    ld_len = b - a
    W = w * ld_len
    centroid = a + ld_len / 2
    R_A = W * (L - centroid) / L

    # M(x) piecewise
    M_expr = sp.Piecewise(
        (R_A * x_sym, x_sym <= a),
        (R_A * x_sym - w * (x_sym - a)**2 / 2, (x_sym > a) & (x_sym <= b)),
        (R_A * x_sym - W * (x_sym - centroid), x_sym > b),
    )

    # EI * y'' = -M  =>  integrate twice
    theta = sp.integrate(-M_expr / (E * I), x_sym)
    y = sp.integrate(theta, x_sym)

    # BCs: y(0) = 0, y(L) = 0
    C1, C2 = sp.symbols("C1 C2")
    y = y + C1 * x_sym + C2
    sol = sp.solve([y.subs(x_sym, 0), y.subs(x_sym, L)], [C1, C2])
    y = y.subs(sol)
    y_func = sp.lambdify(x_sym, y, modules=["numpy"])

    return np.array([float(y_func(xv)) for xv in x_vals])


# ──────────────────────── cantilever beam ────────────────

def _cantilever(L, E, I, point_loads, dist_loads, n=200):
    """Fixed at x=0, free at x=L."""
    x_vals = np.linspace(0, L, n)
    V = np.zeros(n)
    M = np.zeros(n)
    delta = np.zeros(n)

    R_A_total = 0.0
    M_A_total = 0.0

    for ld in point_loads:
        if ld["type"] == "point":
            P = ld["P"]
            a = ld["a"]
            R_A_total += P
            M_A_total += P * a
            for i, x in enumerate(x_vals):
                if x <= a:
                    V[i] += P
                    M[i] += -P * (a - x)
                    delta[i] += (P * x**2) / (6 * E * I) * (3 * a - x)
                else:
                    V[i] += 0
                    M[i] += 0
                    delta[i] += (P * a**2) / (6 * E * I) * (3 * x - a)

    for dl in dist_loads:
        w = dl["w"]
        a = dl.get("a", 0)
        b = dl.get("b", L)
        R_A_total += w * (b - a)
        M_A_total += w * (b - a) * (a + (b - a) / 2)

        if abs(a) < 1e-12 and abs(b - L) < 1e-12:
            for i, x in enumerate(x_vals):
                V[i] += w * (L - x)
                M[i] += -w * (L - x)**2 / 2
                delta[i] += (w * x**2) / (24 * E * I) * (6 * L**2 - 4 * L * x + x**2)
        else:
            for i, x in enumerate(x_vals):
                if x <= a:
                    V[i] += w * (b - a)
                    M[i] += -w * (b - a) * ((a + (b - a) / 2) - x)
                elif x <= b:
                    V[i] += w * (b - x)
                    M[i] += -w * (b - x)**2 / 2
                # deflection via sympy below
            delta += _deflection_sympy_cantilever(L, E, I, w, a, b, x_vals)

    max_V = float(np.max(np.abs(V)))
    max_M = float(np.max(np.abs(M)))
    max_delta = float(np.max(np.abs(delta)))

    return {
        "beam_type": "cantilever",
        "span_m": L,
        "reactions": {"R_A_N": round(R_A_total, 2), "M_A_Nm": round(M_A_total, 2)},
        "max_shear_N": round(max_V, 2),
        "max_moment_Nm": round(max_M, 2),
        "max_deflection_m": round(max_delta, 6),
        "max_deflection_mm": round(max_delta * 1000, 3),
        "tip_deflection_mm": round(abs(float(delta[-1])) * 1000, 3),
        "shear_diagram": {"x_m": x_vals.tolist(), "V_N": V.tolist()},
        "moment_diagram": {"x_m": x_vals.tolist(), "M_Nm": M.tolist()},
        "deflection_curve": {"x_m": x_vals.tolist(), "delta_m": delta.tolist()},
    }


def _deflection_sympy_cantilever(L, E, I, w, a, b, x_vals):
    """Symbolic deflection for partial UDL on cantilever (fixed at x=0)."""
    x_sym = sp.Symbol("x", positive=True)
    ld_len = b - a
    W = w * ld_len
    centroid = a + ld_len / 2

    M_expr = sp.Piecewise(
        (-W * (centroid - x_sym), x_sym <= a),
        (-w * (b - x_sym)**2 / 2, (x_sym > a) & (x_sym <= b)),
        (sp.Integer(0), x_sym > b),
    )

    theta = sp.integrate(-M_expr / (E * I), x_sym)
    y = sp.integrate(theta, x_sym)

    C1, C2 = sp.symbols("C1 C2")
    y = y + C1 * x_sym + C2
    theta_full = theta + C1
    sol = sp.solve([y.subs(x_sym, 0), theta_full.subs(x_sym, 0)], [C1, C2])
    y = y.subs(sol)
    y_func = sp.lambdify(x_sym, y, modules=["numpy"])

    return np.array([float(y_func(xv)) for xv in x_vals])


# ──────────────────────── fixed-fixed beam ──────────────

def _fixed_fixed(L, E, I, point_loads, dist_loads, n=200):
    """Both ends fixed. Uses SymPy for exact solution."""
    x_sym = sp.Symbol("x", positive=True)
    x_vals = np.linspace(0, L, n)

    R_A, R_B, M_A, M_B = sp.symbols("R_A R_B M_A M_B")

    # Build load expressions using Macaulay brackets via Piecewise
    V_expr = R_A
    M_expr = M_A + R_A * x_sym

    for ld in point_loads:
        if ld["type"] == "point":
            P = ld["P"]
            a = ld["a"]
            V_expr += -P * sp.Piecewise((1, x_sym >= a), (0, True))
            M_expr += -P * sp.Piecewise((x_sym - a, x_sym >= a), (0, True))

    for dl in dist_loads:
        w = dl["w"]
        a_pos = dl.get("a", 0)
        b_pos = dl.get("b", L)
        V_expr += -w * sp.Piecewise((x_sym - a_pos, (x_sym >= a_pos) & (x_sym <= b_pos)), (b_pos - a_pos, x_sym > b_pos), (0, True))
        M_expr += -w * sp.Piecewise(((x_sym - a_pos)**2 / 2, (x_sym >= a_pos) & (x_sym <= b_pos)),
                                     ((b_pos - a_pos) * (x_sym - (a_pos + b_pos) / 2), x_sym > b_pos),
                                     (0, True))

    # EI y'' = M => integrate twice with BCs y(0)=0, y'(0)=0, y(L)=0, y'(L)=0
    theta_expr = sp.integrate(M_expr / (E * I), x_sym)
    y_expr = sp.integrate(theta_expr, x_sym)

    # BCs
    eqs = [
        y_expr.subs(x_sym, 0),       # y(0)=0
        theta_expr.subs(x_sym, 0),   # y'(0)=0
        y_expr.subs(x_sym, L),       # y(L)=0
        theta_expr.subs(x_sym, L),   # y'(L)=0
    ]

    # Equilibrium: R_A + R_B = total vertical load
    total_P = sum(ld["P"] for ld in point_loads if ld["type"] == "point")
    total_W = sum(dl["w"] * (dl.get("b", L) - dl.get("a", 0)) for dl in dist_loads)

    sol = sp.solve(eqs + [R_A + R_B - total_P - total_W], [R_A, R_B, M_A, M_B], dict=True)
    if not sol:
        return {"error": "Could not solve fixed-fixed beam. Check inputs."}
    sol = sol[0]

    V_num = sp.lambdify(x_sym, V_expr.subs(sol), modules=["numpy"])
    M_num = sp.lambdify(x_sym, M_expr.subs(sol), modules=["numpy"])
    y_num = sp.lambdify(x_sym, y_expr.subs(sol), modules=["numpy"])

    V_arr = np.array([float(V_num(xv)) for xv in x_vals])
    M_arr = np.array([float(M_num(xv)) for xv in x_vals])
    d_arr = np.array([float(y_num(xv)) for xv in x_vals])

    return {
        "beam_type": "fixed_fixed",
        "span_m": L,
        "reactions": {
            "R_A_N": round(float(sol[R_A]), 2),
            "R_B_N": round(float(sol[R_B]), 2),
            "M_A_Nm": round(float(sol[M_A]), 2),
            "M_B_Nm": round(float(sol[M_B]), 2),
        },
        "max_shear_N": round(float(np.max(np.abs(V_arr))), 2),
        "max_moment_Nm": round(float(np.max(np.abs(M_arr))), 2),
        "max_deflection_m": round(float(np.max(np.abs(d_arr))), 6),
        "max_deflection_mm": round(float(np.max(np.abs(d_arr))) * 1000, 3),
        "shear_diagram": {"x_m": x_vals.tolist(), "V_N": V_arr.tolist()},
        "moment_diagram": {"x_m": x_vals.tolist(), "M_Nm": M_arr.tolist()},
        "deflection_curve": {"x_m": x_vals.tolist(), "delta_m": d_arr.tolist()},
    }


# ──────────────────────── propped cantilever ────────────

def _propped_cantilever(L, E, I, point_loads, dist_loads, n=200):
    """Fixed at x=0, simply-supported at x=L (roller)."""
    x_sym = sp.Symbol("x", positive=True)
    x_vals = np.linspace(0, L, n)

    R_A, R_B, M_A = sp.symbols("R_A R_B M_A")

    V_expr = R_A
    M_expr = M_A + R_A * x_sym

    for ld in point_loads:
        if ld["type"] == "point":
            P = ld["P"]
            a = ld["a"]
            V_expr += -P * sp.Piecewise((1, x_sym >= a), (0, True))
            M_expr += -P * sp.Piecewise((x_sym - a, x_sym >= a), (0, True))

    for dl in dist_loads:
        w = dl["w"]
        a_pos = dl.get("a", 0)
        b_pos = dl.get("b", L)
        V_expr += -w * sp.Piecewise((x_sym - a_pos, (x_sym >= a_pos) & (x_sym <= b_pos)), (b_pos - a_pos, x_sym > b_pos), (0, True))
        M_expr += -w * sp.Piecewise(((x_sym - a_pos)**2 / 2, (x_sym >= a_pos) & (x_sym <= b_pos)),
                                     ((b_pos - a_pos) * (x_sym - (a_pos + b_pos) / 2), x_sym > b_pos),
                                     (0, True))

    theta_expr = sp.integrate(M_expr / (E * I), x_sym)
    y_expr = sp.integrate(theta_expr, x_sym)

    total_P = sum(ld["P"] for ld in point_loads if ld["type"] == "point")
    total_W = sum(dl["w"] * (dl.get("b", L) - dl.get("a", 0)) for dl in dist_loads)

    # BCs: y(0)=0, y'(0)=0, y(L)=0, equilibrium ΣFy=0
    eqs = [
        y_expr.subs(x_sym, 0),
        theta_expr.subs(x_sym, 0),
        y_expr.subs(x_sym, L),
        R_A + R_B - total_P - total_W,
    ]

    sol = sp.solve(eqs, [R_A, R_B, M_A], dict=True)
    if not sol:
        return {"error": "Could not solve propped cantilever. Check inputs."}
    sol = sol[0]

    V_num = sp.lambdify(x_sym, V_expr.subs(sol), modules=["numpy"])
    M_num = sp.lambdify(x_sym, M_expr.subs(sol), modules=["numpy"])
    y_num = sp.lambdify(x_sym, y_expr.subs(sol), modules=["numpy"])

    V_arr = np.array([float(V_num(xv)) for xv in x_vals])
    M_arr = np.array([float(M_num(xv)) for xv in x_vals])
    d_arr = np.array([float(y_num(xv)) for xv in x_vals])

    return {
        "beam_type": "propped_cantilever",
        "span_m": L,
        "reactions": {
            "R_A_N": round(float(sol[R_A]), 2),
            "R_B_N": round(float(sol[R_B]), 2),
            "M_A_Nm": round(float(sol[M_A]), 2),
        },
        "max_shear_N": round(float(np.max(np.abs(V_arr))), 2),
        "max_moment_Nm": round(float(np.max(np.abs(M_arr))), 2),
        "max_deflection_m": round(float(np.max(np.abs(d_arr))), 6),
        "max_deflection_mm": round(float(np.max(np.abs(d_arr))) * 1000, 3),
        "shear_diagram": {"x_m": x_vals.tolist(), "V_N": V_arr.tolist()},
        "moment_diagram": {"x_m": x_vals.tolist(), "M_Nm": M_arr.tolist()},
        "deflection_curve": {"x_m": x_vals.tolist(), "delta_m": d_arr.tolist()},
    }


# ──────────────────────── main ──────────────────────────

BEAM_FUNCS = {
    "simply_supported": _simply_supported,
    "cantilever": _cantilever,
    "fixed_fixed": _fixed_fixed,
    "propped_cantilever": _propped_cantilever,
}


def main():
    parser = argparse.ArgumentParser(description="Beam analysis (SI units)")
    parser.add_argument("--type", required=True, choices=BEAM_FUNCS.keys())
    parser.add_argument("--length", type=float, required=True, help="Beam span (m)")
    parser.add_argument("--E", type=float, default=210e9, help="Young's modulus (Pa)")
    parser.add_argument("--I", type=float, required=True, help="Moment of inertia (m⁴)")
    parser.add_argument("--loads", type=str, default=None, help="JSON array of point loads")
    parser.add_argument("--distributed", type=str, default=None, help="JSON array of distributed loads")
    parser.add_argument("--n", type=int, default=200, help="Number of evaluation points")

    args = parser.parse_args()

    if args.length <= 0:
        print(json.dumps({"error": "Beam length must be positive"}))
        sys.exit(1)
    if args.I <= 0:
        print(json.dumps({"error": "Moment of inertia must be positive"}))
        sys.exit(1)

    point_loads = _parse_loads(args.loads)
    dist_loads = _parse_distributed(args.distributed)

    if not point_loads and not dist_loads:
        print(json.dumps({"error": "At least one load must be specified (--loads or --distributed)"}))
        sys.exit(1)

    func = BEAM_FUNCS[args.type]
    result = func(args.length, args.E, args.I, point_loads, dist_loads, args.n)

    # Strip diagram arrays for cleaner console output (keep summary)
    output = {k: v for k, v in result.items() if k not in ("shear_diagram", "moment_diagram", "deflection_curve")}
    output["note"] = "Full diagram arrays computed. Key results shown above."

    print(json.dumps(output, indent=2, default=str))


if __name__ == "__main__":
    main()

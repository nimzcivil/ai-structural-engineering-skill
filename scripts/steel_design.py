# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "numpy>=1.24",
# ]
# ///
"""
Steel Member Design per Eurocode 3 (EN 1993-1-1).

Modules:
1. Cross-section classification (cl. 5.5)
2. Moment resistance Mc,Rd (cl. 6.2.5)
3. Shear resistance Vpl,Rd (cl. 6.2.6)
4. Lateral-torsional buckling Mb,Rd (cl. 6.3.2)

Inputs: dimensions in mm, strengths in MPa.
"""

import argparse
import json
import math
import sys


# ──────────────────── section properties (reuse logic) ──

def _I_section_props(params):
    h = params["h"]
    b = params["b"]
    tf = params["tf"]
    tw = params["tw"]
    hw = h - 2 * tf

    A = 2 * b * tf + hw * tw
    Iy = (b * h**3 - (b - tw) * hw**3) / 12
    Iz = (2 * tf * b**3 + hw * tw**3) / 12
    Wel_y = Iy / (h / 2)
    Wpl_y = 2 * b * tf * (h / 2 - tf / 2) + tw * hw**2 / 4
    Wel_z = Iz / (b / 2)
    Wpl_z = 2 * (tf * b**2 / 4) + hw * tw**2 / 4
    Iw = (tf * b**3 * (h - tf)**2) / 24
    It = (2 * b * tf**3 + hw * tw**3) / 3
    Av = max(A - 2 * b * tf + (tw + 2 * 0) * tf, hw * tw)  # shear area (approx)

    return {
        "h": h, "b": b, "tf": tf, "tw": tw, "hw": hw,
        "A": A, "Iy": Iy, "Iz": Iz,
        "Wel_y": Wel_y, "Wpl_y": Wpl_y,
        "Wel_z": Wel_z, "Wpl_z": Wpl_z,
        "Iw": Iw, "It": It, "Av": Av,
    }


def _hollow_rect_props(params):
    b_o = params["b"]
    h_o = params["h"]
    t = params["t"]
    b_i = b_o - 2 * t
    h_i = h_o - 2 * t
    A = b_o * h_o - b_i * h_i
    Iy = (b_o * h_o**3 - b_i * h_i**3) / 12
    Iz = (h_o * b_o**3 - h_i * b_i**3) / 12
    Wel_y = Iy / (h_o / 2)
    Wpl_y = b_o * h_o**2 / 4 - b_i * h_i**2 / 4
    Av = 2 * h_o * t
    It = 2 * t * (b_o - t)**2 * (h_o - t)**2 / (b_o + h_o - 2 * t)
    return {
        "h": h_o, "b": b_o, "tf": t, "tw": t, "hw": h_i,
        "A": A, "Iy": Iy, "Iz": Iz,
        "Wel_y": Wel_y, "Wpl_y": Wpl_y,
        "Wel_z": Iz / (b_o / 2), "Wpl_z": h_o * b_o**2 / 4 - h_i * b_i**2 / 4,
        "Iw": 0, "It": It, "Av": Av,
    }


SHAPE_FUNCS = {
    "I_section": _I_section_props,
    "hollow_rect": _hollow_rect_props,
}


# ──────────────────── classification ───────────────────

def classify_section(props, fy):
    """
    EN 1993-1-1 Table 5.2 — simplified classification for I-sections.
    """
    epsilon = math.sqrt(235 / fy)
    hw = props["hw"]
    tw = props["tw"]
    b = props["b"]
    tf = props["tf"]

    # Outstand flange (one side): c = (b - tw) / 2
    c_flange = (b - tw) / 2
    cf_tf = c_flange / tf

    # Web in bending: c = hw
    cw_tw = hw / tw

    # Flange class
    if cf_tf <= 9 * epsilon:
        flange_class = 1
    elif cf_tf <= 10 * epsilon:
        flange_class = 2
    elif cf_tf <= 14 * epsilon:
        flange_class = 3
    else:
        flange_class = 4

    # Web class (pure bending)
    if cw_tw <= 72 * epsilon:
        web_class = 1
    elif cw_tw <= 83 * epsilon:
        web_class = 2
    elif cw_tw <= 124 * epsilon:
        web_class = 3
    else:
        web_class = 4

    section_class = max(flange_class, web_class)

    return {
        "epsilon": round(epsilon, 4),
        "c_flange_mm": round(c_flange, 1),
        "c_flange_over_tf": round(cf_tf, 2),
        "c_web_over_tw": round(cw_tw, 2),
        "flange_class": flange_class,
        "web_class": web_class,
        "section_class": section_class,
    }


# ──────────────────── moment resistance ────────────────

def moment_resistance(props, fy, gamma_M0=1.0):
    """Mc,Rd per EN 1993-1-1 cl. 6.2.5."""
    cls = classify_section(props, fy)
    section_class = cls["section_class"]

    if section_class <= 2:
        Mc_Rd = props["Wpl_y"] * fy / gamma_M0  # N·mm
    elif section_class == 3:
        Mc_Rd = props["Wel_y"] * fy / gamma_M0
    else:
        # Class 4 — effective section needed (simplified: use Wel as conservative)
        Mc_Rd = props["Wel_y"] * fy / gamma_M0

    return {
        "calc": "moment_resistance",
        "code": "EN 1993-1-1 cl.6.2.5",
        "classification": cls,
        "fy_MPa": fy,
        "gamma_M0": gamma_M0,
        "Wpl_y_mm3": round(props["Wpl_y"], 1),
        "Wel_y_mm3": round(props["Wel_y"], 1),
        "Mc_Rd_kNm": round(Mc_Rd / 1e6, 2),
    }


# ──────────────────── shear resistance ─────────────────

def shear_resistance(props, fy, gamma_M0=1.0):
    """Vpl,Rd per EN 1993-1-1 cl. 6.2.6."""
    Av = props["Av"]
    Vpl_Rd = Av * (fy / math.sqrt(3)) / gamma_M0  # N

    return {
        "calc": "shear_resistance",
        "code": "EN 1993-1-1 cl.6.2.6",
        "Av_mm2": round(Av, 1),
        "fy_MPa": fy,
        "gamma_M0": gamma_M0,
        "Vpl_Rd_kN": round(Vpl_Rd / 1e3, 2),
    }


# ──────────────────── lateral-torsional buckling ───────

def ltb_resistance(props, fy, Lcr, C1=1.0, gamma_M1=1.0):
    """
    Mb,Rd per EN 1993-1-1 cl. 6.3.2.

    Lcr: lateral buckling length (mm)
    C1: moment distribution factor (1.0 for uniform moment)
    """
    E = 210000  # MPa
    G = 80769   # MPa

    Iz = props["Iz"]
    It = props["It"]
    Iw = props["Iw"]
    Wpl_y = props["Wpl_y"]
    Wel_y = props["Wel_y"]

    cls = classify_section(props, fy)
    section_class = cls["section_class"]

    # Elastic critical moment for LTB (cl. 6.3.2.2)
    # Mcr = C1 * π²EIz / Lcr² * √(Iw/Iz + Lcr²GIt / (π²EIz))
    pi2_E_Iz = math.pi**2 * E * Iz
    Mcr = C1 * pi2_E_Iz / Lcr**2 * math.sqrt(Iw / Iz + Lcr**2 * G * It / pi2_E_Iz)

    if section_class <= 2:
        W_y = Wpl_y
    else:
        W_y = Wel_y

    # Non-dimensional slenderness
    lam_LT = math.sqrt(W_y * fy / Mcr)

    # Imperfection factor (general case, rolled I-sections)
    # Table 6.3 / 6.4: for rolled I h/b ≤ 2 → curve a (αLT = 0.21)
    #                     h/b > 2  → curve b (αLT = 0.34)
    h_over_b = props["h"] / props["b"]
    if h_over_b <= 2:
        alpha_LT = 0.21  # curve a
        curve = "a"
    else:
        alpha_LT = 0.34  # curve b
        curve = "b"

    # General method (cl. 6.3.2.2)
    lam_LT_0 = 0.2  # plateau length (recommended)
    beta = 1.0       # recommended

    Phi_LT = 0.5 * (1 + alpha_LT * (lam_LT - lam_LT_0) + beta * lam_LT**2)
    chi_LT = 1 / (Phi_LT + math.sqrt(Phi_LT**2 - beta * lam_LT**2))
    chi_LT = min(chi_LT, 1.0)
    # Also cap: chi_LT ≤ 1/λ_LT²
    if lam_LT > 0:
        chi_LT = min(chi_LT, 1 / lam_LT**2)

    Mb_Rd = chi_LT * W_y * fy / gamma_M1  # N·mm

    return {
        "calc": "lateral_torsional_buckling",
        "code": "EN 1993-1-1 cl.6.3.2",
        "classification": cls,
        "Lcr_mm": Lcr,
        "C1": C1,
        "fy_MPa": fy,
        "gamma_M1": gamma_M1,
        "Iz_mm4": round(Iz, 1),
        "It_mm4": round(It, 1),
        "Iw_mm6": round(Iw, 1),
        "Mcr_kNm": round(Mcr / 1e6, 2),
        "W_y_mm3": round(W_y, 1),
        "lambda_LT": round(lam_LT, 4),
        "h_over_b": round(h_over_b, 2),
        "LTB_curve": curve,
        "alpha_LT": alpha_LT,
        "Phi_LT": round(Phi_LT, 4),
        "chi_LT": round(chi_LT, 4),
        "Mb_Rd_kNm": round(Mb_Rd / 1e6, 2),
    }


# ──────────────────── main ─────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Steel design — EN 1993-1-1")
    parser.add_argument("--calc", required=True, choices=["classify", "moment", "shear", "ltb"])
    parser.add_argument("--shape", required=True, choices=SHAPE_FUNCS.keys())
    parser.add_argument("--fy", type=float, required=True, help="Yield strength (MPa)")
    parser.add_argument("--gamma_M0", type=float, default=1.0)
    parser.add_argument("--gamma_M1", type=float, default=1.0)
    parser.add_argument("--Lcr", type=float, default=0, help="Lateral buckling length (mm)")
    parser.add_argument("--C1", type=float, default=1.0, help="Moment distribution factor")

    args, unknown = parser.parse_known_args()

    shape_params = {}
    try:
        for i in range(0, len(unknown), 2):
            key = unknown[i].lstrip("-")
            val = float(unknown[i+1])
            shape_params[key] = val
    except (IndexError, ValueError):
        print(json.dumps({"error": "Dimension arguments must be provided as key-value pairs, e.g. --h 300 --b 150"}))
        sys.exit(1)

    shape_func = SHAPE_FUNCS.get(args.shape)
    if not shape_func:
        print(json.dumps({"error": f"Shape '{args.shape}' not supported for steel design"}))
        sys.exit(1)

    try:
        props = shape_func(shape_params)
    except KeyError as e:
        print(json.dumps({"error": f"Missing required dimension for {args.shape}: {e}"}))
        sys.exit(1)

    if args.calc == "classify":
        result = classify_section(props, args.fy)
    elif args.calc == "moment":
        result = moment_resistance(props, args.fy, args.gamma_M0)
    elif args.calc == "shear":
        result = shear_resistance(props, args.fy, args.gamma_M0)
    elif args.calc == "ltb":
        if args.Lcr <= 0:
            print(json.dumps({"error": "--Lcr (lateral buckling length in mm) is required for LTB check"}))
            sys.exit(1)
        result = ltb_resistance(props, args.fy, args.Lcr, args.C1, args.gamma_M1)

    print(json.dumps(result, indent=2))

if __name__ == "__main__":
    main()

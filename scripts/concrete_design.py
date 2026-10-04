# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "numpy>=1.24",
# ]
# ///
"""
Reinforced Concrete Design per Eurocode 2 (EN 1992-1-1).

Modules:
1. Flexural design — required reinforcement area (As_req)
2. Shear design — required shear reinforcement (Asw/s)
3. Crack width check (wk) per clause 7.3.4

Units: mm, MPa, kN·m for inputs. Outputs clearly labelled.
"""

import argparse
import json
import math
import sys


# ──────────────────── material helpers ──────────────────

def concrete_props(fck: float):
    """Derive concrete design properties from fck (MPa)."""
    gamma_c = 1.5
    alpha_cc = 1.0  # (some national annexes use 0.85)
    fcd = alpha_cc * fck / gamma_c
    fctm = 0.30 * fck**(2 / 3) if fck <= 50 else 2.12 * math.log(1 + (fck + 8) / 10)
    Ecm = 22000 * ((fck + 8) / 10)**0.3  # MPa
    return {"fcd": fcd, "fctm": fctm, "Ecm": Ecm, "gamma_c": gamma_c}


def steel_props(fyk: float):
    gamma_s = 1.15
    fyd = fyk / gamma_s
    Es = 200000  # MPa
    return {"fyd": fyd, "Es": Es, "gamma_s": gamma_s}


# ──────────────────── flexural design ──────────────────

def flexural_design(b: float, d: float, fck: float, fyk: float, M_Ed: float):
    """
    Rectangular section flexural design (singly reinforced).

    Parameters: b (mm), d (mm), fck (MPa), fyk (MPa), M_Ed (kN·m)
    """
    c = concrete_props(fck)
    s = steel_props(fyk)
    fcd = c["fcd"]
    fyd = s["fyd"]

    M_Ed_Nmm = M_Ed * 1e6  # kN·m → N·mm

    # Normalised moment coefficient
    mu = M_Ed_Nmm / (b * d**2 * fcd)

    # Check if compression reinforcement is needed
    # Limiting mu for singly reinforced (ξ_bal ≈ 0.45 for C≤50, x/d limit)
    xi_lim = 0.45  # EN 1992 recommended
    mu_lim = xi_lim * (1 - 0.4 * xi_lim)

    if mu > mu_lim:
        return {
            "calc": "flexure",
            "status": "COMPRESSION_REINFORCEMENT_NEEDED",
            "mu": round(mu, 4),
            "mu_lim": round(mu_lim, 4),
            "note": f"μ = {mu:.4f} > μ_lim = {mu_lim:.4f}. Section needs compression steel or increased depth."
        }

    # Depth of compression zone
    xi = 1 - math.sqrt(1 - 2 * mu)
    x = xi * d  # neutral axis depth
    z = d * (1 - 0.4 * xi)  # lever arm

    # Required tension reinforcement
    As_req = M_Ed_Nmm / (fyd * z)

    # Minimum reinforcement (cl. 9.2.1.1)
    As_min = max(0.26 * c["fctm"] / fyk * b * d, 0.0013 * b * d)

    # Maximum reinforcement
    As_max = 0.04 * b * d  # conservative estimate for total section height ≈ d + cover

    As_design = max(As_req, As_min)

    return {
        "calc": "flexure",
        "code": "EN 1992-1-1",
        "inputs": {
            "b_mm": b, "d_mm": d,
            "fck_MPa": fck, "fyk_MPa": fyk,
            "fcd_MPa": round(fcd, 2), "fyd_MPa": round(fyd, 2),
            "M_Ed_kNm": M_Ed,
        },
        "results": {
            "mu": round(mu, 4),
            "mu_lim": round(mu_lim, 4),
            "xi (x/d)": round(xi, 4),
            "x_mm": round(x, 1),
            "z_mm": round(z, 1),
            "As_req_mm2": round(As_req, 1),
            "As_min_mm2": round(As_min, 1),
            "As_max_mm2": round(As_max, 1),
            "As_design_mm2": round(As_design, 1),
        },
        "status": "OK" if As_design <= As_max else "SECTION_TOO_SMALL",
    }


# ──────────────────── shear design ─────────────────────

def shear_design(b: float, d: float, fck: float, fyk: float,
                  V_Ed: float, As_prov: float):
    """
    Shear design per EN 1992-1-1 clause 6.2.

    V_Ed in kN, As_prov (longitudinal reinforcement) in mm².
    """
    c = concrete_props(fck)
    s = steel_props(fyk)
    fcd = c["fcd"]
    fyd = s["fyd"]
    fywd = min(fyd, 0.8 * fyk)  # design yield of shear reinforcement

    V_Ed_N = V_Ed * 1e3  # kN → N

    # Concrete shear resistance without reinforcement (cl. 6.2.2)
    rho_l = min(As_prov / (b * d), 0.02)
    k = min(1 + math.sqrt(200 / d), 2.0)
    CRd_c = 0.18 / c["gamma_c"]
    k1 = 0.15
    sigma_cp = 0  # no axial force assumed

    V_Rd_c = (CRd_c * k * (100 * rho_l * fck)**(1 / 3) + k1 * sigma_cp) * b * d
    v_min = 0.035 * k**1.5 * fck**0.5
    V_Rd_c_min = (v_min + k1 * sigma_cp) * b * d
    V_Rd_c = max(V_Rd_c, V_Rd_c_min)

    needs_shear_reinf = V_Ed_N > V_Rd_c

    result = {
        "calc": "shear",
        "code": "EN 1992-1-1 cl.6.2",
        "inputs": {
            "b_mm": b, "d_mm": d,
            "fck_MPa": fck, "fyk_MPa": fyk,
            "V_Ed_kN": V_Ed,
            "As_prov_mm2": As_prov,
        },
        "concrete_resistance": {
            "k": round(k, 3),
            "rho_l": round(rho_l, 5),
            "V_Rd_c_kN": round(V_Rd_c / 1e3, 2),
        },
        "needs_shear_reinforcement": needs_shear_reinf,
    }

    if needs_shear_reinf:
        # Variable strut inclination method (cl. 6.2.3)
        # Assume θ = 21.8° (cot θ = 2.5) initially, then check max
        nu = 0.6 * (1 - fck / 250)
        V_Rd_max_cot25 = (nu * fcd * b * 0.9 * d * 2.5) / (1 + 2.5**2)
        V_Rd_max_cot10 = (nu * fcd * b * 0.9 * d * 1.0) / (1 + 1.0**2)

        if V_Ed_N <= V_Rd_max_cot25:
            cot_theta = 2.5
            theta_deg = math.degrees(math.atan(1 / 2.5))
        elif V_Ed_N <= V_Rd_max_cot10:
            # Find theta from V_Ed = V_Rd_max formula
            # V_Ed = nu * fcd * b * z * cot_theta / (1 + cot_theta^2)
            # Solve numerically
            z = 0.9 * d
            for ct in [x / 10 for x in range(25, 9, -1)]:
                v_rd = (nu * fcd * b * z * ct) / (1 + ct**2)
                if v_rd >= V_Ed_N:
                    cot_theta = ct
                    break
            else:
                cot_theta = 1.0
            theta_deg = math.degrees(math.atan(1 / cot_theta))
        else:
            result["status"] = "SECTION_TOO_SMALL_FOR_SHEAR"
            result["V_Rd_max_kN"] = round(V_Rd_max_cot10 / 1e3, 2)
            return result

        z = 0.9 * d
        Asw_over_s = V_Ed_N / (z * fywd * cot_theta)  # mm²/mm

        # Minimum shear reinforcement (cl. 9.2.2)
        rho_w_min = 0.08 * math.sqrt(fck) / fyk
        Asw_s_min = rho_w_min * b  # mm²/mm

        Asw_s_design = max(Asw_over_s, Asw_s_min)

        result["strut_and_tie"] = {
            "theta_deg": round(theta_deg, 1),
            "cot_theta": round(cot_theta, 2),
            "z_mm": round(z, 1),
            "V_Rd_max_kN": round((nu * fcd * b * z * cot_theta / (1 + cot_theta**2)) / 1e3, 2),
            "Asw_s_req_mm2_per_mm": round(Asw_over_s, 4),
            "Asw_s_min_mm2_per_mm": round(Asw_s_min, 4),
            "Asw_s_design_mm2_per_mm": round(Asw_s_design, 4),
        }
        result["status"] = "SHEAR_REINFORCEMENT_DESIGNED"
    else:
        result["status"] = "CONCRETE_ALONE_SUFFICIENT"

    return result


# ──────────────────── crack width check ────────────────

def crack_width_check(b: float, d: float, h: float, fck: float, fyk: float,
                       As_prov: float, M_Ed: float,
                       bar_dia: float, cover: float, spacing: float):
    """
    Crack width calculation per EN 1992-1-1 clause 7.3.4.

    M_Ed in kN·m, dimensions in mm.
    """
    c = concrete_props(fck)
    s = steel_props(fyk)
    Es = s["Es"]
    Ecm = c["Ecm"]
    fctm = c["fctm"]

    M_Ed_Nmm = M_Ed * 1e6

    alpha_e = Es / Ecm  # modular ratio

    # Cracked section neutral axis depth (transformed section)
    rho = As_prov / (b * d)
    xi_cr = alpha_e * rho * (-1 + math.sqrt(1 + 2 / (alpha_e * rho)))
    x_cr = xi_cr * d

    # Cracked moment of inertia
    I_cr = b * x_cr**3 / 3 + alpha_e * As_prov * (d - x_cr)**2

    # Steel stress under service load
    sigma_s = alpha_e * M_Ed_Nmm * (d - x_cr) / I_cr

    # Effective tension area (cl. 7.3.2)
    h_c_eff = min(2.5 * (h - d), (h - x_cr) / 3, h / 2)
    A_c_eff = b * h_c_eff
    rho_p_eff = As_prov / A_c_eff

    # Maximum crack spacing (cl. 7.3.4 eq 7.11)
    k1 = 0.8  # high bond bars
    k2 = 0.5  # bending
    k3 = 3.4  # recommended
    k4 = 0.425  # recommended

    s_r_max = k3 * cover + k1 * k2 * k4 * bar_dia / rho_p_eff

    # Check if spacing > 5(c + φ/2) → use simplified formula
    if spacing > 5 * (cover + bar_dia / 2):
        s_r_max = 1.3 * (h - x_cr)

    # Mean strain difference (eq 7.9)
    kt = 0.4  # long-term loading
    eps_sm_minus_eps_cm = (sigma_s - kt * fctm / rho_p_eff * (1 + alpha_e * rho_p_eff)) / Es
    eps_sm_minus_eps_cm = max(eps_sm_minus_eps_cm, 0.6 * sigma_s / Es)

    # Crack width
    wk = s_r_max * eps_sm_minus_eps_cm

    # Limits per exposure class (typical)
    wk_limit = 0.3  # mm — common for XC2, XC3, XC4

    return {
        "calc": "crack_width",
        "code": "EN 1992-1-1 cl.7.3.4",
        "inputs": {
            "b_mm": b, "d_mm": d, "h_mm": h,
            "fck_MPa": fck, "fyk_MPa": fyk,
            "As_prov_mm2": As_prov,
            "M_Ed_kNm": M_Ed,
            "bar_dia_mm": bar_dia,
            "cover_mm": cover,
            "spacing_mm": spacing,
        },
        "results": {
            "alpha_e": round(alpha_e, 2),
            "x_cr_mm": round(x_cr, 1),
            "I_cr_mm4": round(I_cr, 0),
            "sigma_s_MPa": round(sigma_s, 1),
            "h_c_eff_mm": round(h_c_eff, 1),
            "rho_p_eff": round(rho_p_eff, 5),
            "s_r_max_mm": round(s_r_max, 1),
            "eps_sm_minus_eps_cm": round(eps_sm_minus_eps_cm, 6),
            "wk_mm": round(wk, 3),
            "wk_limit_mm": wk_limit,
        },
        "status": "OK" if wk <= wk_limit else "CRACK_WIDTH_EXCEEDED",
        "utilisation_ratio": round(wk / wk_limit, 3),
    }


# ──────────────────── main ─────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="RC design — EN 1992-1-1")
    parser.add_argument("--calc", required=True, choices=["flexure", "shear", "crack_width"])
    parser.add_argument("--b", type=float, help="Section width (mm)")
    parser.add_argument("--d", type=float, help="Effective depth (mm)")
    parser.add_argument("--h", type=float, help="Total section height (mm)")
    parser.add_argument("--fck", type=float, help="Concrete strength fck (MPa)")
    parser.add_argument("--fyk", type=float, default=500, help="Steel yield fyk (MPa)")
    parser.add_argument("--M_Ed", type=float, default=0, help="Design moment (kN·m)")
    parser.add_argument("--V_Ed", type=float, default=0, help="Design shear force (kN)")
    parser.add_argument("--As_prov", type=float, default=0, help="Provided reinforcement (mm²)")
    parser.add_argument("--bar_dia", type=float, default=16, help="Bar diameter (mm)")
    parser.add_argument("--cover", type=float, default=35, help="Concrete cover (mm)")
    parser.add_argument("--spacing", type=float, default=150, help="Bar spacing (mm)")

    args = parser.parse_args()

    if args.calc == "flexure":
        if not all([args.b, args.d, args.fck, args.M_Ed]):
            print(json.dumps({"error": "Flexure requires: --b, --d, --fck, --M_Ed"}))
            sys.exit(1)
        result = flexural_design(args.b, args.d, args.fck, args.fyk, args.M_Ed)

    elif args.calc == "shear":
        if not all([args.b, args.d, args.fck, args.V_Ed, args.As_prov]):
            print(json.dumps({"error": "Shear requires: --b, --d, --fck, --V_Ed, --As_prov"}))
            sys.exit(1)
        result = shear_design(args.b, args.d, args.fck, args.fyk, args.V_Ed, args.As_prov)

    elif args.calc == "crack_width":
        if not all([args.b, args.d, args.h, args.fck, args.As_prov, args.M_Ed]):
            print(json.dumps({"error": "Crack width requires: --b, --d, --h, --fck, --As_prov, --M_Ed"}))
            sys.exit(1)
        result = crack_width_check(
            args.b, args.d, args.h, args.fck, args.fyk,
            args.As_prov, args.M_Ed, args.bar_dia, args.cover, args.spacing,
        )

    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()

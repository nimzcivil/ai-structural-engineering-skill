# /// script
# requires-python = ">=3.10"
# ///
"""
Interaction Checks per Simplified Methods (Eurocode 3 / EAE).
Adds knowledge extracted from the provided structural mechanics PDF.

Modules:
1. shear_bending: Reduced moment capacity when V_Ed > 0.5 V_pl,Rd
2. buckling_interaction: Simplified combined bending + compression check
"""

import argparse
import json
import sys

def shear_bending(params):
    """
    Calculates reduced moment resistance M_V,Rd if shear is high.
    Expected params: V_Ed (kN), V_pl_Rd (kN), M_c_Rd (kNm), A_w (mm2), t_w (mm), fy (MPa)
    """
    V_Ed = params["V_Ed"]
    V_pl_Rd = params["V_pl_Rd"]
    M_c_Rd = params["M_c_Rd"]
    A_w = params["A_w"]
    t_w = params["t_w"]
    fy = params["fy"]
    gamma_M0 = params.get("gamma_M0", 1.0)

    if V_Ed <= 0.5 * V_pl_Rd:
        return {
            "calc": "shear_bending_interaction",
            "status": "NO_REDUCTION",
            "condition": f"V_Ed ({V_Ed} kN) <= 0.5 * V_pl_Rd ({0.5*V_pl_Rd} kN)",
            "M_V_Rd_kNm": M_c_Rd
        }
    
    # Plastic reduction factor rho
    rho = (2 * V_Ed / V_pl_Rd - 1)**2
    
    # Reduction term in N*mm, converted to kNm
    reduction_kNm = (rho * (A_w**2) / (4 * t_w)) * fy / 1e6 / gamma_M0
    M_V_Rd = M_c_Rd - reduction_kNm

    return {
        "calc": "shear_bending_interaction",
        "status": "REDUCED_MOMENT_CAPACITY",
        "condition": f"V_Ed ({V_Ed} kN) > 0.5 * V_pl_Rd ({0.5*V_pl_Rd} kN)",
        "rho": round(rho, 4),
        "reduction_kNm": round(reduction_kNm, 2),
        "M_c_Rd_original_kNm": M_c_Rd,
        "M_V_Rd_kNm": round(M_V_Rd, 2)
    }

def buckling_interaction(params):
    """
    Simplified Method for members subjected to Bending and Compression.
    Expected params:
    N_Ed, N_b_Rd, N_cr_y, N_cr_z (in kN)
    M_y_Ed, M_z_Ed, M_y_Rk, M_z_Rk (in kNm)
    chi_LT, c_my, c_mz, gamma_M1
    """
    N_Ed = params.get("N_Ed", 0)
    N_b_Rd = params["N_b_Rd"]
    N_cr_y = params.get("N_cr_y", 1e9)  # Default to infinity if not provided (e.g., no buckling risk)
    N_cr_z = params.get("N_cr_z", 1e9)
    
    M_y_Ed = params.get("M_y_Ed", 0)
    M_z_Ed = params.get("M_z_Ed", 0)
    M_y_Rk = params.get("M_y_Rk", 1e9)
    M_z_Rk = params.get("M_z_Rk", 1e9)
    
    chi_LT = params.get("chi_LT", 1.0)
    c_my = params.get("c_my", 1.0)
    c_mz = params.get("c_mz", 1.0)
    gamma_M1 = params.get("gamma_M1", 1.0)

    # 1. Axial term
    term_N = N_Ed / N_b_Rd

    # 2. Bending y-y term (with amplification factor for N_Ed / N_cr_y)
    amp_y = 1.0
    if N_Ed > 0 and N_cr_y < 1e9:
        amp_y = 1 / (1 - N_Ed / N_cr_y)
        if amp_y < 1.0: amp_y = 1.0 # Cannot be less than 1
        
    term_My = amp_y * (c_my * M_y_Ed) / (chi_LT * M_y_Rk / gamma_M1) if M_y_Rk < 1e9 else 0.0

    # 3. Bending z-z term (with amplification factor for N_Ed / N_cr_z)
    amp_z = 1.0
    if N_Ed > 0 and N_cr_z < 1e9:
        amp_z = 1 / (1 - N_Ed / N_cr_z)
        if amp_z < 1.0: amp_z = 1.0

    term_Mz = amp_z * (c_mz * M_z_Ed) / (M_z_Rk / gamma_M1) if M_z_Rk < 1e9 else 0.0

    utilization = term_N + term_My + term_Mz

    return {
        "calc": "simplified_buckling_interaction",
        "terms": {
            "axial_utilization": round(term_N, 4),
            "bending_y_amplified_utilization": round(term_My, 4),
            "bending_z_amplified_utilization": round(term_Mz, 4)
        },
        "factors": {
            "amplification_y_axis": round(amp_y, 4),
            "amplification_z_axis": round(amp_z, 4)
        },
        "total_utilization": round(utilization, 4),
        "status": "PASS" if utilization <= 1.0 else "FAIL"
    }

def main():
    parser = argparse.ArgumentParser(description="Structural Interaction Checks")
    parser.add_argument("--calc", choices=["shear_bending", "buckling_interaction"], required=True)
    
    args, unknown = parser.parse_known_args()
    
    params = {}
    try:
        for i in range(0, len(unknown), 2):
            key = unknown[i].lstrip("-")
            val = float(unknown[i+1])
            params[key] = val
    except (IndexError, ValueError):
        print(json.dumps({"error": "Interaction arguments must be key-value pairs, e.g. --V_Ed 17.5"}))
        sys.exit(1)

    try:
        if args.calc == "shear_bending":
            result = shear_bending(params)
        elif args.calc == "buckling_interaction":
            result = buckling_interaction(params)
        print(json.dumps(result, indent=2))
    except KeyError as e:
        print(json.dumps({"error": f"Missing required parameter: {e}"}))
        sys.exit(1)

if __name__ == "__main__":
    main()

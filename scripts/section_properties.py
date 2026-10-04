# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "numpy>=1.24",
# ]
# ///
"""
Cross-Section Properties Calculator.

Computes: Area, centroid, moments of inertia (Iy, Iz), elastic section modulus
(Wel), plastic section modulus (Wpl), radius of gyration for common structural
shapes.

All inputs in mm, outputs in mm-based units (mm², mm⁴, mm³).
"""

import argparse
import json
import math
import sys

import numpy as np


def rectangle(params):
    b = params["b"]
    h = params["h"]
    A = b * h
    Iy = b * h**3 / 12
    Iz = h * b**3 / 12
    Wel_y = b * h**2 / 6
    Wel_z = h * b**2 / 6
    Wpl_y = b * h**2 / 4
    Wpl_z = h * b**2 / 4
    iy = math.sqrt(Iy / A)
    iz = math.sqrt(Iz / A)
    return {
        "shape": "rectangle",
        "dimensions_mm": params,
        "A_mm2": round(A, 2),
        "Iy_mm4": round(Iy, 2),
        "Iz_mm4": round(Iz, 2),
        "Wel_y_mm3": round(Wel_y, 2),
        "Wel_z_mm3": round(Wel_z, 2),
        "Wpl_y_mm3": round(Wpl_y, 2),
        "Wpl_z_mm3": round(Wpl_z, 2),
        "iy_mm": round(iy, 2),
        "iz_mm": round(iz, 2),
        "centroid_y_mm": round(h / 2, 2),
        "centroid_z_mm": round(b / 2, 2),
    }


def circle(params):
    if "d" in params:
        d = params["d"]
        r = d / 2
    else:
        r = params["r"]
        d = 2 * r
    A = math.pi * r**2
    I = math.pi * d**4 / 64
    Wel = math.pi * d**3 / 32
    Wpl = d**3 / 6
    i_r = d / 4
    return {
        "shape": "circle",
        "diameter_mm": round(d, 2),
        "A_mm2": round(A, 2),
        "I_mm4": round(I, 2),
        "Wel_mm3": round(Wel, 2),
        "Wpl_mm3": round(Wpl, 2),
        "i_mm": round(i_r, 2),
    }


def hollow_circle(params):
    d_out = params["d_out"]
    d_in = params["d_in"]
    if d_in >= d_out:
        return {"error": "Inner diameter must be less than outer diameter"}
    A = math.pi / 4 * (d_out**2 - d_in**2)
    I = math.pi / 64 * (d_out**4 - d_in**4)
    Wel = I / (d_out / 2)
    Wpl = (d_out**3 - d_in**3) / 6
    i_r = math.sqrt(I / A)
    return {
        "shape": "hollow_circle",
        "d_out_mm": d_out,
        "d_in_mm": d_in,
        "A_mm2": round(A, 2),
        "I_mm4": round(I, 2),
        "Wel_mm3": round(Wel, 2),
        "Wpl_mm3": round(Wpl, 2),
        "i_mm": round(i_r, 2),
    }


def I_section(params):
    """Symmetric I/H section."""
    h = params["h"]       # total height
    b = params["b"]       # flange width
    tf = params["tf"]     # flange thickness
    tw = params["tw"]     # web thickness
    hw = h - 2 * tf       # web clear height

    A = 2 * b * tf + hw * tw
    Iy = (b * h**3 - (b - tw) * hw**3) / 12
    Iz = (2 * tf * b**3 + hw * tw**3) / 12
    Wel_y = Iy / (h / 2)
    Wel_z = Iz / (b / 2)

    # Plastic section modulus about strong axis
    # Two flanges + web
    Wpl_y = 2 * b * tf * (h / 2 - tf / 2) + tw * hw**2 / 4
    Wpl_z = 2 * (tf * b**2 / 4) + hw * tw**2 / 4

    iy = math.sqrt(Iy / A)
    iz = math.sqrt(Iz / A)

    # Warping constant and torsional constant (approx)
    Iw = (tf * b**3 * (h - tf)**2) / 24  # warping constant (mm⁶)
    It = (2 * b * tf**3 + hw * tw**3) / 3  # St-Venant torsion (mm⁴)

    return {
        "shape": "I_section",
        "dimensions_mm": params,
        "hw_mm": round(hw, 2),
        "A_mm2": round(A, 2),
        "Iy_mm4": round(Iy, 2),
        "Iz_mm4": round(Iz, 2),
        "Wel_y_mm3": round(Wel_y, 2),
        "Wel_z_mm3": round(Wel_z, 2),
        "Wpl_y_mm3": round(Wpl_y, 2),
        "Wpl_z_mm3": round(Wpl_z, 2),
        "iy_mm": round(iy, 2),
        "iz_mm": round(iz, 2),
        "It_mm4": round(It, 2),
        "Iw_mm6": round(Iw, 2),
        "centroid_y_mm": round(h / 2, 2),
    }


def T_section(params):
    """T-section (flange on top)."""
    b = params["b"]
    h = params["h"]
    tf = params["tf"]
    tw = params["tw"]
    hw = h - tf

    A_f = b * tf
    A_w = hw * tw
    A = A_f + A_w

    # centroid from bottom
    y_bar = (A_w * hw / 2 + A_f * (hw + tf / 2)) / A

    Iy = (tw * hw**3 / 12 + A_w * (y_bar - hw / 2)**2 +
          b * tf**3 / 12 + A_f * (hw + tf / 2 - y_bar)**2)
    Iz = (hw * tw**3 / 12 + tf * b**3 / 12)

    Wel_y_top = Iy / (h - y_bar)
    Wel_y_bot = Iy / y_bar

    iy = math.sqrt(Iy / A)
    iz = math.sqrt(Iz / A)

    return {
        "shape": "T_section",
        "dimensions_mm": params,
        "A_mm2": round(A, 2),
        "centroid_from_bottom_mm": round(y_bar, 2),
        "Iy_mm4": round(Iy, 2),
        "Iz_mm4": round(Iz, 2),
        "Wel_y_top_mm3": round(Wel_y_top, 2),
        "Wel_y_bot_mm3": round(Wel_y_bot, 2),
        "iy_mm": round(iy, 2),
        "iz_mm": round(iz, 2),
    }


def channel(params):
    """C-channel (symmetric about y-axis)."""
    h = params["h"]
    b = params["b"]
    tf = params["tf"]
    tw = params["tw"]
    hw = h - 2 * tf

    A = 2 * b * tf + hw * tw

    # centroid from web face
    z_bar = (2 * b * tf * b / 2 + hw * tw * tw / 2) / A

    Iy = (tw * hw**3 / 12 + 2 * (b * tf**3 / 12 + b * tf * ((h - tf) / 2)**2))
    Iz = (2 * (tf * b**3 / 12 + b * tf * (b / 2 - z_bar)**2) +
          hw * tw**3 / 12 + hw * tw * (z_bar - tw / 2)**2)

    Wel_y = Iy / (h / 2)
    iy = math.sqrt(Iy / A)
    iz = math.sqrt(Iz / A)

    return {
        "shape": "channel",
        "dimensions_mm": params,
        "A_mm2": round(A, 2),
        "centroid_from_web_mm": round(z_bar, 2),
        "Iy_mm4": round(Iy, 2),
        "Iz_mm4": round(Iz, 2),
        "Wel_y_mm3": round(Wel_y, 2),
        "iy_mm": round(iy, 2),
        "iz_mm": round(iz, 2),
    }


def hollow_rect(params):
    """Rectangular hollow section (RHS/SHS)."""
    b = params["b"]
    h = params["h"]
    t = params["t"]
    b_i = b - 2 * t
    h_i = h - 2 * t

    if b_i <= 0 or h_i <= 0:
        return {"error": "Wall thickness too large for given dimensions"}

    A = b * h - b_i * h_i
    Iy = (b * h**3 - b_i * h_i**3) / 12
    Iz = (h * b**3 - h_i * b_i**3) / 12
    Wel_y = Iy / (h / 2)
    Wel_z = Iz / (b / 2)
    iy = math.sqrt(Iy / A)
    iz = math.sqrt(Iz / A)

    return {
        "shape": "hollow_rect",
        "dimensions_mm": params,
        "A_mm2": round(A, 2),
        "Iy_mm4": round(Iy, 2),
        "Iz_mm4": round(Iz, 2),
        "Wel_y_mm3": round(Wel_y, 2),
        "Wel_z_mm3": round(Wel_z, 2),
        "iy_mm": round(iy, 2),
        "iz_mm": round(iz, 2),
    }


SHAPES = {
    "rectangle": rectangle,
    "circle": circle,
    "hollow_circle": hollow_circle,
    "I_section": I_section,
    "T_section": T_section,
    "channel": channel,
    "hollow_rect": hollow_rect,
}


def main():
    parser = argparse.ArgumentParser(description="Cross-section properties calculator")
    parser.add_argument("--shape", required=True, choices=SHAPES.keys())
    
    # Parse known args, everything else is treated as dimension parameters (e.g. --h 300)
    args, unknown = parser.parse_known_args()

    params = {}
    try:
        # Parse unknown args like "--h", "300", "--b", "150"
        for i in range(0, len(unknown), 2):
            key = unknown[i].lstrip("-")
            val = float(unknown[i+1])
            params[key] = val
    except (IndexError, ValueError):
        print(json.dumps({"error": "Dimension arguments must be provided as key-value pairs, e.g. --h 300 --b 150"}))
        sys.exit(1)

    func = SHAPES[args.shape]
    try:
        result = func(params)
        print(json.dumps(result, indent=2))
    except KeyError as e:
        print(json.dumps({"error": f"Missing required dimension for {args.shape}: {e}"}))
        sys.exit(1)

if __name__ == "__main__":
    main()

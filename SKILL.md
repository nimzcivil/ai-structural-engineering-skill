---
name: structural-engineering
description: >-
  Performs structural engineering calculations using Python computation scripts.
  Use when the user asks to analyze beams (reactions, shear, bending moment,
  deflection), design columns (Euler buckling, slenderness), compute
  cross-section properties (area, moment of inertia, section modulus),
  design reinforced concrete elements per Eurocode 2 (EN 1992), or check
  steel members per Eurocode 3 (EN 1993). Never perform these calculations
  yourself — always run the appropriate Python script.
---

# Structural Engineering Calculation Skill

> **CRITICAL**: You must NEVER perform structural calculations yourself.
> Always run the appropriate Python script from this skill's `scripts/`
> directory using `uv run`. The scripts use `numpy`, `sympy`, and `scipy`
> for verified numerical computation.

## Prerequisites

Before running any script, ensure `uv` is installed by following the
[uv skill](../../../.gemini/config/plugins/science/skills/uv/SKILL.md).

## Available Modules

| Module | Script | What it computes |
|--------|--------|------------------|
| Beam Analysis | `scripts/beam_analysis.py` | Reactions, shear diagram values, bending moment diagram values, deflection for simply-supported, cantilever, fixed-fixed, and propped-cantilever beams |
| Cross-Section Properties | `scripts/section_properties.py` | Area, centroid, moment of inertia (Iy, Iz), section modulus (Wel, Wpl), radius of gyration for rectangle, circle, I/H, T, channel, hollow sections |
| Column Design (EC3) | `scripts/column_design.py` | Euler critical load, slenderness ratio, buckling resistance per EN 1993-1-1, buckling curve selection |
| Concrete Design (EC2) | `scripts/concrete_design.py` | Flexural design (As_req), shear design (Asw), crack width check per EN 1992-1-1 |
| Steel Design (EC3) | `scripts/steel_design.py` | Cross-section classification, moment resistance (Mc,Rd), shear resistance (Vpl,Rd), lateral-torsional buckling (Mb,Rd) per EN 1993-1-1 |
| Combined Interaction | `scripts/interaction.py` | Shear-bending moment capacity reduction, Simplified bending-axial buckling interaction |

## How to Use

### 1. Beam Analysis

```powershell
uv run scripts/beam_analysis.py --type simply_supported --length 6.0 --E 210e9 --I 8.356e-5 --loads '[{\"type\":\"point\",\"P\":50000,\"a\":3.0}]' --distributed '[{\"w\":10000,\"a\":0,\"b\":6.0}]'
```

**Beam types**: `simply_supported`, `cantilever`, `fixed_fixed`, `propped_cantilever`

**Load format (--loads)**: JSON array of point loads:
- `{"type": "point", "P": <force_N>, "a": <distance_from_left_m>}`
- `{"type": "moment", "M": <moment_Nm>, "a": <distance_from_left_m>}`

**Distributed load format (--distributed)**: JSON array:
- `{"w": <load_N_per_m>, "a": <start_m>, "b": <end_m>}`

> **PowerShell quoting**: Use `'{\"key\": value}'` (backslash-escaped
> double quotes inside single quotes) for JSON arguments on Windows.

### 2. Cross-Section Properties

```powershell
uv run scripts/section_properties.py --shape I_section --h 300 --b 150 --tf 10.7 --tw 7.1
```

**Shapes & required params** (all dimensions in mm, pass as direct flags e.g., `--h 300 --b 150`):
- `rectangle`: `--b` (width), `--h` (height)
- `circle`: `--d` (diameter) or `--r` (radius)
- `hollow_circle`: `--d_out`, `--d_in`
- `I_section`: `--h`, `--b`, `--tf`, `--tw`
- `T_section`: `--h`, `--b`, `--tf`, `--tw`
- `channel`: `--h`, `--b`, `--tf`, `--tw`
- `hollow_rect`: `--b` (outer width), `--h` (outer height), `--t` (wall thickness)

### 3. Column Design (Eurocode 3)

```powershell
uv run scripts/column_design.py --length 4.0 --A 5380 --I_min 6038000 --fy 355 --E 210000 --buckling_curve b --k 1.0
```

**Parameters**:
- `--length`: Column length (m)
- `--A`: Cross-section area (**mm²**)
- `--I_min`: Minimum moment of inertia (**mm⁴**)
- `--fy`: Yield strength (**MPa**)
- `--E`: Young's modulus (**MPa**) [default: 210000]
- `--buckling_curve`: `a0`, `a`, `b`, `c`, `d` per EN 1993-1-1 Table 6.2
- `--k`: Effective length factor [default: 1.0]

### 4. Concrete Design (Eurocode 2)

```powershell
# Flexural design
uv run scripts/concrete_design.py --calc flexure --b 300 --d 450 --fck 30 --fyk 500 --M_Ed 250

# Shear design
uv run scripts/concrete_design.py --calc shear --b 300 --d 450 --fck 30 --fyk 500 --V_Ed 180 --As_prov 1256

# Crack width
uv run scripts/concrete_design.py --calc crack_width --b 300 --d 450 --h 500 --fck 30 --fyk 500 --As_prov 1256 --M_Ed 150 --bar_dia 16 --cover 35 --spacing 125
```

**Parameters** (flexure/shear/crack):
- `--b`: Section width (mm)
- `--d`: Effective depth (mm)
- `--h`: Total section height (mm) — needed for crack width
- `--fck`: Characteristic concrete strength (MPa)
- `--fyk`: Characteristic steel yield strength (MPa)
- `--M_Ed`: Design bending moment (kN·m)
- `--V_Ed`: Design shear force (kN)
- `--As_prov`: Provided reinforcement area (mm²)

### 5. Steel Design (Eurocode 3)

```powershell
# Moment resistance
uv run scripts/steel_design.py --calc moment --shape I_section --h 300 --b 150 --tf 10.7 --tw 7.1 --fy 355 --gamma_M0 1.0

# Lateral-torsional buckling
uv run scripts/steel_design.py --calc ltb --shape I_section --h 300 --b 150 --tf 10.7 --tw 7.1 --fy 355 --Lcr 6000 --C1 1.0
```

**Calculation types**: `classify`, `moment`, `shear`, `ltb`

### 6. Combined Interaction (Eurocode 3)

Calculates combined interaction effects (Shear+Bending or Bending+Compression buckling) using the simplified methods from EAE / Eurocode 3.

```powershell
# Shear-Bending Interaction
uv run scripts/interaction.py --calc shear_bending --V_Ed 120 --V_pl_Rd 186.7 --M_c_Rd 97.35 --A_w 1176 --t_w 8 --fy 275

# Bending-Compression Buckling Interaction (Simplified Method)
uv run scripts/interaction.py --calc buckling_interaction --N_Ed 100 --N_b_Rd 576 --N_cr_y 7816 --N_cr_z 1150 --M_y_Ed 70 --M_y_Rk 97.35 --chi_LT 0.833 --c_my 1.0
```

## Output

All scripts output structured JSON to stdout with clearly labeled results,
units, and pass/fail status where applicable. Parse and present the results
to the user in a clear tabular format.

## Error Handling

If a script exits with a non-zero code, display the error message to the user
and ask them to check their input parameters. Common issues:
- Negative or zero dimensions
- Incompatible units (mixing mm and m)
- Missing required parameters for the chosen shape/calc type

## References

Detailed Eurocode clause references are documented in
[references/eurocode_clauses.md](./references/eurocode_clauses.md).

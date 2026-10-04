# Eurocode Clause References

This document lists the specific Eurocode clauses implemented in each script
for traceability and verification purposes.

## EN 1992-1-1 (Concrete — `concrete_design.py`)

| Calculation | Clause | Description |
|---|---|---|
| Material properties | 3.1.2, 3.1.3 | fcd, fctm, Ecm derivation |
| Flexural design | 6.1 | ULS bending, rectangular stress block (λ=0.8, η=1.0 for fck≤50) |
| Neutral axis limit | 5.5 | Redistribution limit, ξ_lim = 0.45 for C≤50 |
| Minimum reinforcement | 9.2.1.1 | As,min = max(0.26·fctm/fyk·bt·d, 0.0013·bt·d) |
| Shear — no reinforcement | 6.2.2 | VRd,c with CRd,c = 0.18/γc, k = min(1+√(200/d), 2) |
| Shear — with reinforcement | 6.2.3 | Variable strut inclination, θ ∈ [21.8°, 45°] |
| Minimum shear reinforcement | 9.2.2 | ρw,min = 0.08·√fck/fyk |
| Crack width | 7.3.4 | sr,max (eq 7.11), εsm−εcm (eq 7.9), wk = sr,max·Δε |
| Effective tension area | 7.3.2 | hc,eff = min(2.5(h−d), (h−x)/3, h/2) |

## EN 1993-1-1 (Steel — `steel_design.py`, `column_design.py`)

| Calculation | Clause | Description |
|---|---|---|
| Section classification | 5.5, Table 5.2 | ε = √(235/fy), c/t limits for Class 1–4 |
| Moment resistance | 6.2.5 | Mc,Rd = Wpl·fy/γM0 (Class 1,2), Wel·fy/γM0 (Class 3) |
| Shear resistance | 6.2.6 | Vpl,Rd = Av·(fy/√3)/γM0 |
| Flexural buckling | 6.3.1, eq 6.49 | Φ = 0.5·[1+α(λ̄−0.2)+λ̄²], χ = 1/(Φ+√(Φ²−λ̄²)) |
| Imperfection factors | Table 6.1 | α for curves a0, a, b, c, d |
| LTB — general method | 6.3.2.2 | Mcr elastic critical moment, χLT reduction |
| LTB curves | Table 6.3, 6.4 | Rolled I: curve a (h/b≤2), curve b (h/b>2) |
| Buckling curves | Table 6.2 | Selection based on section type, axis, fy |

## Beam Analysis (`beam_analysis.py`)

The beam analysis script uses classical structural mechanics (no code-specific
design). Calculations are based on:
- Superposition principle for linear elastic analysis
- Euler-Bernoulli beam theory
- Macaulay's method / SymPy symbolic integration for statically indeterminate beams
- Boundary conditions: δ=0, θ=0 at fixed supports; δ=0 at simple supports

## Section Properties (`section_properties.py`)

Standard geometric formulas for:
- Parallel axis theorem (composite sections)
- St-Venant torsional constant (thin-walled approximation)
- Warping constant for doubly-symmetric I-sections

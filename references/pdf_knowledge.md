# Structural Mechanics Reference Knowledge
*(Extracted from provided Structural Steel Design PDF)*

## 1. Limit State Displacements (Deflections)
Maximum allowable vertical deflections based on load combinations and comfort/visual criteria:
*   **Active deflection ($w_{act}$)**: Overall sum of long-term component deflection under permanent loads ($w_2$) and variable loads ($w_3$).
*   **Floors with thin partition walls or rigid pavements (without sealing)**: $L / 500$
*   **Floors with ordinary partition walls or rigid pavements (with sealing)**: $L / 400$
*   **Other cases**: $L / 300$
*   **Comfort criteria (short-term relative deflection)**: $L / 350$
*   **Visual appearance (quasi-permanent loads)**: $L / 300$

## 2. Combined Bending-Axial Buckling (Simplified Method)
For members subjected to transverse loads or sway frames, $c_M = 1.0$.

The combined verification formula for buckling (Class 1, 2, 3) is given by:
$$ \frac{N_{Ed}}{N_{b,Rd}} + \frac{c_{my} M_{y,Ed}}{\left(1 - \frac{N_{Ed}}{N_{cr,y}}\right) \chi_{LT} \frac{M_{y,Rk}}{\gamma_{M1}}} + \frac{c_{mz} M_{z,Ed}}{\left(1 - \frac{N_{Ed}}{N_{cr,z}}\right) \frac{M_{z,Rk}}{\gamma_{M1}}} \le 1.0 $$

## 3. Shear-Bending Interaction (Plastic Criteria)
*   If $V_{Ed} \le 0.5 V_{pl,Rd}$: No reduction in design resistance.
*   If $V_{Ed} > 0.5 V_{pl,Rd}$: A reduced yield strength $(1-\rho)f_y$ is assigned to the shear area.
*   $\rho = \left( \frac{2 V_{Ed}}{V_{pl,Rd}} - 1 \right)^2$
*   For I-sections with equal flanges (major axis bending), reduced moment resistance is:
    $$ M_{V,y,Rd} = \left[ W_{pl,y} - \frac{\rho A_w^2}{4 t_w} \right] \frac{f_y}{\gamma_{M0}} $$

*(Note: The `interaction.py` script automatically evaluates these interaction constraints.)*

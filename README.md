# AI Structural Engineering Skill 🏗️🤖

An open-source AI agent skill designed to bridge the gap between Large Language Models (LLMs) and deterministic structural engineering calculations (Eurocode). 

## The Problem
LLMs are fantastic at reading codes, understanding context, and planning engineering workflows, but they are notoriously terrible at arithmetic. You cannot trust an LLM to calculate the combined bending and axial buckling resistance of a steel column.

## The Solution
This skill separates **reasoning** from **computation**. It equips AI agents with a set of robust, deterministic Python scripts. The AI acts as the engineer (interpreting the problem and deciding which checks are needed) and delegates the actual math to the Python scripts.

## What's Included
* **Beam Analysis:** Calculates reactions, shear, bending moments, and deflections for various boundary conditions.
* **Section Properties:** Computes Area, Moment of Inertia, Section Modulus, and Radius of Gyration for common shapes.
* **Steel Design (EN 1993-1-1):** Cross-section classification, moment resistance, shear resistance, and lateral-torsional buckling.
* **Concrete Design (EN 1992-1-1):** Flexural design, shear design, and crack width verification.
* **Column Design (EN 1993-1-1):** Flexural buckling resistances and slenderness checks.
* **Combined Interaction:** Shear-bending reductions and simplified V-M-N interaction formulas (per EAE/EC3).

## How it Works
The AI reads the `SKILL.md` file to understand how to use the tools. It then executes the Python scripts (via `uv run`) directly in its environment and parses the JSON output to inform you of the structural utilization and pass/fail states.

## Disclaimer ⚠️
**This is an experimental tool.** It is designed to assist and automate calculations but must **never** be used for actual construction or final design without independent verification by a licensed Professional Engineer (PE) or Structural Engineer (SE).

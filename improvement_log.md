# Skill Improvement Log

## Identified Mistakes & Friction Points
1. **JSON Parsing on Windows (PowerShell)**
   - *Error*: Passing `--params '{"h":160}'` failed because PowerShell strips quotes before they reach Python, causing JSON decode errors. We had to use ugly escapes like `'{\"h\":160}'`.
   - *Fix*: Remove the `--params` JSON argument entirely. Scripts should accept flat command-line arguments for dimensions (e.g., `--h 160 --b 160 --tf 13 --tw 8`).

2. **Error-Prone Unit Inputs (Column Design)**
   - *Error*: `column_design.py` expected Area in $m^2$ and Inertia in $m^4$. I manually miscalculated the scientific notation for $m^4$, which could lead to disastrously wrong buckling capacities if unchecked.
   - *Fix*: Standardize ALL scripts to accept structural engineering units: **mm, mm², mm⁴, MPa, kN, kNm**. The Python scripts must handle the conversion to SI base units internally.

3. **Complex List Parsing (Beam Analysis)**
   - *Error*: The `beam_analysis.py` also relies on JSON arrays for loads, which suffers from the same PowerShell escaping issues. 
   - *Fix*: While JSON is hard to avoid for lists of objects (like multiple point loads), we will document the exact PowerShell/Bash escaping rules clearly, or simplify the parser. For now, we will update the documentation to emphasize the quotes.

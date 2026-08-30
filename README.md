# Repair Identifiability and Product-Measurement Complexity of Bond Failures in Graph States

**Ayush Nadiger**

Reproducibility repository for the August 2026 v10 manuscript **“Repair Identifiability and Product-Measurement Complexity of Bond Failures in Graph States.”**

The paper studies fault diagnosis for a **known graph state** when an intended bond is suspect. It separates two questions that ordinary target-state verification can conflate: **which bond failed?** and **which physical repair is required?** The central support theorem shows that, for a severed edge versus complete endpoint `Z`-dephasing, a Pauli observable distinguishes the two mechanisms exactly when

```text
P ∈ ±S(G−e) \ ±S(G).
```

Thus the repair information lives precisely in stabilizers of the failed graph that are absent from the target stabilizer dictionary. The manuscript then develops measurement-context complexity, bounded-Pauli-weight separations, locality criteria, and finite-shot/noise scaling for recovering that information.

## Status

- Manuscript: **v10, August 2026**
- arXiv: **not yet posted**
- Project page: https://ayushnadiger.github.io/projects/graph-state-repair-diagnosis.html

The hidden-cut / StateHSP literature is closely related in motivation but solves a different inference problem: those works learn an unknown factorization, symmetry, or stabilizer structure, whereas this paper assumes the graph topology is known and asks which **physical fault mechanism** occurred at an intended edge under a restricted measurement-and-retention interface.

## Repository layout

```text
paper/
  main.tex              LaTeX manuscript
  main.pdf              compiled v10 manuscript
repro/
  ...                   exact and independent verification scripts / data
ITERATION_AUDIT_v10.md  adversarial freeze-pass audit
CITATION.cff            citation metadata
```

Finite-enumeration claims are labeled **computer-assisted** in the manuscript. The analytic theorems do not depend on software.

## Quick verification

A representative lightweight verification pass can be run from `repro/`:

```bash
python audit_single_edge.py
python audit_endpoint_tomography.py
python audit_compression_gap.py
python audit_w8_weighted_fast.py
python wheel_cert_check.py
python audit_wheel_weight3_symbolic_boundaries.py
python audit_kn_weight_profile.py
```

The Python checks use the standard library plus `numpy` and/or `networkx` where indicated. The exact `W14` obstruction search uses the included C++17/OpenMP sources and is intentionally separate from the lightweight checks because it enumerates `3^14 = 4,782,969` product-Pauli contexts and creates a large intermediate file.

## Reproducibility scope

The repository contains independent implementations of the main finite computations, including:

- single-edge target-diagonal blindness and failed-graph repair witnesses;
- endpoint-tomography classification;
- raw-record versus target-stabilizer compression gaps;
- exact `W8` bounded-weight hierarchy;
- `W10` raw-repair obstruction and cross-checks;
- symbolic low-weight wheel-family checks;
- complete-graph weight profiles;
- connected-graph censuses through eight vertices;
- C++17/OpenMP sources for the exact `W14` raw-repair obstruction.

For the large finite computations, the manuscript and `ITERATION_AUDIT_v10.md` state exactly which conclusions are computer-assisted and which are proved algebraically.

## Citation

Citation metadata is in [`CITATION.cff`](CITATION.cff). Once the arXiv identifier is assigned, this repository and the project page will be updated to use the arXiv record as the preferred paper citation.

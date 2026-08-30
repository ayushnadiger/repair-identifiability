# Reproducibility snapshot — v10

Exact finite checks for **Repair Identifiability and Product-Measurement Complexity of Bond Failures in Graph States**.

The analytic theorems do not depend on software. Computer-assisted propositions and finite censuses are explicitly labeled in the manuscript.

## Quick checks

Run from this directory:

```bash
python audit_single_edge.py
python audit_endpoint_tomography.py
python audit_compression_gap.py
python audit_w8_weighted_fast.py
python wheel_cert_check.py
python audit_kn_weight_profile.py
```

These cover the single-edge blindness/cure identity, the endpoint-pair tomography criterion, raw-record compression examples, the exact W8 bounded-weight hierarchy, symbolic wheel certificates, and the complete-graph weight profile.

## W10 exact obstruction

```bash
python audit_w10_repair.py
```

Expected core output: all `3^10 = 59,049` contexts, `23,698` distinct nonzero separation masks, no two-context repair cover, and 120 independent representation cross-checks.

## W14 exact obstruction

The W14 computation is deliberately split into enumeration and search because the generated intermediate is large:

```bash
g++ -O3 -march=native -fopenmp -std=c++17 audit_w14_repair_enumerate.cpp -o audit_w14_repair_enumerate
./audit_w14_repair_enumerate 14 12 8

g++ -O3 -march=native -fopenmp -std=c++17 audit_w14_repair_search.cpp -o audit_w14_repair_search
./audit_w14_repair_search
```

Expected enumeration: `distinct collision masks=1260824`.
Expected exact search:

```text
NO_TWO_COVER W14 distinct_collision_masks=1260824 groups=36822 checks=217952792
```

The generated `w14_masks.bin` is about 212 MB and is intentionally excluded from Git.

## Dependencies

Python checks use the standard library plus `numpy` and/or `networkx` where imported. C++ checks require a C++17 compiler; OpenMP is used for the large obstruction computation.

This repository is the browseable public verification snapshot. The frozen arXiv source/reproducibility bundle will be linked here once the manuscript receives an arXiv identifier.

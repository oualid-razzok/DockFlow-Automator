# Redocking benchmark results

- complexes attempted: 1
- completed: 1
- successful (best pose RMSD <= 2.0 A): 1/1 (100%)
- exhaustiveness: 16 (fixed seed 2026)

| pdb | ligand | target | best affinity | best RMSD (A) | success |
|---|---|---|---|---|---|
| 1HVR | XK2 | HIV-1 protease (cyclic urea inhibitor) | -11.13 | 1.08 | yes |

Row 1 was produced by the committed example run (identical protocol:
`examples/configs/hiv1_protease_example.yaml`, seed 2026, exhaustiveness
16, pocket-derived box; full artifacts in `examples/results/1HVR/`).
Run `python benchmarks/redocking/run_benchmark.py` to extend the table
with the remaining 23 complexes (~2-4 CPU-hours at exhaustiveness 8).

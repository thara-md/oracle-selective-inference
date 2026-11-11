# Minimum-evidence benchmark (P3, November 2026)

Independent 400-dataset synthetic experiment (100 datasets for each of four block structures; n=120, 36 errors, three alpha budgets and 17 coverage targets; total 20,400 comparisons). No patient-level inputs. The code also emits a four-row **aggregate** clinical summary (no patient rows or images).

From this directory, with Python >=3.10 and `numpy`, `scipy`, `matplotlib`:

```bash
python P3_minimal_evidence_audit.py
python P3_verify_and_figure.py
```

Expected at 5% selective risk: 2,167 case-wise-feasible but information-infeasible dataset–target pairs, 479 certified by one block (22.1%), median 3 blocks. `P3_synthetic_summary.csv` is an archived reference. `P3_experiment_manifest.json` fixes seeds, generator and evaluation conventions. The verification performs separate linear-program spot checks and creates Figure 3.

**Provenance:** recovered from `P3_Oracle_最小証拠量実験_再現パッケージ.zip` inside the private archive `必勝膵oracle原稿.zip`. This is a distinct experiment from the September 16 structural benchmark; their synthetic populations and frequencies are not interchangeable.

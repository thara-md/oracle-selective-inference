# Synthetic structural benchmark (September 16, 2026)

This directory contains a **synthetic-only entrypoint adaptation** of the archived 2026-09-16 study code; neither patient-level records nor raw LLM transcripts are included. The generator and numerical analysis are unchanged; only the command-line entrypoint skips loading a separate institutional clinical input table.

From this directory, with Python >=3.10 and `numpy`, `pandas`, `scipy`:

```bash
python run_synthetic_only.py --phase synthetic --out outputs
```

The fixed `config.json` and archived pre-experiment config SHA in `execution_spec_lock.json` define 400 synthetic datasets: n=100/500, four error/collision structures, 50 seeds per combination. Outputs are written under `outputs/`. The main 5.72% / 5.24% disclosure-limited exclusion results are independently reproduced from these outputs. The reference CSVs are summaries only; regenerate full synthetic datasets locally. Timestamped validation JSON is not expected to match byte-for-byte.

**Provenance:** source `案2_追加検証②③A_実行結果一式_20260916.zip`; clinical-only entrypoint loading disabled for safe distribution. Exact source byte identity is **not** asserted for this adapted script. Synthetic generation and fraction-arithmetic routines remain the archived implementation. This is a conditional fixed-sample stress test, not a clinical deployment validation.

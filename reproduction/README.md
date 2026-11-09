# Reproducibility addendum — Oracle selective inference

**Status:** Verified partial replication package, 2026-11-09. This is **not** a complete reproduction of every experiment in the manuscript. The manuscript and published figure assets are maintained under `paper/` at the repository root.

## Included and tested

- **Clinical aggregate collision oracle (P0)**: `P0_clinical_oracle_200seed_repro.py` reproduces the clinical aggregate-block oracle control with 200 random seeds for four strata and three noise cardinalities (2,400 exact arithmetic checks). It uses only reconstructed *counts of block sizes and errors*, not patient-level files. Current packaged tables match this code exactly.
- **Independent synthetic experiment (P1)**: `_work_p1/synthetic/oracle_negative_control.py` reproduces 200 paired seeds of train/calibration/test controls; the JSON generated with the exact CLI below matches the supplied archived output byte for byte.
- **Figures 1 and 2 audit (P2/P1)**: `P2_figure1_exact_audit.py` and `P1_figure2_and_audit.py` regenerate the described plots and check figure numbers from their aggregate inputs / synthetic results.

## Setup

```bash
python -m pip install -r requirements.txt
```

Run all commands from this `reproduction/` directory:

```bash
python P0_clinical_oracle_200seed_repro.py
python _work_p1/synthetic/oracle_negative_control.py --seeds 200 --alpha 0.05 --output _work_p1/synthetic/oracle_negative_control_results.json
python P2_figure1_exact_audit.py
python P1_figure2_and_audit.py
```

These commands overwrite derived outputs in `reproduction/`, not patient records. The synthetic Python program's original *default* output is `/mnt/data/oracle_negative_control_results.json`; use `--output` as shown to keep results local to the repo.

## Important gaps / provenance

**Not included:**

1. The source generator / exhaustive verification for the first 400-dataset synthetic structural study and the distinct 400-dataset, 20,400 dataset–risk–coverage minimum-certificate study (Section 3.6 / Results 4.2). This addendum cannot reproduce those reported results end-to-end.
2. Raw archived `GPTBbeta.zip` and `GeminiBbeta.zip` files needed to rerun original 174-record label extraction; source-level case provenance must remain institutional/private.
3. Patient-level imaging, clinical data, raw model run transcripts, paths, and linkage tables. Do NOT publish those without institutional authorization.
4. Prospective or held-out *clinical* evidence; the independent test study here is synthetic.

The P0 computation is rederived from audited aggregate equivalence blocks, **not** from the original patient-level sequences. The first 50-seed RNG stream is not recovered; 200 seeds are an independent rerun. The figures have fixed-sample / in-sample interpretation only.

`P0_clinical_oracle_200seed_repro.py` is the original computational script with only two case-number references removed from comments. No algorithms, inputs, or numerical outputs were changed. All other scripts here are retained from the supplied packages. The September and November archival ZIPs remain private and unchanged.

## Safe data sharing

This directory was assembled from scripts, aggregate block statistics, and synthetic results. Patient-level CSVs, raw clinical transcripts, and original input archives have intentionally been excluded. Before converting this repository from Private to Public, review institutional data-sharing and workshop anonymity rules and perform a second repository-wide privacy/secret scan. A `.gitignore` **does not** remove files already committed.

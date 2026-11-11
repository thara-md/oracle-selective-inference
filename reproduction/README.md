# Reproducibility addendum — Oracle selective inference

**Status:** Updated 2026-11-11: the two previously absent 400-dataset synthetic experiments have been recovered and added. Reproduction remains partial for restricted patient-level clinical data. The manuscript and figure assets are maintained under `paper/` at the repository root.

## Included and tested

- **Clinical aggregate collision oracle (P0)**: `P0_clinical_oracle_200seed_repro.py` reproduces the clinical aggregate-block oracle control with 200 random seeds for four strata and three noise cardinalities (2,400 exact arithmetic checks). It uses only reconstructed *counts of block sizes and errors*, not patient-level files. Current packaged tables match this code exactly.
- **Independent synthetic experiment (P1)**: `_work_p1/synthetic/oracle_negative_control.py` reproduces 200 paired seeds of train/calibration/test controls; the JSON generated with the exact CLI below matches the supplied archived output byte for byte.
- **Figures 1 and 2 audit (P2/P1)**: `P2_figure1_exact_audit.py` and `P1_figure2_and_audit.py` regenerate the described plots and check figure numbers from their aggregate inputs / synthetic results.
- **September structural 400-dataset experiment**: [`legacy_structural_400/`](legacy_structural_400/README.md) contains the fixed generator settings, a synthetic-only entrypoint, and reference summaries; generated outputs reproduce the reported 5.72% and 5.24% figures. The original command-line entrypoint was adapted **only to bypass loading restricted clinical inputs**; see its README.
- **November P3 minimum-evidence 400-dataset experiment**: [`minimum_evidence_400/`](minimum_evidence_400/README.md) includes the experiment manifest, original computational code, independent LP checks and reference summary. At 5% risk, it reproduces 2,167 information-infeasible dataset–target pairs, including 479 one-block certificates (22.1%), among 20,400 overall comparisons. 
- **Two-experiment map**: [`README_400_EXPERIMENTS.md`](README_400_EXPERIMENTS.md) describes the differing populations, methods and commands.

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

For the two additional 400-dataset experiments, run from their **respective directories**, with the dependencies stated in each local README:

```bash
cd legacy_structural_400
python run_synthetic_only.py --phase synthetic --out outputs
# Return to reproduction/, then:
cd ../minimum_evidence_400
python P3_minimal_evidence_audit.py
python P3_verify_and_figure.py
```

The corresponding reference CSVs allow comparison with regenerated summaries; the full synthetic records need not be committed.

## Important gaps / provenance

**Not included:**

1. Raw archived `GPTBbeta.zip` and `GeminiBbeta.zip` files needed to rerun original 174-record label extraction; source-level case provenance remains institutional/private.
2. Patient-level imaging, clinical data, raw model run transcripts, paths, and linkage tables. Do NOT publish those without institutional authorization.
3. Prospective or held-out *clinical* evidence; the independent test study here is synthetic.
4. Independent reauthentication of all original patient-level inclusion and clinical GT decisions; synthetic reproducibility does not resolve these separate clinical-provenance limitations.

The P0 computation is rederived from audited aggregate equivalence blocks, **not** from the original patient-level sequences. The first 50-seed RNG stream is not recovered; 200 seeds are an independent rerun. The figures have fixed-sample / in-sample interpretation only.

`P0_clinical_oracle_200seed_repro.py` is the original computational script with only two case-number references removed from comments. Its algorithms, inputs, and numerical outputs were unchanged. The September structural 400-dataset entrypoint is a synthetic-only adaptation of its original code, as clearly documented in that subdirectory; the archived originals remain private and unchanged. Other included scripts retain the provenance documented in their respective READMEs.

## Safe data sharing

This directory was assembled from scripts, aggregate block statistics, and synthetic results. Patient-level CSVs, raw clinical transcripts, and original input archives have intentionally been excluded. Before converting this repository from Private to Public, review institutional data-sharing and workshop anonymity rules and perform a second repository-wide privacy/secret scan. A `.gitignore` **does not** remove files already committed.

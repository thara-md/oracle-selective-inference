# Two reproducible synthetic 400-dataset audits

The repository now has two **different** synthetic experiments in `reproduction/`:

1. [`legacy_structural_400/`](legacy_structural_400/README.md) — 2026-09-16 partial collision-disclosure structural benchmark, n=100/500, 50 seeds × 4 structures each; preserves the 5.72% / 5.24% results for distributed errors. Its archived generator is supplied with a synthetic-only CLI adaptation because the original CLI required a separate restricted institutional clinical input file even with `--phase synthetic`. Source algorithm and JSON configuration are preserved; regenerated outputs byte-match the archived synthetic outputs.
2. [`minimum_evidence_400/`](minimum_evidence_400/README.md) — separate 2026-11 P3 minimum certificate experiment, n=120, 36 errors, four structures ×100 datasets, 3 risk budgets ×17 coverage targets, 20,400 comparisons; reproduces 2,167 information-infeasible pairs at α=5%, with 479 single-block certificates (22.1%).

Both use **artificial/synthetic observations**. These folders exclude patient-level images, GT joins, clinical run transcripts, and institutional data tables. They establish finite-sample synthetic reproducibility, not future-patient risk control. Clinical data provenance has different access restrictions; see `../provenance/README.md`.

Run each experiment from its own directory as described in its README. For institutional sharing/release, recheck repository-wide security, licenses and workshop submission policy.

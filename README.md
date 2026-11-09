# Separating Threshold Failure from Information Limits in Repeated LLM Inference

**Repository preparation package — private draft, 9 November 2026.**

Authors: Takeshi Hara and Hideya Konishi.

This project accompanies a draft for W3PHIAI-27 on retrospective, fixed-sample auditing of selective prediction with repeated diagnostic-label outputs. It separates the failure of a particular uncertainty threshold from information-restricted infeasibility using verified observation-equality blocks. It also distinguishes empirical post-hoc oracle improvements from independent-test performance.

## Contents

- `paper/oracle_aaai27_trial.tex` — manuscript LaTeX source (trial typesetting).
- `paper/oracle_refs.bib` — references.
- `paper/figure1.pdf` and `paper/figure2.png` — manuscript figures as supplied.
- `paper/aaai2027.sty` and `paper/aaai2027.bst` — style files supplied in the author's draft ZIP; use subject to the publisher's instructions.
- `paper/manuscript_review_2026-11-09.pdf` — author-provided review PDF (8 pages).
- `code/partial_collision_certificate.py` — **new, standalone illustrative implementation** of the exact-rational partial-equivalence upper bound for the aggregate numerical example in Section 4.2.
- `tests/test_partial_collision_certificate.py` — tests for selected published aggregate rational values.
- `docs/original_packaging_notes.txt` — original source-ZIP notes, retained for provenance; these describe an earlier 7-page trial and should not be read as current submission status.

## Run the illustrative exact-rational certificate

Python 3.10+ (standard library only):

```bash
python code/partial_collision_certificate.py
python -m unittest discover -s tests -v
```

This verifies the aggregate all-category GPT example at 5% selective risk: the case-wise upper bound (`550/817`), the PDACx20 partial equality certificate (`435/946`), the two-block bound (`210/473`), and minimum candidate-block counts of 1 and 2 for coverage targets 50% and 45%, respectively. These are **not patient-level or prospective tests**.

## Compiling the trial manuscript

Open `paper/oracle_aaai27_trial.tex` with pdfLaTeX + BibTeX (including the adjacent style and figure files), following the conference's current instructions. **Do not treat the `preprint` format option or this document as finalized for submission.**

## Data, provenance, and scope limitations

The underlying retrospective clinical study received IRB approval (No. 7224). Patient-level imaging, clinical data, individual model outputs, and patient identifiers are intentionally **not distributed**. The original experimental analysis scripts and simulation pipelines were **not supplied in the material used to create this package**, so this repository does **not** reproduce all results reported in the manuscript. The separate certificate illustration was written for this repository and is not represented as the archived study implementation.

Results in the draft refer to a fixed labeled cohort and retrospective oracle audits. The repository does not provide prospective clinical risk-control guarantees or a deployable diagnostic system.

## Before publishing or citing

Confirm the final workshop call for papers, review process, public repository policy, appropriate manuscript style/copyright status, bibliography metadata, code provenance, and the permitted scope of institutional sharing. The draft's previous packaging notes identify additional author checks. If this is made public, review all contents for protected health information, nonpublic prompts/logs, API tokens, and account secrets.

No open-source or document redistribution license has yet been assigned to this repository. Permission for reuse should be obtained from the relevant rights holders.

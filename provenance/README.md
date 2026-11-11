# Restricted clinical output provenance (2026-11-09)

This folder **does not contain** individual patient-level data, free-text model outputs, the original run archives, the original case-level ground-truth labels, or any source imaging. The `pancreasstudyBbeta(2).zip` bundle must remain under institutionally approved restricted storage; **do not upload it** to this repository (including a private GitHub repository) without explicit authorization.

## Recovered provenance

The new `pancreasstudyBbeta(2).zip` contains the historic `GPTBbeta.zip` and `GeminiBbeta.zip` and corresponding extracted `.txt` files. The archived provenance `input_provenance.json` dated September 15, 2026 specifies the exact SHA-256 values retained in `verify_source_archives.py`. Both embedded ZIP byte hashes match the archive record, and every extracted text file matches its corresponding embedded ZIP member.

| Property | GPT | Gemini |
| --- | ---: | ---: |
| Archived run files | 1,791 | 1,806 |
| SHA-256 matches historic record | Yes | Yes |

Total: 3,597 raw text files, including 19 flagged files (`mistake` or `missing`) in the combined corpus. Flags are administrative run markers, not licenses to silently alter cohort selection.

Run the metadata-only verification:

```shell
python provenance/verify_source_archives.py /restricted/path/pancreasstudyBbeta.zip
```

## Retrospective reproduction (performed locally; NOT included here)

With the original September 15, 2026 `recover_and_audit.py`, its preserved private `inputs/` and the two SHA-256-verified original ZIP files, the 174-unit / 3,480-label audit was rerun on November 9, 2026. All 16 archived result CSV files matched byte-for-byte. The numerical audit summary matched; its input paths differ because of the local directory structure.

Important qualifications, preserved from the historic audit:

- 52 of 52 previously missing unit sequences were recovered; stored entropy matched all 174 units.
- For 65 units the historic accepted-run list was available. For 108 units, there were exactly 20 unflagged candidates but no independent acceptance manifest certification.
- At the **September 15 archive-only recovery stage**, GPT case005 had 21 unflagged candidate runs sharing the same diagnosis; the recorded PDAC×20 label sequence was invariant to which 20 were selected. The individual accepted-run subset was **not independently identified from those recovered files alone**. A later documentary reconciliation is described below.
- Seven saved modal diagnoses were tied; the stored diagnoses were retained, not re-tie-broken.
- This retrospective verification does not independently certify image submission, backend model identities, final GT pathology, IRB coverage, or prospective risk control.

## Subsequent case005 documentary reconciliation (September 17, 2026)

A later author-compiled `案2_Supplementary_S1_case005と来歴_20260917.docx` reconciled the earlier unresolved *file-subset identification* with existing written run-inclusion rules and an individual case record:

- `Shijisho_v3_3_20260730`, §2.1 R-4, specifies retaining the first 20 valid runs in ascending run-number order when more than 20 valid candidates exist.
- `Amendment4_E2_20260730`, §3.1, records GPT case005 Run21 as an extra valid run to be excluded.
- On that **documentary basis**, Supplementary S1 identifies **Run01–20 as accepted and Run21 as excluded**. All 21 candidates have the same diagnostic label (PDAC), so the 20-label Oracle input is unchanged.

This is a **subsequent documentation-level reconciliation**, not an alteration of the September 15 raw audit or an independent reauthentication of model inputs, images, pathology, or the original per-file manifest. Preserve both dated stages rather than retroactively rewriting historical evidence. The cited original institutional records and Supplementary S1 are not distributed here.

**Do not misrepresent `verify_source_archives.py` as reproducing the full Oracle experiments.** It verifies byte integrity of archived inputs. The original restricted `recover_and_audit.py` and private input tables are separately needed for the full clinical label audit. Recovered synthetic-only generation and verification scripts are documented under [`../reproduction/`](../reproduction/README.md).

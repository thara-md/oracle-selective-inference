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
- For GPT case005, 21 unflagged runs share the same diagnosis. The 20-label sequence is invariant, but the historic selected subset of 20 original run files cannot be identified from the archived evidence.
- Seven saved modal diagnoses were tied; the stored diagnoses were retained, not re-tie-broken.
- This retrospective verification does not independently certify image submission, backend model identities, final GT pathology, IRB coverage, or prospective risk control.

**Do not misrepresent `verify_source_archives.py` as reproducing the full Oracle experiments.** It verifies byte integrity of archived inputs. The original restricted `recover_and_audit.py` and private input tables are separately needed for the full clinical label audit. Other synthetic experiments require their own scripts and seeds.

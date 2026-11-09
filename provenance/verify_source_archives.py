#!/usr/bin/env python3
"""Verify the integrity of restricted original run archives without extracting data.

Usage:
  python provenance/verify_source_archives.py /path/to/pancreasstudyBbeta.zip

Prints only aggregate metadata, never archived clinical output text.
"""
import hashlib
import io
import json
import re
import sys
import zipfile
from pathlib import Path

EXPECTED = {
    'GPT': {'sha256': '2d6b393c62024781d6347c3fe32798debaf072e28bbeb25ae1347a6c80106a78', 'count': 1791},
    'Gemini': {'sha256': '0b3e83fd0802f215ca05960190601724737ef81cff7a427c60a755ee4a0810ed', 'count': 1806},
}


def verify(path: Path) -> dict:
    results = {}
    with zipfile.ZipFile(path) as outer:
        names = outer.namelist()
        for model, expected in EXPECTED.items():
            filename = f'{model}Bbeta.zip'
            matches = [n for n in names if n.endswith('/raw_runB/' + filename)]
            if len(matches) != 1:
                raise ValueError(f'{filename}: expected one nested archive, found {len(matches)}')
            inner_bytes = outer.read(matches[0])
            actual_hash = hashlib.sha256(inner_bytes).hexdigest()
            if actual_hash != expected['sha256']:
                raise ValueError(f'{filename}: SHA-256 differs from archived provenance')
            folder_prefix = matches[0][:-4] + '/'
            direct_files = [n for n in names if n.startswith(folder_prefix) and n.endswith('.txt')]
            with zipfile.ZipFile(io.BytesIO(inner_bytes)) as inner:
                inner_names = [n for n in inner.namelist() if n.endswith('.txt')]
                direct_hashes = {n.rsplit('/', 1)[-1]: hashlib.sha256(outer.read(n)).hexdigest() for n in direct_files}
                inner_hashes = {n.rsplit('/', 1)[-1]: hashlib.sha256(inner.read(n)).hexdigest() for n in inner_names}
            if direct_hashes != inner_hashes:
                raise ValueError(f'{filename}: extracted folder and embedded ZIP do not match')
            if len(inner_hashes) != expected['count']:
                raise ValueError(f'{filename}: unexpected run count {len(inner_hashes)}')
            flagged = sum(bool(re.search(r'_(?:missing|mistake)(?:_missing|_mistake)?\.txt$', name)) for name in inner_hashes)
            results[model] = {'run_files': len(inner_hashes), 'flagged_files': flagged, 'sha256': actual_hash,
                              'embedded_zip_and_folder_identical': True}
    if sum(x['flagged_files'] for x in results.values()) != 19:
        raise ValueError('Flagged file count differs from recorded audit (19)')
    return {'verification': 'PASS', 'total_raw_run_files': sum(x['run_files'] for x in results.values()), 'models': results,
            'notice': 'Integrity check only; does not verify patient-level GT, original imaging submission or analysis approvals.'}


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print(__doc__.strip(), file=sys.stderr)
        raise SystemExit(2)
    print(json.dumps(verify(Path(sys.argv[1])), indent=2))

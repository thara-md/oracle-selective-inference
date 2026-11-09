#!/usr/bin/env python3
"""P0 clinical random-representation control, 200 Monte Carlo seeds.

Does not change patient diagnoses, GT or any original files.
Inputs are fully identified (block size, errors) from 2026-09-15 audit report,
reconciled with 2026-11-04 tables 1, 2, 4 and Appendix A2.
It is an aggregate-block reconstruction, NOT a raw patient-level re-run.
"""
import csv
import hashlib
import json
import math
from fractions import Fraction
from pathlib import Path
from statistics import mean, stdev
import numpy as np

OUTPUT_DIR = Path(__file__).resolve().parent
BASE_SEED = 20261104
N_SEEDS = 200
KS = (2, 4, 16)
ALPHA = Fraction(1, 20)

# Each duplicated 20-output sequence forms a single block; all other
# sequences are singletons. Each tuple is (block size, # incorrect predictions).
# GPT non-unanimous 2-record block: both correct.
# Gemini non-unanimous 2-record block: both incorrect.
BLOCKS = {
    'All-category GPT': [(32, 6), (2, 1), (2, 0)] + [(1, 0)] * 26 + [(1, 1)] * 24,
    'All-category Gemini': [(11, 2), (14, 8), (2, 2)] + [(1, 0)] * 23 + [(1, 1)] * 38,
    'PDAC/NET GPT': [(28, 2), (2, 1), (2, 0)] + [(1, 0)] * 26 + [(1, 1)] * 17,
    'PDAC/NET Gemini': [(11, 2), (13, 7), (2, 2)] + [(1, 0)] * 22 + [(1, 1)] * 29,
}
BLOCK_LABELS = {
    'All-category GPT': ['PDAC×20', 'NET×20', '19PDAC+1NET shared sequence'],
    'All-category Gemini': ['PDAC×20', 'NET×20', '19NET+1PDAC shared sequence'],
    'PDAC/NET GPT': ['PDAC×20', 'NET×20', '19PDAC+1NET shared sequence'],
    'PDAC/NET Gemini': ['PDAC×20', 'NET×20', '19NET+1PDAC shared sequence'],
}
EXPECTED = {
    'All-category GPT': {'n': 86, 'errors': 31, 'fraction': '210/473', 'lambda': '80/11', 'old_50': [50.92, 55.87, 64.85]},
    'All-category Gemini': {'n': 88, 'errors': 50, 'fraction': '115/319', 'lambda': '220/29', 'old_50': [37.11, 39.66, 43.87]},
    'PDAC/NET GPT': {'n': 75, 'errors': 20, 'fraction': '104/135', 'lambda': '20/9', 'old_50': [77.11, 77.16, 77.19]},
    'PDAC/NET Gemini': {'n': 77, 'errors': 40, 'fraction': '80/203', 'lambda': '220/29', 'old_50': [40.79, 44.21, 49.06]},
}

def exact_lp(blocks, alpha=ALPHA):
    """Exact primal and supporting dual for coverage maximization at fixed risk.

    All optimization and primal-dual validity checks use fractions.Fraction;
    floating-point appears only when aggregating/reporting Monte Carlo means.
    """
    N = sum(n for n, e in blocks)
    assert N > 0 and all(type(n) is int and type(e) is int and n >= 1 and 0 <= e <= n for n, e in blocks)
    q = [Fraction(0) for _ in blocks]
    budget = Fraction(0)
    for j, (n, e) in enumerate(blocks):
        cost = Fraction(e) - alpha * n
        if cost <= 0:
            q[j] = Fraction(1)
            budget -= cost
    costly = sorted((j for j, (n, e) in enumerate(blocks) if Fraction(e) - alpha * n > 0),
                    key=lambda j: Fraction(blocks[j][1], blocks[j][0]))
    dual_lambda = Fraction(0)
    for j in costly:
        n, e = blocks[j]
        cost = Fraction(e) - alpha * n
        if budget >= cost:
            q[j] = Fraction(1)
            budget -= cost
            if budget == 0:
                dual_lambda = Fraction(n, 1) / cost
        else:
            if budget > 0:
                q[j] = budget / cost
            dual_lambda = Fraction(n, 1) / cost
            budget = Fraction(0)
            break
    accepted = sum(Fraction(n) * q[j] for j, (n, e) in enumerate(blocks))
    errors = sum(Fraction(e) * q[j] for j, (n, e) in enumerate(blocks))
    assert accepted == 0 or errors <= alpha * accepted
    assert all(Fraction(0) <= x <= 1 for x in q)
    primal = accepted / N
    dual = sum(max(Fraction(0), Fraction(n) - dual_lambda * (Fraction(e) - alpha * n))
               for n, e in blocks) / N
    assert primal == dual, (primal, dual, blocks)
    return primal, dual_lambda, len(blocks)

def refined_blocks(blocks, k, rng):
    """Independent k-ary noise per individual; group by original block and noise.

    Assignments are independent of the fixed outcome status. Grouping correct
    and incorrect items is only an exchangeable computational representation.
    """
    out = []
    for n, e in blocks:
        syms = rng.integers(0, k, size=n)
        # First e are incorrect for exchangeable sampling; not supplied to selector.
        for s in range(k):
            indices = (syms == s)
            ns = int(np.count_nonzero(indices))
            if ns:
                es = int(np.count_nonzero(indices[:e]))
                out.append((ns, es))
    assert sum(n for n, e in out) == sum(n for n, e in blocks)
    assert sum(e for n, e in out) == sum(e for n, e in blocks)
    return out

def write_csv(path, data, columns):
    with path.open('w', newline='', encoding='utf-8-sig') as f:
        w = csv.DictWriter(f, fieldnames=columns, extrasaction='ignore')
        w.writeheader(); w.writerows(data)


def main():
    summary, per_seed, inputs = [], [], []
    details = {}
    for condition_index, (name, blocks) in enumerate(BLOCKS.items()):
        expected = EXPECTED[name]
        n = sum(b[0] for b in blocks)
        ne = sum(b[1] for b in blocks)
        oracle, lam, nblocks = exact_lp(blocks)
        assert (n, ne, str(oracle), str(lam)) == (expected['n'], expected['errors'], expected['fraction'], expected['lambda'])
        uid_oracle, _, _ = exact_lp([(1, 0)] * (n - ne) + [(1, 1)] * ne)
        # A transparent per-block input ledger, including singleton counts.
        for i, (bn, be) in enumerate(blocks[:3]):
            inputs.append(dict(condition=name, block_kind=BLOCK_LABELS[name][i], n=bn, errors=be,
                               is_duplicated=True))
        inputs.extend([
            dict(condition=name, block_kind='unique correct sequences (individual singletons)', n=n-ne-sum(bn-be for bn,be in blocks[:3]), errors=0, is_duplicated=False),
            dict(condition=name, block_kind='unique erroneous sequences (individual singletons)', n=ne-sum(be for bn,be in blocks[:3]), errors=ne-sum(be for bn,be in blocks[:3]), is_duplicated=False)
        ])
        details[name] = {'n': n, 'errors': ne, 'n_base_blocks': nblocks,
                         'full_sequence_fraction': str(oracle),
                         'full_sequence_coverage_pct': float(100*oracle),
                         'lambda': str(lam),
                         'casewise_fraction': str(uid_oracle),
                         'casewise_coverage_pct': float(100*uid_oracle),
                         'noise': {}}
        print(name, 'N=',n,'E=',ne, 'baseline exact=',oracle, 'lambda=',lam, 'UID=',uid_oracle)
        for k_index, k in enumerate(KS):
            values=[]; n_unique=[]
            for seed in range(N_SEEDS):
                rng = np.random.default_rng(np.random.SeedSequence([BASE_SEED, condition_index, k, seed]))
                new_blocks=refined_blocks(blocks,k,rng)
                val, ld, bcount=exact_lp(new_blocks)
                v=float(100*val)
                assert val >= oracle and val <= uid_oracle
                values.append(v); n_unique.append(bcount)
                per_seed.append(dict(condition=name, k=k, seed=seed, coverage_pct=v,
                                     exact_fraction=str(val), dual_multiplier=str(ld), n_refined_blocks=bcount,
                                     baseline_pct=float(100*oracle)))
            mu=mean(values); sd=stdev(values); se=sd / math.sqrt(len(values))
            entry = dict(condition=name, n=n, errors=ne, risk_pct=5, k=k, seeds=N_SEEDS,
                         baseline_pct=float(100*oracle), mean_pct=mu, gain_pct_points=mu-float(100*oracle),
                         sd_pct=sd,se_pct=se,mc_lower95_pct=mu-1.96*se,mc_upper95_pct=mu+1.96*se,
                         min_pct=min(values),max_pct=max(values),casewise_pct=float(100*uid_oracle),
                         old_50_mean_pct=expected['old_50'][k_index],
                         delta_vs_prior50_pct_points=mu-expected['old_50'][k_index],
                         avg_n_refined_blocks=mean(n_unique))
            summary.append(entry)
            details[name]['noise'][str(k)]=entry
            print(' k=',k,'mean=',round(mu,4),'SD=',round(sd,4),'SE=',round(se,4),'95% MC mean CI=',(round(mu-1.96*se,4),round(mu+1.96*se,4)), 'prior50=',entry['old_50_mean_pct'])
    base_sources = ['膵oracle_table1.xlsx', '膵oracle_table2.xlsx', '膵oracle_table3.xlsx', '膵oracle_table4.xlsx', '膵oracle_tableA2.xlsx', 'gemini-code-1793764914887_膵oracle_P0.py']
    hashes={p:hashlib.sha256((OUTPUT_DIR/p).read_bytes()).hexdigest() for p in base_sources if (OUTPUT_DIR/p).exists()}
    manifest={
      'analysis':'P0 clinical random-symbol refinement (aggregate block reconstruction)',
      'timestamp_label':'2026-11-04',
      'n_seeds': N_SEEDS,'seed_range':'0...199','base_seed':BASE_SEED,'rng':'numpy Generator(PCG64 via SeedSequence [base_seed, condition_index, k, seed])',
      'symbol_values':list(KS),'alpha_fraction':str(ALPHA),
      'data_provenance':'Full sequence collision and GT block aggregates from 2026-09-15 archival audit narrative plus user-provided tables, reconciled with 4 exact rational primal-dual certificates. No raw patient-level data or run files reprocessed; no clinical labels modified.',
      'original_50_seed_rng':'not recovered; 200-seed results are independent rerun, not continuation of the original 50 RNG states',
      'mc_uncertainty':'95% normal-approximation Monte Carlo confidence intervals of simulated mean conditional on the fixed clinical cohort; NOT patient-sampling confidence intervals',
      'reference_sha256':hashes,
      'block_counts':{k:v for k,v in BLOCKS.items()},
      'details': details,
      'summary': summary,
    }
    write_csv(OUTPUT_DIR/'P0_200seed_summary.csv',summary,list(summary[0].keys()))
    write_csv(OUTPUT_DIR/'P0_200seed_per_seed.csv',per_seed,list(per_seed[0].keys()))
    write_csv(OUTPUT_DIR/'P0_reconstructed_block_inputs.csv',inputs,list(inputs[0].keys()))
    (OUTPUT_DIR/'P0_200seed_manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    return manifest

if __name__=='__main__':main()

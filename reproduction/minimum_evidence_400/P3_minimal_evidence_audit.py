#!/usr/bin/env python3
"""P3 follow-up to Oracle audit: evidence burden for fixed-sample exclusion.

Independent, fully specified experiment. Uses reconstructed clinical collision
AGGREGATES, not raw patients; synthetic suite is new and is NOT the same as
legacy (unarchived) 400 synthetic sets mentioned in P2.

All risk-coverage optima and comparisons use exact fractions. The task is to
CERTIFY INFEASIBILITY under a fixed input representation. Actual policy
performance/clinical generalization is not evaluated here.
"""
from __future__ import annotations
import csv, json, hashlib
from fractions import Fraction as F
from itertools import combinations
from collections import Counter, defaultdict
from statistics import median, mean
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
BASE_SEED=20261105
N_PER_SCENARIO=100
N=120; ERRORS=36
ALPHAS=(F(1,50),F(1,20),F(1,10))
TARGETS=tuple(F(i,40) for i in range(16,33)) # 40% to 80% by 2.5 points


def exact_lp(blocks:list[tuple[int,int]], alpha:F)->F:
    """Exact optimal coverage, allowing shared random acceptance within blocks.
    Based on fractional knapsack with costs (errors-alpha*size); outcome labels
    are used to compute ORACLE bound, never execution-time selector input.
    """
    total_n=sum(n for n,e in blocks)
    assert total_n>0 and all(isinstance(n,int) and isinstance(e,int) and n>0 and 0<=e<=n for n,e in blocks)
    budget=F(0); acc=F(0)
    costly=[]
    for n,e in blocks:
        cost=F(e)-alpha*n
        if cost<=0:
            budget-=cost; acc+=n
        else:
            costly.append((F(e,n),n,cost))
    costly.sort(key=lambda row:row[0])
    for _,n,cost in costly:
        if budget>=cost:
            acc+=n; budget-=cost
        else:
            if budget:
                acc+=F(n)*budget/cost
            break
    return acc/total_n


def optimistic_bound(verified, N, E, alpha):
    assert all(n>=2 for n,e in verified)
    nr=N-sum(n for n,e in verified)
    er=E-sum(e for n,e in verified)
    assert 0<=er<=nr
    # Every unverified item is casewise separable (optimistic relaxation).
    blocks=list(verified)
    if nr-er:blocks.append((nr-er,0))
    if er:blocks.append((er,er))
    return exact_lp(blocks,alpha)


def synthetic_blocks(scenario:str, rng)->list[tuple[int,int]]:
    """Construct six/seven exact-collision blocks; other records singletons.
    All datasets N=120, ERRORS=36. Only collision structures vary.
    """
    if scenario=='concentrated_mixed':
        sizes=[int(x+rng.integers(-2,3)) for x in [34,12,8]]
        errs=[int(x+rng.integers(-1,2)) for x in [8,4,2]]
    elif scenario=='dispersed_mixed':
        sizes=[int(x+rng.integers(-1,2)) for x in [12,11,10,9,8,7,6]]
        errs=[int(x+rng.integers(-1,2)) for x in [3,3,3,2,2,2,1]]
    elif scenario=='homogeneous':
        # Six blocks that are each all-correct or all-wrong, E total fixed.
        jitter=int(rng.integers(-2,3))
        sizes=[18+jitter,12+int(rng.integers(-2,3)),10,10+int(rng.integers(-2,3)),8-jitter,8+int(rng.integers(-2,3))]
        errs=[sizes[0],0,sizes[2],0,sizes[4],0]
    elif scenario=='high_uniqueness':
        sizes=[int(rng.integers(2,5)) for _ in range(6)]
        errs=[1]*6
    else: raise KeyError(scenario)
    assert len(sizes)==len(errs)
    assert sum(sizes)<N and sum(errs)<=ERRORS
    assert all(0<=e<=n and n>=2 for n,e in zip(sizes,errs))
    return list(zip(sizes,errs))


def all_bounds(verified_blocks:list[tuple[int,int]], alpha:F, N:int, E:int):
    """Subset mask -> optimistic upper-bound fraction, exact."""
    m=len(verified_blocks)
    bounds={}
    for mask in range(1<<m):
        subset=[b for i,b in enumerate(verified_blocks) if mask & (1<<i)]
        bounds[mask]=optimistic_bound(subset,N,E,alpha)
    return bounds


def options_for_target(bounds, blocks, target):
    valid=[mask for mask,b in bounds.items() if b<target]
    if not valid:return None
    # Minimize number of verified equality blocks; then records verified.
    best=min(valid,key=lambda mask:(mask.bit_count(),sum(b[0] for i,b in enumerate(blocks) if mask & (1<<i)),mask))
    return (best.bit_count(),sum(b[0] for i,b in enumerate(blocks) if best & (1<<i)),best)


def prefix_k(bounds,indices,target):
    mask=0
    if bounds[0]<target:return 0,0
    size=0
    for j,i in enumerate(indices):
        mask|=1<<i;size+=1
        if bounds[mask]<target:return size,mask
    return None,mask


def scenario_rows():
    scenario_names=('concentrated_mixed','dispersed_mixed','homogeneous','high_uniqueness')
    experiments=[]; summary=[]; examples=[]; subset_bounds=[]
    for idx,scenario in enumerate(scenario_names):
      for seed in range(N_PER_SCENARIO):
        rng=np.random.default_rng(np.random.SeedSequence([BASE_SEED,idx,seed]))
        blocks=synthetic_blocks(scenario,rng)
        perm=list(map(int,rng.permutation(len(blocks))))
        by_size=sorted(range(len(blocks)),key=lambda i:(-blocks[i][0],i))
        m=len(blocks)
        for alpha in ALPHAS:
            bounds=all_bounds(blocks,alpha,N,ERRORS)
            full=bounds[(1<<m)-1]
            casewise=bounds[0]
            assert full<=casewise
            # The bound can only decrease as verified equality information grows.
            for mask in range(1<<m):
                for i in range(m):
                    if not (mask>>i)&1:
                        assert bounds[mask | (1<<i)]<=bounds[mask]
            if seed==0 and alpha==F(1,20):
                for mask,b in sorted(bounds.items()):
                    subset_bounds.append(dict(scenario=scenario,seed=seed,mask=mask,verified_blocks=mask.bit_count(),
                                              verified_records=sum(n for i,(n,e) in enumerate(blocks) if mask&(1<<i)),
                                              alpha=str(alpha),upper_exact=str(b)))
            for target in TARGETS:
                # The zero-block case is infeasible even with case-wise GT.
                trivial=casewise<target
                is_info=(casewise>=target and full<target)
                feasible=full>=target
                assert trivial or is_info or feasible
                if not is_info: continue
                best=options_for_target(bounds,blocks,target)
                assert best is not None and best[0]>0
                k_size,size_mask=prefix_k(bounds,by_size,target)
                k_rand,random_mask=prefix_k(bounds,perm,target)
                assert k_size is not None and k_rand is not None
                rows=[]
                for frac in (F(1,4),F(1,2),F(3,4)):
                    # Floor, explicit: 25% of 3 blocks means reveal 0 blocks.
                    count=(len(blocks)*frac.numerator)//frac.denominator
                    chosen=by_size[:count]
                    mask=sum(1<<i for i in chosen)
                    rows.append(int(bounds[mask]<target))
                experiments.append(dict(scenario=scenario,seed=seed,alpha=str(alpha),target=str(target),
                        N=N,E=ERRORS,n_blocks=m,full_exact=str(full),casewise_exact=str(casewise),
                        min_blocks=best[0],min_verified_records=best[1],optimal_mask=best[2],
                        largest_first_blocks=k_size,largest_first_records=sum(blocks[i][0] for i in by_size[:k_size]),
                        random_blocks=k_rand,random_records=sum(blocks[i][0] for i in perm[:k_rand]),
                        largest_first_25=rows[0],largest_first_50=rows[1],largest_first_75=rows[2],
                        largest_first_1=int(bounds[1<<by_size[0]]<target),
                        full_structure=str(blocks), random_order=str(perm)))
            summary.append(dict(scenario=scenario,seed=seed,alpha=str(alpha),casewise=str(casewise),full=str(full),
                duplicated_blocks=m,duplicate_records=sum(n for n,e in blocks)))
    return experiments,summary,subset_bounds


def clinical_rows():
    # Only audits of reconstructed 2026-09-15 clinical aggregates, no raw patient reread.
    clinical={
    'All-category GPT':(86,31,[('PDAC×20',32,6),('NET×20',2,1),('non-unanimous duplicate',2,0)]),
    'All-category Gemini':(88,50,[('PDAC×20',11,2),('NET×20',14,8),('non-unanimous duplicate',2,2)]),
    'PDAC/NET GPT':(75,20,[('PDAC×20',28,2),('NET×20',2,1),('non-unanimous duplicate',2,0)]),
    'PDAC/NET Gemini':(77,40,[('PDAC×20',11,2),('NET×20',13,7),('non-unanimous duplicate',2,2)]),
    }
    records=[]
    alpha=F(1,20);target=F(1,2)
    for name,(n,e,items) in clinical.items():
        blocks=[(x[1],x[2]) for x in items]
        bounds=all_bounds(blocks,alpha,n,e)
        best=options_for_target(bounds,blocks,target)
        full=bounds[7];casewise=bounds[0]
        classification='casewise-infeasible' if casewise<target else ('information-infeasible' if full<target else 'in-sample feasible')
        records.append(dict(condition=name,N=n,errors=e,classification=classification,
                            casewise_exact=str(casewise),full_exact=str(full),
                            min_blocks=best[0] if best else 'NA',min_verified_records=best[1] if best else 'NA',
                            minimum_evidence=(', '.join(items[i][0] for i in range(3) if best[2]&(1<<i)) if best else 'not exclusion'),
                            min_bound_exact=str(bounds[best[2]]) if best else 'NA',
                            unverified_items='reconstructed block aggregates; original per-record collision file not reread'))
    return records


def csv_out(name,rows):
    if not rows: return
    with (HERE/name).open('w',newline='',encoding='utf-8-sig') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]))
        w.writeheader();w.writerows(rows)

def results():
    # SciPy independent spot-check will run in a separate audit script.
    data,datasets,subset= scenario_rows()
    clinical=clinical_rows()
    csv_out('P3_synthetic_information_infeasible_targets.csv',data)
    csv_out('P3_synthetic_dataset_manifest.csv',datasets)
    csv_out('P3_subset_exact_audit_seed0.csv',subset)
    csv_out('P3_clinical_minimal_evidence.csv',clinical)
    # Aggregate only over truly information-infeasible target pairs.
    stats=[]
    for alpha in ALPHAS:
      for scenario in ('concentrated_mixed','dispersed_mixed','homogeneous','high_uniqueness','ALL'):
        r=[z for z in data if z['alpha']==str(alpha) and (z['scenario']==scenario or scenario=='ALL')]
        if not r:
            stats.append(dict(scenario=scenario,alpha=str(alpha),pairs=0,affected_datasets=0))
            continue
        affected=len(set((x['scenario'],x['seed']) for x in r))
        def pct(k):return 100*sum(x[k] for x in r)/len(r)
        stats.append(dict(scenario=scenario,alpha=str(alpha),pairs=len(r),affected_datasets=affected,
            min_1block_pct=100*sum(x['min_blocks']==1 for x in r)/len(r),
            min_block_median=median(x['min_blocks'] for x in r),
            min_records_median=median(x['min_verified_records'] for x in r),
            largest_first_median=median(x['largest_first_blocks'] for x in r),
            random_median=median(x['random_blocks'] for x in r),
            largest_1block_detection_pct=pct('largest_first_1'),
            largest_25pct_detection_pct=pct('largest_first_25'),
            largest_50pct_detection_pct=pct('largest_first_50'),
            largest_75pct_detection_pct=pct('largest_first_75'),
            mean_min_blocks=mean(x['min_blocks'] for x in r),
            mean_largest_first_blocks=mean(x['largest_first_blocks'] for x in r),
            mean_random_blocks=mean(x['random_blocks'] for x in r)))
    csv_out('P3_synthetic_summary.csv',stats)
    meta={'experiment':'P3 independent minimal-partial-evidence fixed-sample audit',
     'base_seed':BASE_SEED,'n_per_scenario':N_PER_SCENARIO,'scenario_names':['concentrated_mixed','dispersed_mixed','homogeneous','high_uniqueness'],
     'total_datasets':4*N_PER_SCENARIO,'N':N,'error_count':ERRORS,
     'alpha_exact':list(map(str,ALPHAS)),'coverage_targets_exact':list(map(str,TARGETS)),
     'randomness':'NumPy default_rng SeedSequence [BASE_SEED, scenario_index, seed]',
     'procedure':'Exact optimistic upper bounds for every subset of known duplicate blocks. Unverified records may be perfectly error-separated. Lexicographic minimal certificate by number of verified blocks, then number of records.',
     'casewise_info_infeasible_definition':'casewise coverage >= target and full-sequence coverage < target',
     'exclude_strict':'upper_bound < coverage_target; equality is feasible, not excluded',
     'verification_cost_proxy':'Records in verified equality blocks; this is an abstract workload proxy, not measured analyst time; GT of the entire audited set is still assumed known.',
     'false_certificate':'Impossible when all revealed equality constraints are true; checked monotonicity on every simulated subset.',
     'scope':'Synthetic generator newly created for P3, NOT original unarchived synthetic 400 sets in 2026-11-04 P2, no claim of reproducing earlier numerical percentages.'}
    (HERE/'P3_experiment_manifest.json').write_text(json.dumps(meta,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
    print('P3',len(data),'information-infeasible setting pairs; manifest',len(datasets),'dataset-risk rows')
    for r in stats:
        if r['alpha']=='1/20':print(r)
    print('clinical')
    for c in clinical:print(c)
    return data,stats,clinical

if __name__=='__main__': results()

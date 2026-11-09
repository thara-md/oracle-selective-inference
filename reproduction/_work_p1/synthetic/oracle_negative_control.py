#!/usr/bin/env python3
"""Oracle inflation versus held-out utility; synthetic audit, not clinical evidence.

Usage: python oracle_negative_control.py --seeds 200 --alpha 0.05
Outputs /mnt/data/oracle_negative_control_results.json by default.

Data-generating process:
 X=0 (40%) has P(error)=0.01; X=1 (60%) has P(error)=0.35.
 Z is an independent nuisance symbol in {0,...,k-1}, k=2,4,16.
Positive control: S is a binary noisy indicator of error with a specified
conditional accuracy (artificially observable; not a clinical feature).

Each replicate uses paired independent train/calibration/test samples.
Oracle is optimized *using training labels* and may not generalize.
Deployable score: train-only smoothed group error; train-selected score order;
calibration-only risk/coverage threshold with randomized boundary; test untouched.
Calibration is empirical and is NOT a finite-sample risk guarantee.
"""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.optimize import linprog


def draw(rng, n):
    x = (rng.random(n) >= 0.4).astype(np.int8)
    err = (rng.random(n) < np.where(x == 0, 0.01, 0.35)).astype(np.int8)
    return x, err


def nuisance(rng, n, k):
    return rng.integers(0, k, n, dtype=np.int16)


def positive_proxy(rng, e, q):
    # q=Pr(S=E): research-only artificial signal generated with outcome.
    return np.where(rng.random(len(e)) < q, e, 1-e).astype(np.int8)


def aggregate(x, e, feature=None):
    key = x.astype(np.int32) if feature is None else (x.astype(np.int32)*100 + feature.astype(np.int32))
    uniq, inv = np.unique(key, return_inverse=True)
    n = np.bincount(inv).astype(float)
    m = np.bincount(inv, weights=e).astype(float)
    return uniq, n, m


def oracle_coverage(x, e, alpha, feature=None):
    _, n, m = aggregate(x, e, feature)
    res = linprog(-n, A_ub=[m-alpha*n], b_ub=[0.0], bounds=(0,1), method='highs')
    if not res.success:
        raise RuntimeError(res.message)
    a = res.x
    if np.dot(m-alpha*n, a) > 1e-6:
        raise ValueError('Oracle constraint violation')
    return float(np.dot(n, a)/sum(n))


def scores_from_training(x_train, e_train, x_eval, z_train=None, z_eval=None, lam=20.0):
    # Shared baseline group empirical risk, Laplace smoothed.
    base = {}
    for b in (0,1):
        ix=(x_train==b)
        base[b]=(float(e_train[ix].sum())+1)/(int(ix.sum())+2)
    if z_train is None:
        return np.array([base[int(b)] for b in x_eval])
    group={}
    for b in (0,1):
        for z in np.unique(z_train):
            ix=(x_train==b)&(z_train==z)
            group[b,int(z)]=(float(e_train[ix].sum()) + lam*base[b])/(int(ix.sum())+lam)
    return np.array([group.get((int(b),int(z)),base[int(b)]) for b,z in zip(x_eval,z_eval)])


def calibrate_prefix_randomized(score, err, alpha):
    # Threshold rules with randomized inclusion at boundary.
    # maximize empirical coverage over risk <=alpha, based solely on calibration.
    levels=np.unique(score)
    n=np.array([(score==v).sum() for v in levels],dtype=float)
    m=np.array([err[score==v].sum() for v in levels],dtype=float)
    best=np.zeros(len(levels)); bestcov=0.0
    accum_n=0.0; accum_m=0.0
    for j in range(len(levels)):
        a=np.zeros(len(levels)); a[:j]=1.0
        q=m[j]-alpha*n[j]
        budget=accum_m-alpha*accum_n
        if budget > 1e-10:
            # Earlier whole groups are infeasible; no eligible prefix at this j.
            f=0.0
            feasible=False
            if q < -1e-12:
                f_min=budget/(-q)
                feasible=(f_min<=1+1e-12)
                f=1.0 if feasible else 0.0
        elif q > 1e-12:
            f=min(1.0,max(0.0,-budget/q))
            feasible=True
        else:
            f=1.0
            feasible=True
        if feasible:
            a[j]=f
            cov=(accum_n+f*n[j])/sum(n)
            if cov>bestcov+1e-12:
                best=a;bestcov=cov
        accum_n+=n[j];accum_m+=m[j]
    # Numerical defensive check.
    if np.dot(m-alpha*n,best)>1e-6:
        raise ValueError('Calibration risk constraint violation')
    return dict(zip(levels,best))


def test_risk_coverage(score, err, coverages=(0.25,0.4,0.5,0.6,0.75)):
    """Held-out risk at exact coverage, ties randomized in expectation.

    This uses training-fitted score order but does not use test labels to
    select an operating point. It is a *diagnostic curve*, not deployment.
    """
    levels=np.unique(score)
    counts=np.array([(score==v).sum() for v in levels],dtype=float)
    errors=np.array([err[score==v].sum() for v in levels],dtype=float)
    cum_n=np.r_[0,np.cumsum(counts)]
    cum_e=np.r_[0,np.cumsum(errors)]
    risk=[]
    for c in coverages:
        count_target=c*len(err)
        j=int(np.searchsorted(cum_n,count_target,side='left')-1)
        j=max(0,min(j,len(levels)-1))
        frac=(count_target-cum_n[j])/counts[j]
        risk.append(float((cum_e[j]+frac*errors[j])/count_target))
    return list(risk)


def evaluate_policy(score, err, probs):
    # Numeric scores created from training tables; floating keys match exactly.
    a=np.array([probs.get(s,0.0) for s in score],dtype=float)
    cover=float(a.mean())
    risk=float(np.dot(a,err)/a.sum()) if a.sum()>0 else None
    return cover,risk


def one_seed(seed, alpha, n_train, n_cal, n_test, ks, qs):
    rng=np.random.default_rng(seed)
    xt,et=draw(rng,n_train)
    xc,ec=draw(rng,n_cal)
    xv,ev=draw(rng,n_test)
    out={}
    out['baseline_oracle']=oracle_coverage(xt,et,alpha)
    tr=scores_from_training(xt,et,xt)
    ca=scores_from_training(xt,et,xc)
    te=scores_from_training(xt,et,xv)
    policy=calibrate_prefix_randomized(ca,ec,alpha)
    out['baseline_test']=evaluate_policy(te,ev,policy)
    out['baseline_curve']=test_risk_coverage(te,ev)
    for k in ks:
        zt=nuisance(rng,n_train,k);zc=nuisance(rng,n_cal,k);zv=nuisance(rng,n_test,k)
        out[f'null_{k}_oracle']=oracle_coverage(xt,et,alpha,zt)
        ca=scores_from_training(xt,et,xc,zt,zc)
        te=scores_from_training(xt,et,xv,zt,zv)
        policy=calibrate_prefix_randomized(ca,ec,alpha)
        out[f'null_{k}_test']=evaluate_policy(te,ev,policy)
        out[f'null_{k}_curve']=test_risk_coverage(te,ev)
    for q in qs:
        st=positive_proxy(rng,et,q);sc=positive_proxy(rng,ec,q);sv=positive_proxy(rng,ev,q)
        out[f'positive_{q}_oracle']=oracle_coverage(xt,et,alpha,st)
        ca=scores_from_training(xt,et,xc,st,sc)
        te=scores_from_training(xt,et,xv,st,sv)
        policy=calibrate_prefix_randomized(ca,ec,alpha)
        out[f'positive_{q}_test']=evaluate_policy(te,ev,policy)
        out[f'positive_{q}_curve']=test_risk_coverage(te,ev)
    return out


def ci_mean(d, rng, resamples=1500):
    d=np.asarray(d,float)
    idx=rng.integers(0,len(d),size=(resamples,len(d)))
    return [float(v) for v in np.quantile(d[idx].mean(axis=1),[0.025,0.975])]


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--seeds',type=int,default=200)
    p.add_argument('--alpha',type=float,default=0.05)
    p.add_argument('--n-train',type=int,default=500)
    p.add_argument('--n-cal',type=int,default=500)
    p.add_argument('--n-test',type=int,default=5000)
    p.add_argument('--output',type=Path,default=Path('/mnt/data/oracle_negative_control_results.json'))
    args=p.parse_args()
    ks=[2,4,16];qs=[0.90,0.98]
    runs=[one_seed(20261104+i,args.alpha,args.n_train,args.n_cal,args.n_test,ks,qs) for i in range(args.seeds)]
    bootstrap_rng=np.random.default_rng(4815)
    results={}
    bo=np.array([r['baseline_oracle'] for r in runs]); bc=np.array([r['baseline_test'][0] for r in runs]);
    br=np.array([r['baseline_test'][1] if r['baseline_test'][1] is not None else np.nan for r in runs]);
    results['baseline']={'oracle_train_coverage_mean':float(bo.mean()),'test_coverage_mean':float(bc.mean()),'test_risk_mean':float(np.nanmean(br)),'test_risk_violation_rate':float(np.mean(br>args.alpha)),
                         'test_risk_at_cover50_mean':float(np.mean([r['baseline_curve'][2] for r in runs])),
                         'test_curve_mean_risk':float(np.mean([np.mean(r['baseline_curve']) for r in runs]))}
    for kind,values in [('null',ks),('positive',qs)]:
        for z in values:
            key=f'{kind}_{z}'
            ore=np.array([r[f'{key}_oracle'] for r in runs]); tc=np.array([r[f'{key}_test'][0] for r in runs]);
            risks=np.array([r[f'{key}_test'][1] if r[f'{key}_test'][1] is not None else np.nan for r in runs]);
            curve_base=np.array([r['baseline_curve'] for r in runs])
            curve_this=np.array([r[f'{key}_curve'] for r in runs])
            risk50_delta=curve_this[:,2]-curve_base[:,2]
            curve_delta=curve_this.mean(axis=1)-curve_base.mean(axis=1)
            results[key]={
                'oracle_train_coverage_mean':float(ore.mean()),
                'oracle_minus_baseline_mean_pp':float(100*(ore-bo).mean()),
                'oracle_minus_baseline_95ci_pp':[100*x for x in ci_mean(ore-bo,bootstrap_rng)],
                'test_coverage_mean':float(tc.mean()),
                'test_minus_baseline_mean_pp':float(100*(tc-bc).mean()),
                'test_minus_baseline_95ci_pp':[100*x for x in ci_mean(tc-bc,bootstrap_rng)],
                'test_risk_mean':float(np.nanmean(risks)),
                'test_risk_violation_rate':float(np.mean(risks>args.alpha)),
                'test_feasible_rate':float(np.mean((risks<=args.alpha)&(tc>0))),
                'test_risk_at_cover50_minus_baseline_mean_pp':float(100*risk50_delta.mean()),
                'test_risk_at_cover50_minus_baseline_95ci_pp':[100*x for x in ci_mean(risk50_delta,bootstrap_rng)],
                'test_curve_mean_risk_minus_baseline_pp':float(100*curve_delta.mean()),
                'test_curve_mean_risk_minus_baseline_95ci_pp':[100*x for x in ci_mean(curve_delta,bootstrap_rng)],
            }
    out={'study':'Synthetic paired negative/positive controls for oracle inflation',
         'design':{'seeds':args.seeds,'alpha':args.alpha,'train':args.n_train,'calibration_n':args.n_cal,'test':args.n_test,
         'base':'X=0 with p=0.4, error=0.01; X=1 with p=0.6, error=0.35',
         'negative':'Z uniform independent of (X,error) for k=2/4/16',
         'positive':'S with P(S=error)=q, q=0.90/0.98, synthetic hypothetical observed proxy',
         'calibration':'empirical risk constraint, randomized boundary, not guaranteed population risk'},
         'results':results}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(out,ensure_ascii=False,indent=2))
    print(json.dumps(out,ensure_ascii=False,indent=2))

if __name__=='__main__':
    main()

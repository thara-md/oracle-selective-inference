#!/usr/bin/env python3
"""Independent numerical LP spot check and Figure 3 from P3 synthetic CSV."""
from pathlib import Path
from fractions import Fraction as F
from collections import defaultdict
import csv, json, random, sys
import numpy as np
from scipy.optimize import linprog
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from P3_minimal_evidence_audit import optimistic_bound, synthetic_blocks, exact_lp, BASE_SEED, N, ERRORS

P=Path(__file__).resolve().parent
r=[]
with (P/'P3_synthetic_information_infeasible_targets.csv').open(encoding='utf-8-sig') as f:r=list(csv.DictReader(f))
r5=[x for x in r if x['alpha']=='1/20']
# Unit tests: all true feasibility states and strict boundaries
assert len(r5)==2167
assert sum(int(z['min_blocks'])==1 for z in r5)==479
assert sum(int(z['largest_first_50']) for z in r5)==840
assert sum(int(z['largest_first_75']) for z in r5)==1596

rng=random.Random(20261105)
tests=0
for i in range(100):
    scenario=['concentrated_mixed','dispersed_mixed','homogeneous','high_uniqueness'][i%4]
    seed=rng.randrange(100)
    j=['concentrated_mixed','dispersed_mixed','homogeneous','high_uniqueness'].index(scenario)
    blocks=synthetic_blocks(scenario,np.random.default_rng(np.random.SeedSequence([BASE_SEED,j,seed])))
    verified=[(n,e) for k,(n,e) in enumerate(blocks) if rng.randrange(2)]
    alpha=[F(1,50),F(1,20),F(1,10)][i%3]
    exact=optimistic_bound(verified,N,ERRORS,alpha)
    nr=N-sum(n for n,e in verified)
    er=ERRORS-sum(e for n,e in verified)
    opt_blocks=list(verified)
    if nr-er:opt_blocks.append((nr-er,0))
    if er:opt_blocks.append((er,er))
    costs=np.array([float(e-alpha*n) for n,e in opt_blocks])
    ns=np.array([n for n,e in opt_blocks],dtype=float)
    result=linprog(-ns, A_ub=costs.reshape(1,-1), b_ub=np.array([0.]), bounds=[(0.,1.)]*len(ns),method='highs')
    assert result.success
    numeric=(-result.fun)/N
    assert abs(float(exact)-numeric)<1e-9,(exact,numeric)
    tests+=1

colors=['#31547d','#b56a35','#4f876a']
scenarios=['concentrated_mixed','dispersed_mixed','high_uniqueness']
labels=['Concentrated mixed','Dispersed mixed','High uniqueness']
maxbs=[3,7,6]
fig,axs=plt.subplots(1,2,figsize=(11.5,4.4),gridspec_kw={'width_ratios':[1.22,1]})
ax=axs[0]
for s,label,color,maxb in zip(scenarios,labels,colors,maxbs):
    vals=[z for z in r5 if z['scenario']==s]
    vals_k=[100*sum(int(v['largest_first_blocks'])<=k for v in vals)/len(vals) for k in range(maxb+1)]
    ax.plot(range(maxb+1),vals_k,marker='o',lw=2,ms=5,label=f'{label} (m={maxb})',color=color)
ax.set_xlim(-0.15,7.15);ax.set_ylim(0,104)
ax.set_xticks(range(8));ax.set_yticks(range(0,101,20));ax.set_xlabel('Verified collision blocks (largest first)')
ax.set_ylabel('Infeasible targets excluded (%)');ax.grid(alpha=.24)
ax.legend(loc='lower right',fontsize=8,frameon=False)
ax.set_title('A. Certificate detection as evidence grows',weight='bold',loc='left',fontsize=11)
ax.text(0.0,-0.28,'Conditional on targets: case-wise feasible,\nfull-representation infeasible (risk 5%).',transform=ax.transAxes,fontsize=8,va='top')

ax=axs[1]
from collections import Counter
for i,(s,label,color,maxb) in enumerate(zip(scenarios,labels,colors,maxbs)):
    vals=[int(z['min_blocks']) for z in r5 if z['scenario']==s]
    counter=Counter(vals)
    for k,c in counter.items():
        ax.barh(i+.34*(k/(maxb+1)-.5),100*c/len(vals),height=.30/(maxb+1),color=color,alpha=.8)
# simpler summary as median and 1-block proof rate, use overlay table text
ax.clear()
heights=[np.median([int(z['min_blocks']) for z in r5 if z['scenario']==s]) for s in scenarios]
ax.barh(range(3),heights,color=colors,height=.55)
ax.invert_yaxis();ax.set_yticks(range(3),labels=labels);ax.set_xlim(0,7)
ax.set_xticks(range(8));ax.set_xlabel('Minimum collision blocks (median)')
ax.grid(axis='x',alpha=.24)
for j,(h,s,maxb) in enumerate(zip(heights,scenarios,maxbs)):
    vals=[z for z in r5 if z['scenario']==s]
    one=100*sum(int(v['min_blocks'])==1 for v in vals)/len(vals)
    ax.text(h+.12,j,f'{h:g} / {maxb};  1-block: {one:.1f}%',va='center',fontsize=8)
ax.set_title('B. Best possible proof sparsity',weight='bold',loc='left',fontsize=11)
ax.text(0.0,-0.28,'Exhaustive search over all subsets (m≤7).\nHomogeneous blocks have no information-specific gap.',transform=ax.transAxes,fontsize=8,va='top')
fig.suptitle('How much verified input equivalence is enough?',fontsize=14,fontweight='bold',y=1.03)
fig.subplots_adjust(left=.08,right=.99,bottom=.30,top=.83,wspace=.45)
for ext in ['png','pdf']:
    fig.savefig(P/f'Figure3_Minimal_Evidence_P3.{ext}',dpi=220,bbox_inches='tight')
plt.close(fig)
checks={'scipy_lp_independent_checked':tests,'max_numeric_difference_tolerance':1e-9,
 'info_infeasible_alpha5_pairs':len(r5),'minimum_single_block':479,
 'size50_excluded':840,'size75_excluded':1596,
 'synthetic_data_source':'P3_minimal_evidence_audit.py: new independent 400 datasets, not legacy 400'}
(P/'P3_independent_validation.json').write_text(json.dumps(checks,indent=2)+'\n',encoding='utf-8')
print(checks)

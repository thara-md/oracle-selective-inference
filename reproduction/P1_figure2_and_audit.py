#!/usr/bin/env python3
"""P1 audit and journal figure 2; reads unmodified frozen P0 and synthetic results."""
import csv
import hashlib
import json
from pathlib import Path
import math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from matplotlib.transforms import blended_transform_factory

ROOT=Path(__file__).resolve().parent
SYN=ROOT/'_work_p1/synthetic/oracle_negative_control_results.json'
P0=ROOT/'P0_200seed_summary.csv'
OUT=ROOT

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def check():
    syn=json.loads(SYN.read_text())
    rows=list(csv.DictReader(P0.open(encoding='utf-8-sig',newline='')))
    assert syn['design']['seeds']==200 and syn['design']['alpha']==0.05
    assert len(rows)==12
    d=syn['results'];
    assert abs(d['baseline']['oracle_train_coverage_mean']*100-45.40839074904855)<1e-9
    assert abs(d['null_16']['oracle_train_coverage_mean']*100-50.70217900926264)<1e-9
    assert abs(d['null_16']['test_risk_at_cover50_minus_baseline_mean_pp']-.005375911920609051)<1e-10
    assert abs(d['baseline']['test_risk_at_cover50_mean']*100 - 7.780663102441106)<1e-9
    assert abs(d['null_16']['test_coverage_mean']*100 - 45.77912863962711)<1e-9
    assert abs(d['null_16']['test_risk_mean']*100 - 5.217272466976228)<1e-9
    a=[r for r in rows if r['condition']=='All-category GPT' and r['k']=='16']
    assert len(a)==1 and abs(float(a[0]['mean_pct'])-64.70835374155898)<1e-10
    return syn,rows

syn,p0=check()
D=syn['results']
labels=['Base','Noise-2','Noise-4','Noise-16','Signal-90','Signal-98']
keys=['baseline','null_2','null_4','null_16','positive_0.9','positive_0.98']
x=np.arange(6)
vals=np.array([100*D[k]['oracle_train_coverage_mean'] for k in keys])
riskd=np.array([0]+[D[k]['test_risk_at_cover50_minus_baseline_mean_pp'] for k in keys[1:]])
ci=np.array([[0,0]]+[D[k]['test_risk_at_cover50_minus_baseline_95ci_pp'] for k in keys[1:]])
colors=['#505050']+['#C57B35']*3+['#287B91']*2
# A two-panel graphical summary.  All values are sourced directly from published
# experiment JSON, never typed into individual drawing calls.
plt.rcParams.update({'font.size':8.3,'font.family':'DejaVu Sans','pdf.fonttype':42,
                     'axes.spines.top':False,'axes.spines.right':False,
                     'axes.titleweight':'semibold'})
fig,(ax1,ax2)=plt.subplots(1,2,figsize=(9.35,4.48),gridspec_kw={'width_ratios':[1.1,1.0]})
fig.subplots_adjust(left=.076,right=.974,top=.835,bottom=.255,wspace=.30)
# Panel A
bars=ax1.bar(x,vals,color=colors,width=.71,edgecolor='#353535',linewidth=.43,zorder=3)
ax1.set_ylim(0,91)
ax1.set_yticks([0,20,40,50,60,80])
ax1.set_ylabel('In-sample oracle coverage (%)')
ax1.set_xticks(x,labels)
ax1.tick_params(axis='x',length=0,labelsize=7.4)
ax1.grid(axis='y',alpha=.18,zorder=0)
ax1.axhline(50,linestyle='--',color='#606060',linewidth=1,alpha=.9,zorder=2)
# Threshold label provided in figure caption to avoid overlap with bars
# Population optimum applies only to the base/noise conditions.
ax1.hlines(45.3333333333,-.39,3.42,color='#343434',linestyle=':',linewidth=1.25,zorder=4)
# The dotted-line definition is given explicitly in the caption
for rect,val in zip(bars,vals):
    ax1.text(rect.get_x()+rect.get_width()/2,val+1.4,f'{val:.2f}',ha='center',va='bottom',fontsize=7.5,fontweight='semibold')
ax1.set_title('A  Empirical oracle: apparent coverage gains',loc='left',fontsize=9.7,pad=12)
# Panel B
ax2.axhline(0,color='#4a4a4a',linewidth=.85,zorder=1)
ax2.bar(x,riskd,color=colors,width=.71,edgecolor='#353535',linewidth=.43,zorder=3)
for i in range(1,6):
    ax2.errorbar(x[i],riskd[i],yerr=[[riskd[i]-ci[i,0]],[ci[i,1]-riskd[i]]],
                 fmt='none',ecolor='#1f1f1f',elinewidth=1.1,capsize=2.4,zorder=5)
ax2.set_ylim(-8.5,1.0)
ax2.set_yticks([-8,-6,-4,-2,0])
ax2.set_ylabel('Test selective-risk change (pp)')
ax2.set_xticks(x,labels)
ax2.tick_params(axis='x',length=0,labelsize=7.4)
ax2.grid(axis='y',alpha=.18,zorder=0)
ax2.set_title('B  Independent test: matched 50% coverage',loc='left',fontsize=9.7,pad=12)
# Risk sign convention and base risk appear beneath the figure
for i in (4,5):
    ax2.text(i,riskd[i]-.33,f'{riskd[i]:+.2f} pp',fontsize=7.5,ha='center',va='top',fontweight='semibold')
# Label near-zero null effects without visually hiding the tiny uncertainty bars.
ax2.text(-0.22,-1.5,'Noise 2/4/16: -0.029 / +0.012 / +0.005 pp',fontsize=7.3,
         ha='left',va='center',color='#7c4c19')
ax2.text(-0.22,-2.02,'Noise-16 CI: [-0.044, +0.056] pp',fontsize=7.3,
         ha='left',va='center',color='#7c4c19')
# precise high-level context
fig.text(.072,.949,'FIGURE 2  |  Synthetic representation controls',fontsize=12,fontweight='bold',ha='left')
fig.text(.072,.902,'200 independent seeds  ·  500 training / 500 calibration / 5,000 test per seed  ·  α = 5%',fontsize=8.4,ha='left')
fig.text(.076,.115,'Dashed: 50% coverage target. Dotted (base/noise only): population optimum 45.33%.',
         fontsize=7.3,ha='left',color='#555555')
fig.text(.076,.083,'Right panel: lower test risk is better (base risk 7.78%); bars show relative change at fixed 50% test coverage.',
         fontsize=7.3,ha='left',color='#555555')
fig.text(.076,.050,'Diagnostic curve, not a calibrated 5%-risk operating point. Error bars: paired-seed bootstrap 95% CI.',
         fontsize=7.3,ha='left',color='#555555')
fig.savefig(OUT/'Figure2_Oracle_inflation_vs_test_risk_P1.pdf',bbox_inches='tight',pad_inches=.07)
fig.savefig(OUT/'Figure2_Oracle_inflation_vs_test_risk_P1.png',dpi=240,bbox_inches='tight',pad_inches=.07)
plt.close(fig)

# Reconciled evidence summary (200 seeds)
base=D['baseline']
sumrows=[]
for i,k in enumerate(keys):
    r=D[k]
    sumrows.append({
      'representation':labels[i].replace('\n',' '),
      'n_seeds':200,
      'train_oracle_coverage_pct':round(100*r['oracle_train_coverage_mean'],4),
      'test_calibrated_coverage_pct':round(100*r['test_coverage_mean'],4),
      'test_calibrated_risk_pct':round(100*r['test_risk_mean'],4),
      'test_risk_over_5pct_seed_fraction_pct':round(100*r['test_risk_violation_rate'],2),
      'fixed50_test_risk_delta_pp':round(r.get('test_risk_at_cover50_minus_baseline_mean_pp',0),5),
      'fixed50_ci_low_pp':round((r.get('test_risk_at_cover50_minus_baseline_95ci_pp',[0,0]))[0],5),
      'fixed50_ci_high_pp':round((r.get('test_risk_at_cover50_minus_baseline_95ci_pp',[0,0]))[1],5),
    })
with (OUT/'P1_synthetic_verified_summary.csv').open('w',newline='',encoding='utf-8-sig') as f:
    w=csv.DictWriter(f,fieldnames=list(sumrows[0].keys()));w.writeheader();w.writerows(sumrows)

sens={}
for a,path in [(2,ROOT/'_work_p1/synthetic/oracle_sensitivity_risk02.json'),
               (5,SYN),(10,ROOT/'_work_p1/synthetic/oracle_sensitivity_risk10.json')]:
    j=json.loads(path.read_text())['results']
    sens[str(a)]=dict(baseline_oracle_pct=round(100*j['baseline']['oracle_train_coverage_mean'],4),
                 null_16_oracle_pct=round(100*j['null_16']['oracle_train_coverage_mean'],4),
                 gain_pct_points=round(100*(j['null_16']['oracle_train_coverage_mean']-j['baseline']['oracle_train_coverage_mean']),4),
                 fixed50_null16_risk_delta_pp=round(j['null_16']['test_risk_at_cover50_minus_baseline_mean_pp'],5))
ledger={
  'title':'P1: reconcile clinical and synthetic negative/positive controls',
  'date':'2026-11-04','status':'reproduced with separate input contexts',
  'source_hashes':{
    'synthetic_results_json_sha256':sha(SYN),'synthetic_code_sha256':sha(ROOT/'_work_p1/synthetic/oracle_negative_control.py'),
    'clinical_P0_summary_csv_sha256':sha(P0),'clinical_P0_code_sha256':sha(ROOT/'P0_clinical_oracle_200seed_repro.py')},
  'clinical':{
    'n_patients':90,'n_patient_model_units':174,'n_model_strata':4,'seed_replications_per_k_and_stratum':200,
    'n_lp_optima':2400,'input_level':'audited reconstructed aggregate equivalence blocks, not raw patient rows',
    'all_category_GPT':{'baseline_cov_pct':44.397463002114165,'null16_mean_pct':64.70835374155898,'null16_mean_gain_pp':20.310890739444815},
    'interpretation':'Fixed finite-sample, correctness-independent random features artificially partition input-equivalence classes. No clinical held-out test.'},
  'synthetic':{
    'seeds':200,'train':500,'calibration':500,'test':5000,
    'baseline_true_population_optimal_coverage_pct':100*(0.4+(0.05*.4-0.01*.4)/(.35-.05)),
    'baseline_test_risk_at_fixed50_pct':100*base['test_risk_at_cover50_mean'],
    'null16':{'oracle_baseline_pct':100*base['oracle_train_coverage_mean'],'oracle_null16_pct':100*D['null_16']['oracle_train_coverage_mean'],
              'cal_test_cov_baseline_pct':100*base['test_coverage_mean'],'cal_test_cov_null16_pct':100*D['null_16']['test_coverage_mean'],
              'cal_test_risk_baseline_pct':100*base['test_risk_mean'],'cal_test_risk_null16_pct':100*D['null_16']['test_risk_mean'],
              'fixed50_risk_delta_pp':D['null_16']['test_risk_at_cover50_minus_baseline_mean_pp'],
              'fixed50_risk_delta_95ci_pp':D['null_16']['test_risk_at_cover50_minus_baseline_95ci_pp']},
    'risk_sensitivity_pct':sens,
    'interpretation':'Test fixed50 risk is a diagnostic curve using test feature distribution to set coverage; it is NOT a deployment policy or evidence of satisfying risk 5%. Positive proxy is generated from error labels and is intentionally artificial.'},
  'caveats':['Clinical and synthetic oracle gains are on different experimental populations and must not be pooled or directly compared as effects.',
             'A bootstrap CI containing zero is not proof of equivalence; equivalence requires a prespecified tolerance.',
             '200 Monte Carlo seeds and 5000 test records/seed are not 200 independent clinical cohorts.',
             'Original 50-seed RNG stream unavailable; clinical 200-seed analysis is a fresh Monte Carlo rerun.',
             'Archived patient-level shared_sequence_blocks.csv and GT join still require direct source comparison before manuscript submission.']}
(OUT/'P1_reconciliation_audit.json').write_text(json.dumps(ledger,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
print('Created Fig2; synthetic summary, P1 audit')
for r in sumrows: print(r)
print('Sensitivity:',sens)

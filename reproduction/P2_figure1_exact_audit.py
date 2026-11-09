#!/usr/bin/env python3
"""P2 Figure 1 reconstruction from verified aggregate-equivalence blocks.

Scope: archived 2026-09-15 derived audit counts, independently reconstructed
in P0 2026-11-04; NOT patient-level rerun. Only the single documented best
low-H threshold point is plotted; a full threshold scatter requires per-unit H.

A selector assigns a common acceptance probability to identical (S,modal)
blocks. All oracle quantities are finite-sample, GT-informed expected policies,
not out-of-sample statistical guarantees.
"""
import csv, hashlib, json
from fractions import Fraction as F
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

ROOT=Path(__file__).resolve().parent
P0=ROOT/'P0_reconstructed_block_inputs.csv'
ALPHA=F(1,20)
N=86
# Distinct sequence blocks: shared exact labels must share q.
BLOCKS=[('PDAC×20',32,6),('NET×20',2,1),('19PDAC+1NET',2,0)] + [(f'correct singleton{i}',1,0) for i in range(26)] + [(f'error singleton{i}',1,1) for i in range(24)]
CASE=[('correct',55,0),('error',31,31)]
BEST_LOW_H=(F(34,86), F(7,34))  # historical per-unit H audit: one verified point


def check_inputs():
    rows=list(csv.DictReader(P0.open('r',encoding='utf-8-sig',newline='')))
    expected={r['block_kind']:(int(r['n']),int(r['errors'])) for r in rows if r['condition']=='All-category GPT'}
    assert expected=={
        'PDAC×20':(32,6), 'NET×20':(2,1),
        '19PDAC+1NET shared sequence':(2,0),
        'unique correct sequences (individual singletons)':(26,0),
        'unique erroneous sequences (individual singletons)':(24,24),
    },expected
    assert sum(n for _,n,e in BLOCKS)==86
    assert sum(e for _,n,e in BLOCKS)==31
    return hashlib.sha256(P0.read_bytes()).hexdigest()

def min_risk_exact(blocks, coverage):
    c=F(coverage)
    if not 0 < c <= 1: return None
    mass=c*N; remaining=mass; err=F(0)
    for label,n,e in sorted(blocks,key=lambda z:(F(z[2],z[1]),z[0])):
        take=min(F(n),remaining)
        if take>0:
            err+=take*F(e,n)
            remaining-=take
        if remaining==0:break
    assert remaining==0
    return err/mass

def maximum_coverage_exact(blocks,alpha):
    ordered=sorted(blocks,key=lambda z:(F(z[2],z[1]),z[0]))
    accept=F(0); err=F(0); lamb=F(0)
    for name,n,e in ordered:
        if err + e <= alpha*(accept+n):
            accept+=n; err+=e
        else:
            budget=alpha*accept-err
            assert budget>=0
            # If n>0 all blocks positive rate relative to alpha here.
            part=budget/(F(e,n)-alpha)
            assert 0<=part<=n
            accept+=part;err+=part*F(e,n)
            lamb=F(n)/(F(e)-alpha*n)
            break
    primal=accept/N
    dual=sum((max(F(0),F(n)-lamb*(F(e)-alpha*n)) for _,n,e in blocks),F(0))/N
    assert primal==dual and err<=alpha*accept
    return {'fraction':primal,'dual_fraction':dual,'lambda':lamb}

def vector_risk(blocks,coverage):
    ordered=sorted(blocks,key=lambda z:(z[2]/z[1],z[0]))
    bn=np.asarray([v[1] for v in ordered],dtype=float)
    be=np.asarray([v[2] for v in ordered],dtype=float)
    m=N*np.asarray(coverage)
    e=np.zeros_like(m)
    for n,er in zip(bn,be):
        take=np.minimum(n,np.maximum(0,m-(n if False else 0))) # placeholder overwritten below
    remaining=m.copy();e[:]=0
    for n,er in zip(bn,be):
        take=np.minimum(n,np.maximum(0,remaining))
        e+=take*(er/n)
        remaining-=take
    assert np.all(np.abs(remaining)<1e-7)
    return e/m

def main():
    data_hash=check_inputs()
    seq=maximum_coverage_exact(BLOCKS,ALPHA)
    case=maximum_coverage_exact(CASE,ALPHA)
    c_t=F(1,2)
    r_seq=min_risk_exact(BLOCKS,c_t)
    r_case=min_risk_exact(CASE,c_t)
    assert seq['fraction']==F(210,473),seq
    assert seq['lambda']==F(80,11)
    assert case['fraction']==F(550,817),case
    assert r_seq==F(45,688),(r_seq,F(45,688))
    assert r_case==0
    # All non-unanimous repeated blocks are outcome homogeneous; exact relaxation
    # has the same risk–coverage frontier after splitting them into individual rows.
    relaxed=[('PDAC×20',32,6),('NET×20',2,1)] + [(f'correct other {i}',1,0) for i in range(28)] + [(f'incorrect other {i}',1,1) for i in range(24)]
    pdac_only=[('PDAC×20',32,6)]+[(f'correct outside {i}',1,0) for i in range(29)]+[(f'error outside {i}',1,1) for i in range(25)]
    pdac_bound=maximum_coverage_exact(pdac_only,ALPHA)
    assert pdac_bound['fraction']==F(435,946)
    assert maximum_coverage_exact(relaxed,ALPHA)['fraction']==seq['fraction']
    x=np.unique(np.r_[np.linspace(0.00001,1,2500),
                       np.array([float(F(n,N)) for n in (28,55,60,62,86)]),
                       float(seq['fraction']),float(case['fraction']),0.5])
    rs=vector_risk(BLOCKS,x);rc=vector_risk(CASE,x);rr=vector_risk(relaxed,x)
    assert np.max(abs(rs-rr))<1e-12
    assert np.min(rs-rc)>-1e-12
    idx=np.argmin(abs(x-0.5))
    assert abs(rs[idx]-float(r_seq))<1e-10
    assert abs(rs[-1]-31/86)<1e-11
    # Write exact audit summary with independent mathematical checkpoints.
    summary={
       'source':'P0 2026-11-04 aggregate-equivalence reconstruction from archival 2026-09-15 audit',
       'input_sha256':data_hash,'scope':'All-category GPT; n=86; 31 errors; fixed original predictions',
       'oracle_representation':'ordered 20 labels + saved modal; equal input => equal adoption probability',
       'alpha_fraction':'1/20','target_coverage_fraction':'1/2',
       'sequence_oracle_coverage_fraction':str(seq['fraction']),
       'casewise_oracle_coverage_fraction':str(case['fraction']),
       'primal_dual_multiplier':str(seq['lambda']),
       'minimum_risk_at_50pct_coverage_fraction':str(r_seq),
       'casewise_minimum_risk_at_50pct_coverage_fraction':str(r_case),
       'single_pdac_collision_coverage_upper_bound_fraction':str(pdac_bound['fraction']),
       'best_low_H_point':{'accepted':34,'errors':7,'coverage_fraction':'17/43','risk_fraction':'7/34'},
       'max_curve_gap_full_vs_unanimity_relaxation_float':float(np.max(abs(rs-rr))),
       'clinical_validation':'NONE; no test-set risk guarantee',
       'full_threshold_scatter':'NOT available without individual entropy/GT; single archived verified best point shown'
    }
    (ROOT/'P2_Figure1_math_audit.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    with (ROOT/'P2_Figure1_curve_data.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.writer(f);w.writerow(['coverage_pct','sequence_oracle_minrisk_pct','casewise_oracle_minrisk_pct','unanimity_relaxation_minrisk_pct'])
        for xx,yy,zz,uu in zip(x,rs,rc,rr):w.writerow([f'{100*xx:.8f}',f'{100*yy:.8f}',f'{100*zz:.8f}',f'{100*uu:.8f}'])
    font={'font.family':'DejaVu Sans','pdf.fonttype':42,'ps.fonttype':42,'axes.spines.top':False,'axes.spines.right':False}
    with plt.rc_context(font):
        fig,ax=plt.subplots(figsize=(10.3,6.0),dpi=140)
        fig.subplots_adjust(left=.091,right=.983,top=.82,bottom=.23)
        X=x*100;RS=rs*100;RC=rc*100
        ax.fill_between(X,0,RC,color='#e1e4e8',label='_nolegend_',zorder=0)
        ax.fill_between(X,RC,RS,facecolor='#f7ecec',edgecolor='#b75e63',hatch='///',lw=0, label='_nolegend_',zorder=1)
        line1,=ax.plot(X,RS,color='#9f1722',lw=2.4,label='Full-sequence oracle',zorder=4)
        line2,=ax.plot(X,RC,color='#215184',lw=2.1,ls=(0,(5,2)),label='Case-wise oracle',zorder=4)
        best=BEST_LOW_H
        pt=ax.scatter([float(best[0])*100],[float(best[1])*100],s=55,facecolors='white',edgecolors='#2b7651',linewidths=1.9,zorder=8,label='Best observed low-H threshold')
        ax.plot([0,100],[5,5],color='#666666',lw=1,ls=':',zorder=2)
        ax.scatter([float(seq['fraction'])*100,float(case['fraction'])*100],[5,5],s=48,c=['#9f1722','#215184'],zorder=10)
        ax.scatter([50],[5],marker='X',s=115,color='#202124',zorder=11)
        ax.scatter([50],[float(r_seq)*100],marker='s',s=40,color='#9f1722',zorder=10)
        ax.annotate('Minimum risk at 50%: 6.54%',xy=(50,float(r_seq)*100),xytext=(57,13.0),
                    fontsize=10,color='#9f1722',fontweight='bold',
                    arrowprops=dict(arrowstyle='-',color='#9f1722',lw=1.1),
                    bbox=dict(facecolor='white',edgecolor='none',alpha=.85,pad=1.5))
        ax.annotate('Target (50% coverage, 5% risk)',xy=(50,5),xytext=(56,1.3),
                    fontsize=9.5,color='#26282d',
                    arrowprops=dict(arrowstyle='-',color='#444444',lw=.8),
                    bbox=dict(facecolor='white',edgecolor='none',alpha=.85,pad=1.2))
        ax.annotate('44.40%',xy=(float(seq['fraction'])*100,5),xytext=(30.1,9.1),
                    fontsize=10,color='#9f1722',weight='bold',
                    arrowprops=dict(arrowstyle='-',color='#9f1722',lw=.9))
        ax.annotate('67.32%',xy=(float(case['fraction'])*100,5),xytext=(72,9.1),
                    fontsize=10,color='#215184',weight='bold',
                    arrowprops=dict(arrowstyle='-',color='#215184',lw=.9))
        ax.annotate('Low-H: 34 accepted / 7 errors\n(risk 20.59%)',xy=(float(best[0])*100,float(best[1])*100),xytext=(48,27.0),
                    fontsize=9.6,color='#2b7651',
                    arrowprops=dict(arrowstyle='-',color='#2b7651',lw=.9),
                    bbox=dict(facecolor='white',edgecolor='none',alpha=.9,pad=1.5))
        ax.set(xlim=(0,100),ylim=(0,39),xlabel='Expected coverage (%)',ylabel='Expected selective risk (%)')
        ax.set_xticks(np.arange(0,101,10));ax.set_yticks(np.arange(0,40,5))
        ax.grid(color='#cdd0d7',ls=':',lw=.65,zorder=0)
        ax.set_axisbelow(True)
        fig.suptitle('Empirical risk–coverage limits under fixed label information',y=.963,fontsize=15.0,weight='bold')
        fig.text(.5,.908,'All-category GPT   |   n = 86   |   T = 20 repeated labels   |   ground-truth-informed in-sample oracles',ha='center',fontsize=10.1,color='#444444')
        fig.legend(handles=[line1,line2,pt],loc='upper center',bbox_to_anchor=(.5,.869),ncol=3,frameon=False,fontsize=9.6)
        patches=[Patch(facecolor='#f7ecec',hatch='///',edgecolor='#b75e63',label='Additionally excluded by label-information restriction'),
                 Patch(facecolor='#e1e4e8',edgecolor='#e1e4e8',label='Excluded even with case-wise identification')]
        fig.legend(handles=patches,loc='lower center',bbox_to_anchor=(.5,.113),ncol=2,frameon=False,fontsize=9.1)
        fig.text(.5,.052,'Risk at zero coverage is undefined. Curves allow randomized acceptance; all bounds condition on this fixed clinical cohort.',ha='center',fontsize=8.7,color='#555555')
        fig.savefig(ROOT/'Figure1_Oracle_risk_coverage_P2.pdf',bbox_inches='tight')
        fig.savefig(ROOT/'Figure1_Oracle_risk_coverage_P2.png',dpi=230,bbox_inches='tight')
        plt.close(fig)
        # Zoom for supplementary check - omit low-H threshold marker outside range.
        fig,ax=plt.subplots(figsize=(9.0,4.9),dpi=140)
        fig.subplots_adjust(left=.10,right=.97,top=.83,bottom=.20)
        ax.fill_between(X,0,RC,color='#e1e4e8')
        ax.fill_between(X,RC,RS,facecolor='#f7ecec',edgecolor='#b75e63',hatch='///',lw=0)
        ax.plot(X,RS,color='#9f1722',lw=2.2,label='Full-sequence oracle')
        ax.plot(X,RC,color='#215184',lw=2.0,ls=(0,(5,2)),label='Case-wise oracle')
        ax.axhline(5,color='#666666',lw=.85,ls=':')
        ax.scatter([float(seq['fraction'])*100,float(case['fraction'])*100],[5,5],s=38,c=['#9f1722','#215184'])
        ax.scatter([50],[5],marker='X',s=105,color='#202124',zorder=10)
        ax.scatter([50],[float(r_seq)*100],marker='s',s=34,color='#9f1722',zorder=10)
        ax.annotate('Full-sequence boundary: 6.54% risk\nat 50% coverage',xy=(50,float(r_seq)*100),xytext=(53,9.0),fontsize=9,color='#9f1722',arrowprops={'arrowstyle':'-','lw':.9,'color':'#9f1722'})
        ax.text(13,10.4,'Best low-H threshold: 20.59% risk (above this zoom)',fontsize=9,color='#2b7651')
        ax.set(xlim=(0,100),ylim=(0,12),xlabel='Expected coverage (%)',ylabel='Expected selective risk (%)')
        ax.set_xticks(np.arange(0,101,10));ax.set_yticks(np.arange(0,13,2))
        ax.grid(color='#cdd0d7',ls=':',lw=.65);ax.set_axisbelow(True)
        fig.suptitle('Close-up: 5% risk / 50% coverage requirement',fontsize=14,weight='bold',y=.97)
        fig.legend(loc='upper center',bbox_to_anchor=(.5,.925),ncol=2,frameon=False,fontsize=9.2)
        fig.text(.5,.048,'Fixed-cohort oracle; ground-truth-informed, randomized acceptance. No prospective risk guarantee.',ha='center',fontsize=8.6)
        fig.savefig(ROOT/'Figure1_Oracle_risk_coverage_zoom_P2.pdf',bbox_inches='tight')
        fig.savefig(ROOT/'Figure1_Oracle_risk_coverage_zoom_P2.png',dpi=225,bbox_inches='tight')
        plt.close(fig)
    print(json.dumps(summary,ensure_ascii=False,indent=2))
if __name__=='__main__': main()

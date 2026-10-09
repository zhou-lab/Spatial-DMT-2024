#!/usr/bin/env python3
"""Shared-CpG/depth matching and local beta-binomial scar sensitivity analysis."""
from pathlib import Path
import subprocess
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.special import betaln, digamma, expit, logit
from scipy.stats import t
from scipy.stats import false_discovery_control
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
RUN = Path('/scr1/users/suc3/20260528_runNewSpatialPipe')
OUT = ROOT / 'results/figure2/scar_matched_coverage'


def block_ci(values, positions, seed=42):
    # Resample contiguous physical blocks, retaining gaps in the coordinate axis.
    groups = (np.asarray(positions)-1)//5
    blocks = [np.asarray(values)[groups == g] for g in np.unique(groups)]
    rng = np.random.default_rng(seed)
    estimates = [np.median(np.concatenate([blocks[i] for i in
                 rng.integers(len(blocks), size=len(blocks))])) for _ in range(2000)]
    return np.quantile(estimates, [.025, .975])


def bb_model(m, n, positions):
    """Triplet fixed intercept, common signed flank slope and scar effect.

    beta-binomial overdispersion; uncertainty uses contiguous five-position
    spatial blocks with a small-cluster t reference. Exploratory within-section
    inference, not independent biological replication.
    """
    k = len(m)
    X = np.zeros((3*k, k+2))
    X[np.arange(3*k), np.repeat(np.arange(k), 3)] = 1
    X[:, k] = np.tile([0, 1, 0], k)
    X[:, k+1] = np.tile([-1, 0, 1], k)
    y, total = m.ravel(), n.ravel()

    def calc(theta, per_obs=False):
        p = expit(X @ theta[:-1])
        concentration = np.exp(theta[-1])
        a, b = p*concentration, (1-p)*concentration
        ll = betaln(y+a, total-y+b)-betaln(a,b)
        da = digamma(y+a)-digamma(total+concentration)-digamma(a)+digamma(concentration)
        db = digamma(total-y+b)-digamma(total+concentration)-digamma(b)+digamma(concentration)
        score = np.column_stack((X*(concentration*p*(1-p)*(da-db))[:,None],
                                 concentration*(p*da+(1-p)*db)))
        return score if per_obs else (-ll.mean(), -score.mean(axis=0))

    start = np.r_[logit((m.sum(axis=1)+.5)/(n.sum(axis=1)+1)), 0., 0., np.log(100)]
    fit = minimize(calc, start, jac=True, method='L-BFGS-B',
                   bounds=[(-12,12)]*(k+2)+[(-4,16)],
                   options={'maxiter':3000, 'ftol':1e-12, 'gtol':1e-5})
    if not fit.success:
        raise RuntimeError(fit.message)
    theta=fit.x
    step=1e-4
    H=np.column_stack([(calc(theta+np.eye(len(theta))[j]*step)[1]-
                        calc(theta-np.eye(len(theta))[j]*step)[1])/(2*step)
                       for j in range(len(theta))])
    H=(H+H.T)/2 * len(y)
    inv=np.linalg.pinv(H)
    scores=calc(theta, True)
    groups=np.repeat((np.asarray(positions)-1)//5,3)
    unique=np.unique(groups)
    S=np.array([scores[groups==g].sum(axis=0) for g in unique])
    cov=inv @ (S.T@S) @ inv * len(unique)/(len(unique)-1)
    se=np.sqrt(max(cov[k,k],0))
    crit=t.ppf(.975,len(unique)-1)
    low,high=theta[k]+np.array([-1,1])*crit*se
    base=expit(theta[:k])
    effect=lambda g: float(np.mean(expit(theta[:k]+g)-base))
    return dict(model_log_odds=theta[k],model_log_odds_ci_low=low,
                model_log_odds_ci_high=high,model_beta_delta=effect(theta[k]),
                model_beta_ci_low=effect(low),model_beta_ci_high=effect(high),
                model_p=2*t.sf(abs(theta[k]/se),len(unique)-1),
                model_rho=1/(1+np.exp(theta[-1])),spatial_blocks=len(unique),
                optimizer_success=fit.success,model_max_gradient=float(abs(calc(theta)[1]).max()))


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    binary=OUT/'extract_scar_sites'
    subprocess.run(['cc','-O3',str(ROOT/'tools/figure2/extract_scar_sites.c'),'-lz','-o',str(binary)],check=True)
    previous=pd.read_csv(ROOT/'results/figure2/scar_methylation_orthogonal/orthogonal_scar_methylation_metrics.tsv',sep='\t')
    results=[]; pixels=[]
    for row in previous.itertuples():
        sample=row.sample
        print('Processing',sample,flush=True)
        cg=RUN/'final_output'/sample/f'{sample}.cg'
        if sample.startswith('SpMET'):
            idx=pd.read_csv(str(cg)+'.idx',sep='\t',header=None,dtype=str)
            assert idx[1].is_unique and pd.to_numeric(idx[1]).is_monotonic_increasing
            frame=pd.read_csv(ROOT/'results/figure2/scar_methylation_orthogonal/cache'/f'{sample}_stream.tsv',sep='\t')
            assert len(idx)==len(frame)
            frame['x']=idx[0].str[:2].astype(int);frame['y']=idx[0].str[2:].astype(int)
        else:
            frame=pd.read_csv(ROOT/'results/figure2/smcseq_paper_comparison'/f'{sample}_matched.tsv',sep='\t')
            frame=frame.loc[frame.retained].copy()
            frame['x']=frame.paper_B;frame['y']=frame.paper_A
        axis='x' if row.scar_axis=='vertical' else 'y'
        other='y' if axis=='x' else 'x'
        line=row.scar_index
        lookup={(int(r[axis]),int(r[other])):int(r['record']) for _,r in frame.iterrows()}
        trips=[]
        for pos in sorted(frame[other].unique()):
            keys=[(line+d,int(pos)) for d in [-1,0,1]]
            if all(key in lookup for key in keys):trips.append((int(pos),[lookup[key] for key in keys]))
        selected=sorted({rec for _,trip in trips for rec in trip})
        ids=OUT/f'{sample}_records.txt';ids.write_text('\n'.join(map(str,selected))+'\n')
        cache=OUT/f'{sample}_sites.bin'
        if not cache.exists():
            temp=cache.with_suffix('.tmp')
            with temp.open('wb') as f:subprocess.run([str(binary),str(cg),str(ids)],stdout=f,check=True)
            temp.rename(cache)
        data={}
        with cache.open('rb') as f:
            while True:
                header=np.fromfile(f,dtype='<u8',count=2)
                if len(header)==0:break
                rec,size=map(int,header)
                data[rec]=np.fromfile(f,dtype='<u8',count=3*size).reshape(-1,3).astype(np.int64)
        # Validate record identities and M/U decoding against independently cached summaries.
        reference=frame.set_index('record')
        for rec,a in data.items():
            assert len(a)==int(reference.loc[rec,'covered_cpgs'])
            assert a[:,1].sum()==int(reference.loc[rec,'methylated_calls'])
            assert a[:,1:].sum()==int(reference.loc[rec,'total_calls'])
            assert np.all(np.diff(a[:,0])>0)
        rng=np.random.default_rng(20260918)
        ms=[];ns=[];positions=[];shared_delta=[];equal_delta=[];mc=[];equal_ms=[];equal_ns=[]
        for pos,trip in trips:
            arrays=[data[r] for r in trip]
            sites=np.intersect1d(np.intersect1d(arrays[0][:,0],arrays[1][:,0]),arrays[2][:,0])
            if len(sites)<100:continue
            matched=[a[np.searchsorted(a[:,0],sites)] for a in arrays]
            M=np.array([a[:,1] for a in matched]);N=np.array([a[:,1]+a[:,2] for a in matched])
            minimum=N.min(axis=0)
            ratios=M.sum(axis=1)/N.sum(axis=1)
            delta=ratios[1]-(ratios[0]+ratios[2])/2
            # Expected hypergeometric downsample gives exact equal-depth estimate.
            equal=(M/N*minimum).sum(axis=1)/minimum.sum()
            ed=equal[1]-(equal[0]+equal[2])/2
            draws=[]
            for draw in range(100):
                dm=np.array([rng.hypergeometric(M[j],N[j]-M[j],minimum).sum() for j in range(3)])
                if draw==0:
                    equal_ms.append(dm)
                    equal_ns.append(np.repeat(minimum.sum(),3))
                draws.append((dm[1]-(dm[0]+dm[2])/2)/minimum.sum())
            ms.append(M.sum(axis=1));ns.append(N.sum(axis=1));positions.append(pos)
            shared_delta.append(delta);equal_delta.append(ed);mc.append(draws)
            pixels.append(dict(sample=sample,position=pos,shared_cpgs=len(sites),
                matched_calls_per_pixel=int(minimum.sum()),shared_delta=delta,equal_depth_delta=ed,
                downsample_mean=np.mean(draws),downsample_sd=np.std(draws),
                before_m=int(M[0].sum()),scar_m=int(M[1].sum()),after_m=int(M[2].sum()),
                before_n=int(N[0].sum()),scar_n=int(N[1].sum()),after_n=int(N[2].sum()),
                equal_before_m=int(equal_ms[-1][0]),equal_scar_m=int(equal_ms[-1][1]),
                equal_after_m=int(equal_ms[-1][2])))
        ci=block_ci(equal_delta,positions)
        result=dict(sample=sample,scar_axis=row.scar_axis,scar_index=line,triplets=len(positions),
            original_delta=row.median_paired_beta_delta,shared_cpg_delta=np.median(shared_delta),
            equal_depth_delta=np.median(equal_delta),equal_depth_ci_low=ci[0],equal_depth_ci_high=ci[1],
            downsample_median_mean=np.median(np.array(mc),axis=0).mean(),
            **bb_model(np.array(ms),np.array(ns),positions),
            **{'depth_matched_'+k:v for k,v in bb_model(np.array(equal_ms),np.array(equal_ns),positions).items()})
        results.append(result)
        pd.DataFrame(results).to_csv(OUT/'matched_coverage_results.tsv',sep='\t',index=False)
        pd.DataFrame(pixels).to_csv(OUT/'matched_triplets.tsv',sep='\t',index=False)
        print(result,flush=True)
    qvalues=false_discovery_control([r['depth_matched_model_p'] for r in results])
    for r,q in zip(results,qvalues):r['depth_matched_model_q']=q
    pd.DataFrame(results).to_csv(OUT/'matched_coverage_results.tsv',sep='\t',index=False)
    fig,ax=plt.subplots(figsize=(9,5),layout='constrained')
    for j,r in enumerate(results):
        ax.scatter(r['original_delta']*100,j+.18,c='#888888',s=25,label='Original' if j==0 else None)
        for offset,prefix,color,label in [(0,'equal_depth','#2678A5','Shared CpGs + equal depth'),
                                         (-.18,'depth_matched_model_beta','#B64645','Depth-matched beta-binomial')]:
            delta=r[prefix+'_delta'] if prefix=='equal_depth' else r[prefix+'_delta']
            lo,hi=r[prefix+'_ci_low'],r[prefix+'_ci_high']
            ax.errorbar(delta*100,j+offset,xerr=[[max(0,delta-lo)*100],[max(0,hi-delta)*100]],
                        fmt='o',color=color,ms=4,label=label if j==0 else None)
    ax.axvline(0,color='black',lw=.8);ax.axvspan(-2,2,color='#dddddd',alpha=.35)
    ax.set_yticks(range(len(results)),[r['sample'] for r in results]);ax.invert_yaxis()
    ax.set_xlabel('Scar − adjacent methylation (percentage points; 95% interval)')
    ax.set_title('Scar methylation after matching loci and coverage')
    ax.legend(fontsize=8);fig.savefig(OUT/'Scar_matched_coverage.png',dpi=300)

if __name__=='__main__':main()

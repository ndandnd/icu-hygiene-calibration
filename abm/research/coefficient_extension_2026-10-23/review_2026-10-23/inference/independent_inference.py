"""Independent read-only numerical audit; imports no production analysis modules.

Reuses only immutable record/seed/assay files. Outputs stay inside this review.
"""
from pathlib import Path
import hashlib
import json
import math
import numpy as np

HERE = Path(__file__).resolve().parent
STUDY = HERE.parents[1]
RESEARCH = STUDY.parent
SEEDBOOK = RESEARCH / 'full_icu_confirmation_2026-10-08/bridge/seedbook.json'
NBOOT = 20000

def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def percentile(x):
    finite = np.asarray(x)[np.isfinite(x)]
    return {'interval': np.quantile(finite, [.025, .975]).tolist() if len(finite) else None,
            'defined': len(finite), 'undefined': len(x)-len(finite)}

def renewal(h, rho, e):
    # Direct infinite-series recurrence, stopping only after remaining terms
    # are below 1e-15 relative to the accumulated sum.
    P = np.array([[1-(1-rho)*h, (1-rho)*h],
                  [(1-rho)*(1-h), 1-(1-rho)*(1-h)]])
    D = np.diag([1, 1-e])
    v = np.array([1-h,h]) @ D @ P @ D
    step = P @ D @ P @ D
    ans = 0.
    for k in range(100000):
        add = v.sum()
        ans += add
        if add <= 1e-15 * ans:
            return ans, k+1
        v = v @ step
    raise ArithmeticError('renewal recurrence failed to converge')

def exact_feasible_domain(lim):
    """Analytic domain check for the continuous q1,q2 confidence set."""
    l1,u1,l2,u2=lim
    h=.9
    a=max((1-u1)/h, 0., (1-math.sqrt(u2))/h)
    b=min((1-l1)/h, 1.)
    if l2>1-h:
        b=min(b, 1-math.sqrt(max(0., (l2-(1-h))/h)))
    if a>b:
        return {'empty':True, 'out':False}
    ec=(1-u2)/h
    vals=[a,b]+([ec] if a<=ec<=b else [])
    rmax=min(1., max((u2-(1-h*e)**2)/(e*e*h*(1-h)) for e in vals))
    return {'empty':False, 'out':a<.5 or rmax>.95, 'e_low':a,'e_high':b,'rho_max':rmax}

def main():
    seeds=json.loads(SEEDBOOK.read_text())['all_seeds']
    manifest=json.loads((STUDY/'records/manifest.json').read_text())
    params=manifest['science']['pairs']
    ids=[p['id'] for p in params]
    loc={pid:i for i,pid in enumerate(ids)}
    K=np.empty((len(seeds),len(ids),2))
    hashes={}
    for k,seed in enumerate(seeds):
        path=STUDY/f'records/{seed}.json'
        row=json.loads(path.read_text())
        assert row['seed']==seed and row['surface']=='all' and row['gate']['passed']
        assert set(row['pairs'])==set(ids)
        for j,p in enumerate(params):
            arm=row['pairs'][p['id']]
            assert all(arm['parameters'][x]==p[x] for x in ('h0','h1','rho','e'))
            K[k,j]=arm['baseline']['K'],arm['policy']['K']
        hashes[str(seed)]=digest(path)
    means=np.array([[math.fsum(K[:,j,a])/len(seeds) for a in (0,1)] for j in range(len(ids))])
    F=1-means[:,1]/means[:,0]
    index=np.random.default_rng(20261023).integers(len(seeds),size=(NBOOT,len(seeds)))
    # Independently use direct indexed sample sums, not the production weighted
    # matrix multiplication. 182 small iterations avoid a multi-GB temporary.
    BF=np.empty((NBOOT,len(ids)))
    for j in range(len(ids)):
        sums=K[:,j,:][index].sum(axis=1)
        BF[:,j]=1-sums[:,1]/sums[:,0]
    def lookup(h0,h1,r,e):
        for p in params:
            if (p['h0'],p['h1'],p['rho'],p['e'])==(h0,h1,r,e):
                return loc[p['id']]
        raise KeyError((h0,h1,r,e))
    def diff(a,b):
        xa,ya=K[:,a,0],K[:,a,1]
        xb,yb=K[:,b,0],K[:,b,1]
        psi=-(ya-(1-F[a])*xa)/means[a,0]+(yb-(1-F[b])*xb)/means[b,0]
        se=np.std(psi,ddof=1)/math.sqrt(len(seeds))
        return {'D':float(F[a]-F[b]),'F_a':float(F[a]),'F_b':float(F[b]),
                'normal_95':[(F[a]-F[b])-1.96*se,(F[a]-F[b])+1.96*se],
                'bootstrap_95':percentile(BF[:,a]-BF[:,b]),'se':float(se)}
    eg=sorted({p['e'] for p in params if p['h0']==.9 and p['h1']==.95 and p['rho'] in (0,.5)})
    eg=[e for e in eg if all(any((p['h0'],p['h1'],p['rho'],p['e'])==(.9,.95,r,e) for p in params) for r in (0,.5))]
    ix0=[lookup(.9,.95,0,e) for e in eg]
    ix5=[lookup(.9,.95,.5,e) for e in eg]
    D=F[ix5]-F[ix0]
    def roots(v):
        out=[]
        for i in range(len(eg)-1):
            if (v[i]>0)!=(v[i+1]>0):
                out.append(eg[i]+(eg[i+1]-eg[i])*v[i]/(v[i]-v[i+1]))
        return out
    bootroot=[]; count={}
    for d in BF[:,ix5]-BF[:,ix0]:
        rr=roots(d)
        count[len(rr)]=count.get(len(rr),0)+1
        bootroot.append(rr[0] if rr else np.nan)
    pred1={'e_grid':eg,'D':D.tolist(),'roots':roots(D),'root_bootstrap_95':percentile(bootroot),
           'bootstrap_number_of_roots_histogram':count,'holds':.91<=roots(D)[0]<=.935}
    pred2={}
    for h0,h1 in ((.6,.7),(.7,.8),(.8,.9),(.9,.95),(.95,.99)):
        for e in (.8,1.):
            pred2[f'{h0:g}->{h1:g};e={e:g}']=diff(lookup(h0,h1,.5,e),lookup(h0,h1,0,e))
    pred3={str(e):diff(loc[f'x_trans_c0.9_to_0.95_rho0.5_e{e:g}'],loc[f'x_trans_c0.9_to_0.95_rho0_e{e:g}']) for e in (1.,.8)}
    renewal_rows={}
    for p in params:
        if p['group']!='fine_grid':continue
        s0,n0=renewal(p['h0'],p['rho'],p['e']);s1,n1=renewal(p['h1'],p['rho'],p['e'])
        j=loc[p['id']]
        renewal_rows[p['id']]={'F_geo':1-s1/s0,'F_star':F[j], 'error':1-s1/s0-F[j], 'terms':max(n0,n1)}
    worst=max(renewal_rows,key=lambda p:abs(renewal_rows[p]['error']))
    pred4={'worst_row':worst,'max_abs_error_pp':100*abs(renewal_rows[worst]['error']),
           'holds':abs(renewal_rows[worst]['error'])<=.035,'rows':renewal_rows}
    R=sorted({p['rho'] for p in params if p['group']=='fine_grid'})
    E=sorted({p['e'] for p in params if p['group']=='fine_grid'})
    Z=np.array([[F[lookup(.9,.95,r,e)] for e in E] for r in R])
    def interp(pts):
        # Sequential 1-D interpolation is independently equivalent to bilinear
        # interpolation and uses no production surface implementation.
        out=np.empty(len(pts))
        for e in np.unique(pts[:,1]):
            chosen=pts[:,1]==e
            vertical=np.array([np.interp(e,E,row) for row in Z])
            out[chosen]=np.interp(pts[chosen,0],R,vertical)
        return out
    holdout=[]
    for p in params:
        if p['group']=='crossing' and p['e'] not in E:
            value=float(interp(np.array([[p['rho'],p['e']]]))[0])
            holdout.append({'id':p['id'],'interpolated_F':value,'direct_F':F[loc[p['id']]],
                            'interpolation_error_pp':100*(value-F[loc[p['id']]])})
    uq=RESEARCH/'robustness_uq_2026-10-17/uq/run_001/raw_benchmark.npz'
    assay=np.load(uq)
    projected={}
    for truth in ('independent','clustered','imperfect','combined'):
        lims=assay[truth+'_n8000_confidence_limits']
        widths=[];status=[];analytic=[]
        for lim in lims:
            l1,u1,l2,u2=lim
            a=max((1-u1)/.9,1e-9); b=min((1-l1)/.9,1.)
            rows=[]
            if a<=b:
                for e in np.linspace(a,b,81):
                    base=(1-.9*e)**2; scale=e**2*.9*.1
                    lo=max((l2-base)/scale,0.);hi=min((u2-base)/scale,1.)
                    if lo<=hi:
                        rows.append(np.column_stack([np.linspace(lo,hi,61),np.full(61,e)]))
            an=exact_feasible_domain(lim)
            analytic.append('empty' if an['empty'] else ('out' if an['out'] else 'defined'))
            if not rows:
                status.append('empty');continue
            pts=np.concatenate(rows)
            if pts[:,0].max()>.95 or pts[:,1].min()<.5:
                status.append('out');continue
            z=interp(pts)
            widths.append(float(z.max()-z.min()));status.append('defined')
        projected[truth]={'median_width_pp':100*float(np.median(widths)),
            'counts':{s:status.count(s) for s in ('defined','empty','out')},
            'analytic_domain_counts':{s:analytic.count(s) for s in ('defined','empty','out')},
            'domain_classification_disagreements':sum(x!=y for x,y in zip(status,analytic))}
    checks={}
    summary_path=STUDY/'summary.json'
    if summary_path.exists():
        orig=json.loads(summary_path.read_text())
        checks['all_182_pooled_response_max_abs_error']=max(abs(F[loc[p]]-v['limiting_response']) for p,v in orig['pairs'].items())
        checks['all_182_bootstrap_endpoint_max_abs_error']=max(abs(x-y) for p,v in orig['pairs'].items()
            for x,y in zip(percentile(BF[:,loc[p]])['interval'],v['F_bootstrap_95']['interval']))
        checks['crossing_point_abs_error']=abs(pred1['roots'][0]-orig['predictions']['prediction_1_crossing']['e_star'])
        checks['crossing_bootstrap_max_endpoint_error']=max(abs(x-y) for x,y in zip(pred1['root_bootstrap_95']['interval'],orig['predictions']['prediction_1_crossing']['e_star_bootstrap_95']['interval']))
        checks['renewal_formula_max_abs_error']=max(abs(v['F_geo']-orig['predictions']['prediction_4_renewal_comparator']['rows'][p]['F_geo']) for p,v in renewal_rows.items())
        checks['projection_median_max_abs_pp_error']=max(abs(v['median_width_pp']-orig['secondary_exact_icu_projection']['by_truth'][p]['median_width_pp']) for p,v in projected.items())
    result={'n_seeds':len(seeds),'n_pairs':len(ids),'bootstrap_index_sha256':hashlib.sha256(index.astype('<i8').tobytes()).hexdigest(),
            'prediction_1':pred1,'prediction_2':pred2,'prediction_3':pred3,'prediction_4':pred4,
            'projection':projected,'interpolation_holdout_checks':holdout,
            'comparison_with_production_summary':checks,
            'inputs':{'records':hashes,'seedbook_sha256':digest(SEEDBOOK),'assay_sha256':digest(uq),'manifest_sha256':digest(STUDY/'records/manifest.json')},
            'runtime':{'numpy':np.__version__}}
    np.savez_compressed(HERE/'independent_arrays.npz',K=K,F=F,bootstrap_F=BF,ids=ids,seeds=seeds)
    default=lambda obj:obj.item() if isinstance(obj,np.generic) else str(obj)
    (HERE/'independent_results.json').write_text(json.dumps(result,indent=2,allow_nan=False,default=default)+'\n')
    print(json.dumps({k:result[k] for k in ('prediction_1','prediction_3','projection','interpolation_holdout_checks','comparison_with_production_summary')},indent=2,default=default))

if __name__=='__main__':main()

"""Saved coefficient algebra only. No ABM imports, simulation, or resampling."""
from pathlib import Path
from collections import defaultdict, Counter
from fractions import Fraction as F
from itertools import combinations
import hashlib, json, math

ROOT=Path(__file__).resolve().parents[3]
OUT=Path(__file__).resolve().parent
HASHES={}
def sha(p):
    p=Path(p); value=hashlib.sha256(p.read_bytes()).hexdigest();HASHES[str(p)]=value;return value
def canonical(x): return json.dumps(x,sort_keys=True,separators=(',',':'),allow_nan=False)
def digest(x): return hashlib.sha256(canonical(x).encode()).hexdigest()
def read(p): sha(p);return json.loads(Path(p).read_text())
def sealed(r): assert digest({k:v for k,v in r.items() if k!='payload_sha256'})==r['payload_sha256']
def frac(x): return F(str(x))
SEEDBOOK=ROOT/'research/full_icu_confirmation_2026-10-08/bridge/seedbook.json'
book=read(SEEDBOOK); seeds=book['all_seeds'];assert len(seeds)==len(set(seeds))==256
banks={s:b for b,ss in book['seed_blocks'].items() for s in ss}
EXT=ROOT.parent/'research/author_notes_2026-10-10/collections/bridge_20261011_061144_utc/results/bridge'
EXTM=EXT/'part_00/manifest.json'; extended_manifest=read(EXTM)
assert sha(EXTM)=='e50eaefb3f8aac31600abcd4fda1cf3b5d357f6e79e701bdcf59560d54bc0ed9'
PROFILES={'independent_perfect':(0,1),'clustered_perfect':(.5,1),'independent_imperfect':(0,.8),'combined':(.5,.8)}
PAIRS={name:f'h0.9_to_0.95_rho{rho:g}_e{e:g}' for name,(rho,e) in PROFILES.items()}
validation={}; results={}; worker_data={}
for workflow,base in [('original','coefficient_extension_2026-10-23'),('alternative','workflow_coefficient_2026-10-21')]:
    directory=ROOT/'research'/base/'records'; manifest=read(directory/'manifest.json')
    assert digest({k:v for k,v in manifest.items() if k!='identity_sha256'})==manifest['identity_sha256']
    assert digest(manifest['science'])==manifest['science_sha256']
    config=manifest['science']['resolved_config'];assert set(config['seeds'])==set(seeds)
    assert config['warmup_ticks']==90*2592 and config['observation_ticks']==56*2592
    assert config['h0']==.9 and config['h1']==.95 and config['opportunities']=='zone_separate'
    assert config['model_overrides']['prob_patient_sick']==.03
    source_map=manifest['source_sha256']
    for p,h in source_map.items():
        path=EXTM if p=='extended_manifest(part_00/manifest.json)' else ROOT/p
        assert sha(path)==h,(path,'source drift')
    files=list(directory.glob('*.json')); records=[p for p in files if p.name!='manifest.json']
    assert {int(p.stem) for p in records}==set(seeds) and len(records)==256
    pool={name:{arm:defaultdict(F) for arm in ('baseline','policy')} for name in PAIRS}
    sums=defaultdict(float);max_rel=0.;counts=Counter();refs_verified=0
    missing=[]
    for p in sorted(records):
        r=read(p);sealed(r);seed=int(p.stem)
        assert r['seed']==seed and r['status']=='complete' and r['bank']==banks[seed]
        assert r['identity_sha256']==manifest['identity_sha256'] and r['science_sha256']==manifest['science_sha256']
        assert r['workflow']['alternative']==(workflow=='alternative')
        assert r['executed_ticks']==146*2592 and r['model_kwargs']['seed']==seed
        assert r['model_kwargs']['initial_environment_contaminated'] is False
        expected_sources=dict(source_map)
        if workflow=='alternative':
            extra='research/workflow_coefficient_2026-10-21/code/run_batch.py';expected_sources[extra]=sha(ROOT/extra)
        assert r['source_sha256']==expected_sources
        if workflow=='original':
            assert r['gate']['passed'] and r['gate']['event_hash']['equal']
            ref=r['gate']['archived_extended_record'];rp=Path(ref['path']); old=read(rp);sealed(old)
            assert sha(rp)==ref['file_sha256'] and old['seed']==seed and old['status']=='complete'
            assert old['ordered_bridge_events_sha256']==r['ordered_bridge_events_sha256']
            for key,h in r['code_sha256'].items(): assert sha(ROOT/key)==h
            refs_verified+=1
        for name,pair in PAIRS.items():
            assert pair in r['pairs'],(workflow,pair)
            par=r['pairs'][pair]['parameters'];rho,e=PROFILES[name]
            assert par['h0']==.9 and par['h1']==.95 and par['rho']==rho and par['e']==e
            for arm in ('baseline','policy'):
                row=r['pairs'][pair][arm]; workers=row.get('K_by_worker')
                if not isinstance(workers,dict):missing.append([seed,name,arm]);continue
                assert workers and all(math.isfinite(v) and v>=0 for v in workers.values())
                total=math.fsum(workers.values());K=row['K'];rel=abs(total-K)/max(1.,abs(K));max_rel=max(max_rel,rel)
                assert math.isclose(total,K,rel_tol=2e-10,abs_tol=1e-8),(p,name,arm,total,K)
                if workflow=='original': assert old['pairs'][pair][arm]['K_by_worker']==workers
                for worker,value in workers.items():pool[name][arm][worker]+=frac(value)
                sums[(name,arm)]+=K; counts[(name,arm)]+=1
    if missing:
        validation[workflow]={'status':'STOP_MISSING_DICTIONARIES','missing':missing};continue
    assert all(v==256 for v in counts.values())
    validation[workflow]={'status':'PASS','records':256,'banks':dict(Counter(banks.values())),
      'manifest_identity':manifest['identity_sha256'],'sources_verified':len(source_map),
      'sealed_payloads_verified':256,'original_reference_records_verified':refs_verified,
      'max_worker_sum_relative_error':max_rel,'same_seed_set_as_seedbook':True,
      'reference_binding': 'Original: per-record archived SHA/payload/event hash and selected worker dictionaries match; alternative: source/config/seed identity and archived extended-manifest source binding, intentionally different event histories.'}
    results[workflow]={};worker_data[workflow]={}
    for name,arms in pool.items():
        workers=sorted(set(arms['baseline'])|set(arms['policy']))
        a=[arms['baseline'][u] for u in workers];b=[arms['policy'][u] for u in workers]
        assert all(x>0 or y==0 for x,y in zip(a,b))
        data=[(u,x,y) for u,x,y in zip(workers,a,b) if x>0]; data.sort(key=lambda z:z[2]/z[1])
        k0=sum(a);k1=sum(b);nominal=1-k1/k0;rmin=data[0][2]/data[0][1];rmax=data[-1][2]/data[-1][1]
        worker_data[workflow][name]={u:{'baseline':str(x),'policy':str(y)} for u,x,y in data}
        entry={'pair':PAIRS[name],'workers':len(data),'nominal_F':float(nominal),'nominal_F_exact':str(nominal),
               'K0_mean_from_workers':float(k0/256),'K1_mean_from_workers':float(k1/256),
               'worker_ratio_min':float(rmin),'worker_ratio_max':float(rmax),'boxes':[]}
        for m,M in [(F(1,2),F(2)),(F(1,4),F(4))]:
            candidates=[]
            # Any extremum uses endpoint factors sorted by ratio around the optimum.
            for cut in range(len(data)+1):
                for low_heavy in (True,False):
                    factors=[(M if ((i<cut)==low_heavy) else m) for i in range(len(data))]
                    den=sum(z[1]*t for z,t in zip(data,factors));num=sum(z[2]*t for z,t in zip(data,factors))
                    candidates.append((num/den, factors))
            mn,tmin=min(candidates,key=lambda z:z[0]);mx,tmax=max(candidates,key=lambda z:z[0])
            # Verify exact fractional first-order certificate min_t Σt(b−x a)=0.
            assert sum((m if y-mn*x>=0 else M)*(y-mn*x) for _,x,y in data)==0
            assert sum((M if y-mx*x>=0 else m)*(y-mx*x) for _,x,y in data)==0
            lower,upper=1-mx,1-mn
            # These chosen symmetric boxes have sqrt(M/m)=M, so TV radius=(M−1)/(M+1).
            radius=(rmax-rmin)*(M-1)/(M+1)
            assert lower<=nominal<=upper and max(nominal-lower,upper-nominal)<=radius
            entry['boxes'].append({'m':float(m),'M':float(M),'F_min':float(lower),'F_max':float(upper),
              'F_min_exact':str(lower),'F_max_exact':str(upper),'radius':float(radius),
              'F_min_worker_factors':dict(zip([z[0] for z in data],map(str,tmax))),
              'F_max_worker_factors':dict(zip([z[0] for z in data],map(str,tmin)))})
        results[workflow][name]=entry
    contrasts=[]
    for p,q in combinations(PAIRS,2):
        ep,eq=results[workflow][p],results[workflow][q]
        for bp,bq in zip(ep['boxes'],eq['boxes']):
            lo=bp['F_min']-bq['F_max'];hi=bp['F_max']-bq['F_min'];gap=ep['nominal_F']-eq['nominal_F']
            contrasts.append({'P':p,'Q':q,'m':bp['m'],'M':bp['M'],'nominal_P_minus_Q':gap,
               'conservative_shared_t_contrast_min':lo,'conservative_shared_t_contrast_max':hi,
               'sign':'positive' if lo>0 else ('negative' if hi<0 else 'unresolved'),
               'radius_sign_certificate': abs(gap)>bp['radius']+bq['radius']})
    results[workflow]['contrasts']=contrasts
output={'status':'PASS' if all(v['status']=='PASS' for v in validation.values()) else 'PARTIAL',
        'scope':'Post hoc conditional finite-horizon leading-coefficient worker reweighting only; no Monte Carlo uncertainty, new backgrounds, resampling, or clinical validation.',
        'validation':validation,'results':results}
(OUT/'results.json').write_text(json.dumps(output,indent=2)+'\n')
(OUT/'pooled_worker_coefficients.json').write_text(json.dumps(worker_data,indent=2)+'\n')
sha(OUT/'PLAN.md');sha(Path(__file__));
(OUT/'input_hashes.json').write_text(json.dumps(HASHES,sort_keys=True,indent=2)+'\n')
lines=['# Saved-data worker heterogeneity robustness','',output['scope'],'',
       'Exact fractional optimization treats the saved decimal worker coefficients as exact inputs. Numerical precision of the original coefficient engine is not bounded by this calculation. All seed worker sums agree with recorded K to the tolerance reported in results.json.','',
       '| Workflow | Profile | Nominal F (%) | t in [.5,2] F (%) | t in [.25,4] F (%) |',
       '|---|---|---:|---:|---:|']
for workflow,entries in results.items():
 for name in PAIRS:
  if name not in entries:continue
  e=entries[name];bs=e['boxes'];lines.append(f"| {workflow} | {name} | {100*e['nominal_F']:.4f} | [{100*bs[0]['F_min']:.4f}, {100*bs[0]['F_max']:.4f}] | [{100*bs[1]['F_min']:.4f}, {100*bs[1]['F_max']:.4f}] |")
lines+=['','## Conservative shared-multiplier contrasts','','P−Q denotes prevention difference in percentage points. Each mechanism is optimized separately, so the interval contains every shared-t contrast but need not be attainable as a joint interval. A zero-containing interval is unresolved, not evidence of reversal.','', '| Workflow | P | Q | Bounds | Contrast interval (pp) | Guaranteed sign |','|---|---|---|---|---:|---|']
for workflow,entries in results.items():
 for c in entries.get('contrasts',[]):lines.append(f"| {workflow} | {c['P']} | {c['Q']} | [{c['m']}, {c['M']}] | [{100*c['conservative_shared_t_contrast_min']:.4f}, {100*c['conservative_shared_t_contrast_max']:.4f}] | {c['sign']} |")
lines+=['','## Mathematics and limits','','For worker u, let A_u and B_u be pooled baseline and policy coefficients. F(t)=1−Σt_u B_u/Σt_u A_u. Sorting B_u/A_u and evaluating all threshold endpoint assignments yields the exact linear-fractional box extrema; the script verifies the root certificates with rational arithmetic. Multipliers apply to the same worker ID across all 256 seeds and both arms; their overall scale cancels. This is a realizable subclass of directional heterogeneity: assign a worker-specific pickup factor a_u and transfer factor b_u with a_u b_u=t_u. It does not vary a_u/b_u separately, because that directional ratio cancels at this order.','',
       'Workers with zero pooled baseline and policy coefficient contribute nothing and are omitted. A positive policy coefficient with zero baseline is rejected. The saved worker labels are used literally; no role interpretation is required. The boxes are mathematical stress ranges with no empirical probability interpretation. Their extrema do not incorporate background Monte Carlo uncertainty, parameter estimation, finite-transmission feedback, environmental routes, endpoint-specific heterogeneity, or arbitrary variation within a worker.','',
       'This is a post hoc analysis under the saved PLAN.md. It establishes conditional robustness only over the specified worker-multiplier class on these saved backgrounds. The leading-limit invariance to global direction multipliers remains an algebraic fact; these nontrivial bounds concern heterogeneous products. Positive or negative contrast certificates use conservative per-profile interval differences. No exact shared-t contrast optimization is claimed.','',
       'Input hashes include both manifests, all 512 analysis records, all 256 original archived reference records, source/config/seed files, this script and the saved plan. The original extension is linked per seed to archived source data using archived file hash, sealed payload, event hash and exact selected worker-dictionary equality. The alternative workflow has different event histories by construction; it is linked to the same seedbook and frozen scientific/source/configuration identities, not falsely asserted to have matching histories.']
(OUT/'RESULTS.md').write_text('\n'.join(lines)+'\n')
print('\n'.join(lines[:18])); print('hashed inputs',len(HASHES))

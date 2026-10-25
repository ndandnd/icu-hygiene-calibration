#!/usr/bin/env python3
"""Portable exact replay of pooled algebra; does not validate raw provenance."""
import argparse
from fractions import Fraction as F
import hashlib
from itertools import combinations
import json
from pathlib import Path


def read(path):
    raw = path.read_bytes()
    return json.loads(raw), hashlib.sha256(raw).hexdigest()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def residual_root(a, b, m, M, minimize):
    """Exact Dinkelbach updates; each step solves the box linear subproblem."""
    x = sum(b) / sum(a)
    seen = set()
    while True:
        require(x not in seen, 'Unexpected fractional iteration cycle')
        seen.add(x)
        factors = []
        for ai, bi in zip(a, b):
            positive = bi-x*ai >= 0
            factors.append((m if positive else M) if minimize else (M if positive else m))
        residual = sum(t*(bi-x*ai) for t, ai, bi in zip(factors, a, b))
        if residual == 0:
            return x, len(seen)
        x = sum(t*bi for t, bi in zip(factors, b)) / sum(t*ai for t, ai in zip(factors, a))


def close_float(actual, exact, label):
    require(abs(actual-float(exact)) <= 1e-13, label)


def replay(pooled, reported):
    profiles = ['independent_perfect', 'clustered_perfect', 'independent_imperfect', 'combined']
    require(set(pooled) == {'original', 'alternative'}, 'Unexpected workflows')
    checked = {}
    contrast_count = box_count = 0
    for workflow in ['original', 'alternative']:
        entries = pooled[workflow]
        require(set(entries) == set(profiles), 'Incomplete profile set')
        checked[workflow] = {}
        for profile in profiles:
            workers = entries[profile]
            a = [F(v['baseline']) for v in workers.values()]
            b = [F(v['policy']) for v in workers.values()]
            require(len(a) > 0 and all(x > 0 for x in a) and all(y >= 0 for y in b), 'Invalid pooled coefficients')
            row = reported['results'][workflow][profile]
            nominal = 1-sum(b)/sum(a)
            require(nominal == F(row['nominal_F_exact']), 'Nominal exact mismatch')
            close_float(row['nominal_F'], nominal, 'Nominal float mismatch')
            boxes = []
            require(len(row['boxes']) == 2, 'Expected two boxes')
            for expected, box in zip([(F(1,2), F(2)), (F(1,4), F(4))], row['boxes']):
                m, M = F(str(box['m'])), F(str(box['M']))
                require((m, M) == expected, 'Unexpected box')
                rmin, nmin = residual_root(a, b, m, M, True)
                rmax, nmax = residual_root(a, b, m, M, False)
                lower, upper = 1-rmax, 1-rmin
                require(lower == F(box['F_min_exact']) and upper == F(box['F_max_exact']), 'Exact box mismatch')
                close_float(box['F_min'], lower, 'Lower box float mismatch')
                close_float(box['F_max'], upper, 'Upper box float mismatch')
                boxes.append({'m': str(m), 'M': str(M), 'F_min_exact': str(lower), 'F_max_exact': str(upper), 'iterations': [nmin, nmax]})
                box_count += 1
            checked[workflow][profile] = {'nominal_F_exact': str(nominal), 'workers': len(a), 'boxes': boxes}
        expected_keys = {(p,q,str(m),str(M)) for p,q in combinations(profiles,2) for m,M in [(F(1,2),F(2)),(F(1,4),F(4))]}
        seen = set()
        contrasts = []
        for row in reported['results'][workflow]['contrasts']:
            p,q = row['P'],row['Q']
            m,M = F(str(row['m'])), F(str(row['M']))
            key=(p,q,str(m),str(M))
            require(key in expected_keys and key not in seen, 'Unexpected/duplicate contrast')
            seen.add(key)
            index = 0 if m == F(1,2) else 1
            pb,qb = checked[workflow][p]['boxes'][index], checked[workflow][q]['boxes'][index]
            lower=F(pb['F_min_exact'])-F(qb['F_max_exact'])
            upper=F(pb['F_max_exact'])-F(qb['F_min_exact'])
            close_float(row['conservative_shared_t_contrast_min'],lower,'Contrast lower mismatch')
            close_float(row['conservative_shared_t_contrast_max'],upper,'Contrast upper mismatch')
            sign='positive' if lower>0 else 'negative' if upper<0 else 'unresolved'
            require(row['sign']==sign,'Contrast sign mismatch')
            contrasts.append({'P':p,'Q':q,'m':str(m),'M':str(M),'conservative_min_exact':str(lower),'conservative_max_exact':str(upper),'sign':sign})
            contrast_count+=1
        require(seen==expected_keys,'Incomplete contrast family')
        checked[workflow]['contrasts']=contrasts
    require(box_count==16 and contrast_count==24,'Incomplete results')
    return {'status':'PASS','method':'Exact Fraction residual-root/Dinkelbach updates, independent of threshold enumeration',
            'scope':'Conditional pooled-algebra replay only; no raw record/source/provenance validation, ABM import, simulation, new seeds or bootstrap',
            'profiles_checked':8,'box_extrema_pairs_checked':box_count,'conservative_contrasts_checked':contrast_count,'replayed':checked}


def main():
    here=Path(__file__).resolve().parent
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--pooled',type=Path,default=here/'pooled_worker_coefficients.json')
    parser.add_argument('--results',type=Path,default=here/'results.json')
    parser.add_argument('--output',type=Path,required=True,help='New receipt path; existing files are never overwritten')
    args=parser.parse_args()
    require(not args.output.exists(),'Output already exists; choose a new path')
    protected={args.pooled.resolve(),args.results.resolve(),Path(__file__).resolve()}
    require(args.output.resolve() not in protected,'Output must not replace an input or script')
    pooled,ph=read(args.pooled);results,rh=read(args.results)
    receipt=replay(pooled,results)
    receipt['inputs']={'pooled':{'name':args.pooled.name,'sha256':ph},'results':{'name':args.results.name,'sha256':rh}}
    receipt['script_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    with args.output.open('x') as stream:
        json.dump(receipt,stream,indent=2);stream.write('\n')
    print(f"PASS: 8 profiles, 16 boxes, 24 contrasts; receipt: {args.output}")


if __name__=='__main__':
    main()

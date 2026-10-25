"""Compare cached independent arrays and metrics to the completed primary summary."""
from pathlib import Path
import json
import hashlib
import numpy as np

HERE=Path(__file__).resolve().parent
STUDY=HERE.parents[1]
z=np.load(HERE/'independent_arrays.npz')
K,F,BF=z['K'],z['F'],z['bootstrap_F']
ids=z['ids'].tolist(); loc={p:i for i,p in enumerate(ids)}
ind=json.loads((HERE/'independent_results.json').read_text())
path=STUDY/'summary.json'; s=json.loads(path.read_text())
errors={}; exact={}

def error(name,a,b):
    a,b=np.asarray(a),np.asarray(b)
    if a.shape!=b.shape:raise ValueError((name,a.shape,b.shape))
    errors[name]=max(errors.get(name,0),float(np.max(np.abs(a-b))))

exact['summary_complete']=s['status']=='complete'
exact['summary_self_checks_pass']=s['self_checks_pass']
exact['256_seeds_and_182_pairs']=len(z['seeds'])==s['n_passing']==256 and len(ids)==len(s['pairs'])==182
exact['all_record_hashes']=ind['inputs']['records']==s['inputs']['records_sha256']
exact['bootstrap_index_hash']=ind['bootstrap_index_sha256']==s['bootstrap']['all']['index_matrix_sha256_int64LE']
for pid,j in loc.items():
    p=s['pairs'][pid]
    error('all_182_F',F[j],p['limiting_response'])
    error('all_182_bootstrap_endpoints',np.quantile(BF[:,j],[.025,.975]),p['F_bootstrap_95']['interval'])
    assert p['F_bootstrap_95']['undefined_draws']==0
    x,y=K[:,j,0],K[:,j,1]; residual=(y-(1-F[j])*x)/np.mean(x)
    se=np.std(residual,ddof=1)/np.sqrt(len(x))
    error('all_182_delta_se',se,p['limiting_response_delta_monte_carlo_se'])
    error('all_182_normal_endpoints',F[j]+np.array([-1.96,1.96])*se,p['limiting_response_normal_95_interval'])
    error('all_182_K_means',[np.mean(x),np.mean(y)],[p['K0']['mean'],p['K1']['mean']])

pred=s['predictions']
p1=pred['prediction_1_crossing'];i1=ind['prediction_1']
error('crossing_root',i1['roots'][0],p1['e_star'])
error('crossing_bootstrap_endpoints',i1['root_bootstrap_95']['interval'],p1['e_star_bootstrap_95']['interval'])
exact['crossing_draws_all_one']=p1['bootstrap_sign_change_counts']['draws_one_sign_change']==20000 and ind['prediction_1']['bootstrap_number_of_roots_histogram']=={'1':20000}
for e,d in zip(i1['e_grid'],i1['D']):error('crossing_grid_D',d,p1['D_by_e'][f'{e:g}']['D'])

def compare_diff(prefix,i,p):
    error(prefix+'_D',i['D'],p['D'])
    error(prefix+'_F',[i['F_a'],i['F_b']],[p['F_a'],p['F_b']])
    error(prefix+'_delta_se',i['se'],p['delta_monte_carlo_se'])
    error(prefix+'_normal_endpoints',i['normal_95'],p['normal_95_interval'])
    error(prefix+'_bootstrap_endpoints',i['bootstrap_95']['interval'],p['bootstrap_95']['interval'])
    assert p['bootstrap_95']['undefined_draws']==i['bootstrap_95']['undefined']==0

for key,i in ind['prediction_2'].items():
    span,e=key.split(';e=');h0,h1=span.split('->')
    compare_diff('completion',i,pred['prediction_2_completion_ranges']['D_rows'][f'h{h0}_to_{h1}_e{e}'])
for key,i in ind['prediction_3'].items():
    compare_diff('transition',i,pred['prediction_3_transition_policy']['D'][f'e{float(key):g}'])
p4=pred['prediction_4_renewal_comparator']
for pid,i in ind['prediction_4']['rows'].items():
    error('renewal_series_vs_matrix',i['F_geo'],p4['rows'][pid]['F_geo'])
error('renewal_max_error_pp',ind['prediction_4']['max_abs_error_pp'],p4['max_abs_error_pp'])
exact['renewal_worst_row']=ind['prediction_4']['worst_row']==p4['worst_row']
for truth,i in ind['projection'].items():
    p=s['secondary_exact_icu_projection']['by_truth'][truth]
    error('projection_width_pp',i['median_width_pp'],p['median_width_pp'])
    exact['projection_counts_'+truth]=i['counts']=={'defined':p['defined'],'empty':p['empty'],'out':p['out_of_domain']['total']}
    exact['analytic_domain_'+truth]=i['domain_classification_disagreements']==0
exact['prediction_pass_fail']=list(v['holds'] for v in pred.values())==[True,False,True,False]

# Floating precision tolerances distinguish coefficients' larger physical scale
# from ratios/contrasts. No criterion is selected from observed differences.
tolerances={key:(1e-9 if key=='all_182_K_means' else 1e-10 if key.endswith('_pp') else 1e-12) for key in errors}
accepted=all(exact.values()) and all(errors[key]<=tolerances[key] for key in errors)
out={'accepted':accepted,'exact_checks':exact,'maximum_absolute_errors':errors,'absolute_tolerances':tolerances,
     'summary_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
     'independent_results_sha256':hashlib.sha256((HERE/'independent_results.json').read_bytes()).hexdigest(),
     'scope':'Numerical inference and specified predictions conditional on valid records; not clinical validation, finite-p validation, or certification of continuous ICU interpolation.'}
(HERE/'comparison.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
assert accepted

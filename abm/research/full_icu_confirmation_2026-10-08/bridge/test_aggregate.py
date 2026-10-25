"""Bounded synthetic 38-shard fixtures; never read the running bridge outcomes."""
from copy import deepcopy
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import tempfile
import unittest

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('bridge_analysis',HERE/'aggregate.py')
analysis=importlib.util.module_from_spec(spec);spec.loader.exec_module(analysis)
reviewed=analysis.backend()


def save(path,value):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value,indent=2,sort_keys=True,allow_nan=False)+'\n')


def arm(k,policy=False):
    return {'K':k,'K_from_warmup_pickups':.1*k,'K_by_worker':{'-1':.7*k,'-2':.3*k},
        'K_by_opportunities':{'0':.1*k,'1':0.,'2':.6*k,'>=3':.3*k},
        'raw_pairs_without_reset':100,
        'raw_pairs_without_reset_by_opportunities':{'0':5,'1':0,'2':65,'>=3':30},
        'event_counts':{'imported_carrier_contacts':10,'susceptible_contacts':80,
                        'opportunities_warmup':200,'opportunities_observed':70,'implicit_resets':20},
        'completion_law_endpoint':{'-1':[0.,.05 if policy else .1,.95 if policy else .9]}}


def row(seed,config,identity):
    x=1+seed%5
    profiles={}
    for i,p in enumerate(config['profiles']):
        k0=2*(i+1)*x;k1=k0*(.2+.2*i+.01*(seed%3))
        profiles[p['id']]={'baseline':arm(k0),'policy':arm(k1,True)}
    first=config['profiles'][0]
    return {'seed':seed,'status':'complete','identity_sha256':identity,'elapsed_seconds':1.,
        'model_kwargs':dict(reviewed.FIXED_MODEL,**config['model_overrides'],seed=seed,
           hygiene_opportunities=config['opportunities'],hygiene_autocorrelation=first['rho'],hygiene_efficacy=first['e']),
        'executed_ticks':config['warmup_ticks']+config['observation_ticks'],
        'observation_patient_days':1560.,'observation_susceptible_patient_days':1500.,
        'observation_admissions':120,'observation_imported_colonizations':4,
        'total_admissions':230,'total_imported_colonizations':7,
        'event_counts':{'patient':90,'environment':100,'opportunity':270,'reset':20},
        'ordered_bridge_events_sha256':'a'*64,'profiles':profiles}


def fixture(root):
    real=HERE.parents[2];prefix=f'research/{HERE.parent.name}'
    book=json.loads((HERE/'seedbook.json').read_text())
    names=[f'{prefix}/array_worker.py',f'{prefix}/seed_audit.json',f'{prefix}/seedbook.json',f'{prefix}/runtime_validation.json',
       f'{prefix}/bridge/seedbook.json',f'{prefix}/bridge/PROTOCOL.md',
       *[f'{prefix}/bridge/{n}' for n in ('array_worker.py','run.py','coefficient.py','background.py','test_bridge.py')],
       *[f'src/{n}' for n in ('CPEmodel.py','agents.py','routes.py')],
       *[str(p.relative_to(real)) for p in (HERE/'configs').glob('*.json')]]
    for name in names:
        p=root/name;p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes((real/name).read_bytes())
    freeze={'schema_version':1,'study':HERE.parent.name,'phase':'bridge',
       'files':{name:analysis.digest(root/name) for name in names}}
    freeze['identity_sha256']=hashlib.sha256(analysis.canonical(freeze).encode()).hexdigest()
    freeze_path=root/f'{prefix}/test_freeze_bridge.json';save(freeze_path,freeze)
    freeze_sha=analysis.digest(freeze_path)
    sources={n:freeze['files'][f'{prefix}/bridge/{n}'] for n in ('run.py','coefficient.py','background.py','test_bridge.py')}
    sources.update({f'src/{n}':freeze['files'][f'src/{n}'] for n in ('CPEmodel.py','agents.py','routes.py')})
    scientific={f'{prefix}/bridge/{n}':sources[n] for n in ('run.py','coefficient.py','background.py')}
    scientific.update({f'src/{n}':sources[f'src/{n}'] for n in ('CPEmodel.py','agents.py','routes.py')})
    runtimes=json.loads((root/f'{prefix}/runtime_validation.json').read_text())
    for site,parts in book['partitions'].items():
        for part in parts:
            directory=root/part['output_scope'];raw=json.loads((root/part['config_path']).read_text());cfg=reviewed.normalize(raw)
            manifest={'schema_version':1,'study':HERE.parent.name+'_bridge','input_config':raw,
                      'resolved_config':cfg,'fixed_model':reviewed.FIXED_MODEL,'source_sha256':sources,
                      'versions':runtimes[f'{site}_runtime']}
            manifest['identity_sha256']=hashlib.sha256(analysis.canonical(manifest).encode()).hexdigest()
            save(directory/'manifest.json',manifest)
            outer={'schema_version':1,'study':HERE.parent.name,'module':'bridge',
                'freeze_identity':freeze['identity_sha256'],'freeze_sha256':freeze_sha,
                'site':site,'partition_index':part['index'],'partition_count':len(parts),
                'config_path':part['config_path'],'config_sha256':analysis.digest(root/part['config_path']),
                'protocol_path':f'{prefix}/bridge/PROTOCOL.md','protocol_sha256':analysis.digest(root/f'{prefix}/bridge/PROTOCOL.md'),
                'seed_audit_path':f'{prefix}/seed_audit.json','seed_audit_sha256':analysis.digest(root/f'{prefix}/seed_audit.json'),
                'bridge_seedbook_path':f'{prefix}/bridge/seedbook.json','bridge_seedbook_sha256':analysis.digest(root/f'{prefix}/bridge/seedbook.json'),
                'seeds':part['seeds'],'output_scope':part['output_scope'],'backend':f'{prefix}/bridge/run.py',
                'scientific_sources':scientific,'dispatch_config_path':f'{prefix}/bridge/configs/bridge_{site}.json',
                'dispatch_config_sha256':analysis.digest(root/f'{prefix}/bridge/configs/bridge_{site}.json'),
                'dispatch_output_scope':f'results/{HERE.parent.name}/bridge_{site}'}
            save(directory/'bridge_binding.json',outer)
            records=[row(s,cfg,manifest['identity_sha256']) for s in part['seeds']]
            for r in records:save(directory/'completed'/f"{r['seed']}.json",r)
            summary=reviewed.summarize(records,cfg);summary.update(identity_sha256=manifest['identity_sha256'],status='complete')
            save(directory/'summary.json',summary)
    return book,freeze_path


class RowTests(unittest.TestCase):
    def setUp(self):
        self.cfg=reviewed.normalize(json.loads((HERE/'configs/bridge_all.json').read_text()))
        self.record=row(self.cfg['seeds'][0],self.cfg,'b'*64)

    def test_histogram_worker_state_and_exposure_corruption_rejected(self):
        mutations=[lambda r:r['profiles']['iid_perfect']['baseline']['K_by_opportunities'].update({'2':999.}),
                   lambda r:r['profiles']['iid_perfect']['baseline']['K_by_worker'].update({'-1':999.}),
                   lambda r:r['profiles']['iid_perfect']['baseline']['completion_law_endpoint'].update({'-1':[0.,.1,.1]}),
                   lambda r:r.update(observation_susceptible_patient_days=2000.),
                   lambda r:r.update(observation_patient_days=1559.0001),
                   lambda r:r['model_kwargs'].update(prob_patient_sick=.05)]
        for mutate in mutations:
            r=deepcopy(self.record);mutate(r)
            with self.assertRaises(ValueError):analysis.validate_row(r,r['seed'],self.cfg,'b'*64,reviewed.FIXED_MODEL)

    def test_zero_coefficient_is_retained(self):
        r=deepcopy(self.record)
        for arms in r['profiles'].values():
            for name in ('baseline','policy'):arms[name]=arm(0.,name=='policy')
        analysis.validate_row(r,r['seed'],self.cfg,'b'*64,reviewed.FIXED_MODEL)
        answer=reviewed.summarize([r],dict(self.cfg,seeds=[r['seed']]))
        self.assertIsNone(answer['profiles']['iid_perfect']['limiting_response'])


class CompleteTests(unittest.TestCase):
    def test_all38_exact256_pool_banks_receipts_and_output_binding(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp).resolve();book,freeze=fixture(root)
            rows,cfg,actual_book,provenance=analysis.collect(root,freeze)
            self.assertEqual(len(rows),256);self.assertEqual(len(provenance['banks']),38)
            self.assertEqual(sum(b['jobs'] for b in provenance['banks']),256)
            answer=analysis.summarize_complete(rows,cfg,actual_book,provenance,'c'*64)
            direct=reviewed.summarize(rows,cfg)
            self.assertEqual(answer['profiles'],direct['profiles'])
            self.assertEqual(answer['seeds'],book['all_seeds'])
            self.assertEqual(answer['statistical_banks']['A']['seeds_completed'],128)
            for key in answer['profiles']:
                a=answer['statistical_banks']['A']['profiles'][key];b=answer['statistical_banks']['B']['profiles'][key]
                diff=answer['bank_B_minus_A'][key]['limiting_response']
                self.assertAlmostEqual(diff['estimate'],b['limiting_response']-a['limiting_response'])
                self.assertAlmostEqual(diff['monte_carlo_se'],math.hypot(a['limiting_response_delta_monte_carlo_se'],b['limiting_response_delta_monte_carlo_se']))
            output=root/'synthetic_report';analysis.publish(answer,output)
            binding=json.loads((output/'binding.json').read_text())
            self.assertEqual(binding['summary_sha256'],analysis.digest(output/'summary.json'))
            self.assertEqual(binding['seeds'],book['all_seeds'])
            with self.assertRaises(FileExistsError):analysis.publish(answer,output)

    def test_incomplete_gate_precedes_even_corrupt_outcome_read(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp).resolve();book,freeze=fixture(root)
            first=book['partitions']['local'][0];last=book['partitions']['cluster'][-1]
            (root/first['output_scope']/'completed'/f"{first['seeds'][0]}.json").write_text('not valid json')
            (root/last['output_scope']/'summary.json').unlink()
            with self.assertRaisesRegex(ValueError,'Missing regular artifact'):analysis.collect(root,freeze)

    def test_corrupt_row_outer_runtime_and_shard_summary_are_rejected(self):
        changes=[('row',lambda p,d:d['profiles']['iid_perfect']['policy'].update(K=999.)),
                 ('outer',lambda p,d:d.update(freeze_identity='wrong')),
                 ('manifest',lambda p,d:d['versions'].update(numpy='different')),
                 ('summary',lambda p,d:d['profiles']['iid_perfect']['K0'].update(mean=999.))]
        for kind,mutate in changes:
            with self.subTest(kind=kind),tempfile.TemporaryDirectory() as tmp:
                root=Path(tmp).resolve();book,freeze=fixture(root);part=book['partitions']['local'][0];directory=root/part['output_scope']
                filename={'row':f"completed/{part['seeds'][0]}.json",'outer':'bridge_binding.json','manifest':'manifest.json','summary':'summary.json'}[kind]
                p=directory/filename;data=json.loads(p.read_text());mutate(p,data)
                if kind=='manifest':
                    payload={k:v for k,v in data.items() if k!='identity_sha256'}
                    data['identity_sha256']=hashlib.sha256(analysis.canonical(payload).encode()).hexdigest()
                save(p,data)
                with self.assertRaises(ValueError):analysis.collect(root,freeze)


if __name__=='__main__':unittest.main()

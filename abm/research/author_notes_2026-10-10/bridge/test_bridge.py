"""Tiny independent history enumeration, original replay and publication guards."""
import importlib.util
import json
import math
from pathlib import Path
import sys
import tempfile
import time
import unittest

try:
    from .coefficient import BridgeCoefficient, RawBridgeHistogram
    from .background import ABM, FIXED_MODEL, make_background_model
    from .run import PAIRS, SCIENCE, atomic_json, normalize, simulate_seed, validate_replay, seal_seed, verify_seed_payload
except ImportError:
    from coefficient import BridgeCoefficient, RawBridgeHistogram
    from background import ABM, FIXED_MODEL, make_background_model
    from run import PAIRS, SCIENCE, atomic_json, normalize, simulate_seed, validate_replay, seal_seed, verify_seed_payload


def brute_force(events,h0,h1,rho,e):
    # Scalar branch enumeration, independent of matrix/vector recurrence.
    branches=[(1.,None,[],0.)]
    for kind,observed,details in events:
        next_branches=[]
        for prob,previous,sources,reward in branches:
            if kind=='opportunity':
                h=h1 if observed else h0
                completion=h if previous is None else (h+rho*(1-h) if previous else h*(1-rho))
                for yes,chance in ((False,1-completion),(True,completion)):
                    for clear,effchance in ((False,1-e),(True,e)) if yes else ((False,1.),):
                        next_branches.append((prob*chance*effchance,yes,[] if clear else list(sources),reward))
            elif kind=='reset':
                next_branches.append((prob,previous,[],reward))
            elif details[0]:
                next_branches.append((prob,previous,sources+[details[1]],reward))
            else:
                next_branches.append((prob,previous,list(sources),reward+(sum(sources)*details[1] if observed else 0.)))
        branches=next_branches
    assert math.isclose(sum(p for p,_,_,_ in branches),1.,abs_tol=1e-11)
    return sum(p*k for p,_,_,k in branches)


class IndependentCoefficientTests(unittest.TestCase):
    def test_all_22_pairs_against_brute_force_both_arms(self):
        self.assertEqual(len(PAIRS),22)
        self.assertEqual(len({p['id'] for p in PAIRS}),22)
        self.assertEqual(sum('legacy_id' in p for p in PAIRS),3)
        events=[('patient',False,(True,.6)),('opportunity',False,()),
                ('patient',False,(True,1.)),('opportunity',False,()),
                ('patient',True,(False,.6)),('opportunity',True,()),
                ('patient',True,(False,1.)),('reset',True,()),
                ('patient',True,(True,.6)),('opportunity',True,()),
                ('patient',True,(False,1.))]
        for pair in PAIRS:
            for h1 in (pair['h0'],pair['h1']):
                with self.subTest(pair=pair['id'],h1=h1):
                    arm=BridgeCoefficient(pair['h0'],h1,pair['rho'],pair['e'],None)
                    for kind,observed,details in events:
                        if kind=='opportunity':arm.opportunity('worker',observed)
                        elif kind=='reset':arm.reset('worker')
                        else:arm.patient_contact('worker','admission',*details,observed)
                    expected=brute_force(events,pair['h0'],h1,pair['rho'],pair['e'])
                    self.assertAlmostEqual(arm.result()['K'],expected,places=12)

    def test_raw_histogram_fine_bins_tail_and_worker_resets(self):
        histogram=RawBridgeHistogram(12)
        histogram.patient_contact('A',True,False)
        histogram.patient_contact('A',False,False)  # warmup targets not counted
        histogram.patient_contact('B',False,True)   # cannot join workers
        for n in range(15):
            histogram.patient_contact('A',False,True)
            histogram.opportunity('A')
        histogram.reset('A')
        histogram.patient_contact('A',False,True)
        result=histogram.result()
        self.assertEqual(result['raw_pairs_without_reset'],15)
        self.assertEqual(result['raw_pairs_without_reset_by_opportunities'],
                         dict([(str(n),1) for n in range(13)]+[('>=13',2)]))


class GuardTests(unittest.TestCase):
    def test_frozen_science_and_unique_seed_bank(self):
        seeds=json.loads((Path(__file__).parent/'config_all.json').read_text())['seeds']
        self.assertEqual(len(seeds),256)
        self.assertEqual(normalize(dict(seeds=seeds))['science'],SCIENCE)
        for raw in (dict(seeds=[seeds[0],seeds[0]]),dict(seeds=[10801]),
                    dict(seeds=[seeds[0]],warmup_ticks=1),dict(seeds=[seeds[0]],h0=.8)):
            with self.assertRaises(ValueError):normalize(raw)

    def test_completed_payload_digest_detects_nonlegacy_arm_drift(self):
        row=seal_seed(dict(seed=42,pairs={'new_arm':{'K':1.}}))
        verify_seed_payload(row)
        row['pairs']['new_arm']['K']=2.
        with self.assertRaises(ValueError):verify_seed_payload(row)

    def test_atomic_immutable_json_refuses_conflict(self):
        with tempfile.TemporaryDirectory() as temporary:
            path=Path(temporary)/'seed.json'
            atomic_json(path,dict(status='complete',K=1.),immutable=True)
            atomic_json(path,dict(status='complete',K=1.),immutable=True)
            with self.assertRaises(ValueError):atomic_json(path,dict(status='complete',K=2.),immutable=True)
            self.assertEqual(json.loads(path.read_text())['K'],1.)
            self.assertFalse(path.with_suffix('.json.tmp').exists())


@unittest.skipUnless(importlib.util.find_spec('mesa'),'Mesa unavailable')
class TinyReplayTests(unittest.TestCase):
    def test_old_adapter_same_seed_and_role_observation_counters(self):
        old=ABM/'research/full_icu_confirmation_2026-10-08/bridge'
        # Import original coefficient in isolation so this check executes frozen
        # old code, including the old coarse histogram, as its own dependency.
        def load(name,path):
            spec=importlib.util.spec_from_file_location(name,path)
            module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module
        previous=sys.modules.get('coefficient')
        sys.modules['coefficient']=load('_author_test_old_coefficient',old/'coefficient.py')
        try:original=load('_author_test_old_background',old/'background.py')
        finally:
            if previous is None:sys.modules.pop('coefficient',None)
            else:sys.modules['coefficient']=previous
        old_factory=original.make_background_model()
        new_factory=make_background_model()
        profiles=[dict(id=p['legacy_id'],rho=p['rho'],e=p['e']) for p in PAIRS if p.get('legacy_id')]
        kwargs=dict(FIXED_MODEL,seed=10801,prob_patient_sick=.03,hygiene_opportunities='zone_separate',
                    hygiene_autocorrelation=0.,hygiene_efficacy=1.)
        old_model=old_factory(warmup_ticks=128,h0=.9,h1=.95,profiles=profiles,histogram_cap=2,**kwargs)
        new_model=new_factory(warmup_ticks=128,pairs=PAIRS,**kwargs)
        fields=('hygiene_opportunity_count','patient_contacts','environment_contacts')
        def worker_counts(model):
            return {a.unique_id:(a.hygiene_role,*[getattr(a,f) for f in fields])
                    for a in model.schedule.agents if a.isHCW}
        before=None
        for tick in range(1024):
            old_model.step();new_model.step()
            if tick==127:before=worker_counts(new_model)
        old_result=old_model.bridge_result();new_result=new_model.bridge_result()
        self.assertEqual(old_result['ordered_bridge_events_sha256'],new_result['ordered_bridge_events_sha256'])
        self.assertEqual(old_result['event_counts'],new_result['event_counts'])
        self.assertEqual(json.dumps(old_model.rng.bit_generator.state),json.dumps(new_model.rng.bit_generator.state))
        self.assertEqual(old_model.random.getstate(),new_model.random.getstate())
        for pair in PAIRS:
            if pair.get('legacy_id'):
                self.assertEqual(old_result['profiles'][pair['legacy_id']],new_result['pairs'][pair['id']])
        after=worker_counts(new_model)
        for role in ('Nurse','Dr','XrayDr'):
            for index,label in enumerate(('opportunity','patient','surface'),start=1):
                expected=sum(after[uid][index]-before[uid][index] for uid in after if after[uid][0]==role)
                self.assertEqual(new_result['observation_event_counts_by_role'][role][label],expected)
        old_raw=old_result['profiles']['iid_perfect']['baseline']['raw_pairs_without_reset_by_opportunities']
        new_raw=new_result['raw_pairs_without_reset_by_opportunities']
        self.assertEqual([new_raw[str(i)] for i in range(3)],[old_raw[str(i)] for i in range(3)])
        self.assertEqual(sum(new_raw[str(i)] for i in range(3,13))+new_raw['>=13'],old_raw['>=3'])
        for role in ('Nurse','Dr','XrayDr'):
            self.assertGreater(new_result['observation_event_counts_by_role'][role]['opportunity'],0)

    def test_expired_deadline_restarts_without_steps(self):
        with self.assertRaises(TimeoutError):
            simulate_seed(10801,make_background_model(),time.monotonic()-1,warmup_ticks=128,observation_ticks=896)


if __name__=='__main__':unittest.main()

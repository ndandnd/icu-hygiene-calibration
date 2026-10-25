"""Read-only frozen ICU adapter for a clean-reservoir p=0 background.

Import isolation keeps hooks out of canonical modules. Only the loaded private
classes are patched. All frozen p=0 RNG calls are retained, including impossible
pickup trials and imported-carrier screening; hooks consume no random draws.
"""
from collections import Counter
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

try:
    from .coefficient import BridgeCoefficient, RawBridgeHistogram
except ImportError:
    from coefficient import BridgeCoefficient, RawBridgeHistogram

HERE = Path(__file__).resolve().parent
ABM = HERE.parents[2]
CORE_SHA256 = {
    'CPEmodel.py': '5659a207818cdc8f1dcd40dd41b2900fa9d02a302a92b9dc6c5f6b376213ff53',
    'agents.py': 'bb4d52252f2c4d4fba599245fa1abc07254afe12fc2165e692ad1279504fc656',
    'routes.py': '48c1ce8440cf71a82996afd92ede83b1237744db56b439bd9163ce7b60f0dcd5',
}
FIXED_MODEL = dict(num_HCWs=10, num_Patients=30, num_Goo=30,
                   prob_patient_sick=.05, prob_new_patient=.003,
                   isolation_factor=.6, cleaningDay=1., isolate_sick=True,
                   height=11, width=32, stay_extension_days=0,
                   los_mode='synthetic', susceptible_los_mean_days=5.,
                   colonized_los_mean_days=14.,
                   initial_environment_contaminated=False,
                   data_collection_interval=2592)


def load_core(abm_root=ABM, tag='_bridge'):
    """Return private classes. Serial loading only (one process per worker)."""
    src = Path(abm_root) / 'src'
    for name, expected in CORE_SHA256.items():
        if hashlib.sha256((src / name).read_bytes()).hexdigest() != expected:
            raise RuntimeError('Frozen scientific source changed: ' + name)

    def load(name, path):
        spec = importlib.util.spec_from_file_location(name, path)
        module = importlib.util.module_from_spec(spec)
        sys.modules[name] = module
        spec.loader.exec_module(module)
        return module

    agents = load(tag + '_agents', src / 'agents.py')
    routes = load(tag + '_routes', src / 'routes.py')
    old = {name: sys.modules.get(name) for name in ('agents', 'routes')}
    sys.modules.update(agents=agents, routes=routes)
    try:
        core = load(tag + '_core', src / 'CPEmodel.py')
    finally:
        for name, module in old.items():
            if module is None:
                sys.modules.pop(name, None)
            else:
                sys.modules[name] = module
    return agents, core


def make_background_model(abm_root=ABM):
    agents, core = load_core(abm_root)
    original_wash = agents.HCW.handwash

    def wash(worker):
        worker.model.bridge_event('opportunity', worker.unique_id)
        original_wash(worker)

    def spread(worker):
        model = worker.model
        if model.hygiene_opportunities == 'zone_separate':
            worker._update_zone_hygiene()
        worker.reached = False
        eligible = [a for a in model.grid.get_cell_list_contents([worker.pos]) if a.isPatient or a.isGoo]
        if worker._contact_position != worker.pos:
            worker._contacted_objects.clear()
        worker._contact_position = worker.pos
        worker._contacted_objects.intersection_update({a.unique_id for a in eligible})
        for other in eligible:
            if other.unique_id in worker._contacted_objects:
                continue
            worker._contacted_objects.add(other.unique_id)
            worker.reached = True
            counter = 'patient_contacts' if other.isPatient else 'environment_contacts'
            setattr(worker, counter, getattr(worker, counter) + 1)
            setattr(model, counter, getattr(model, counter) + 1)
            worker.last_contact_tick = model.total_ticks
            worker.last_contact_pos = worker.pos
            if model.hygiene_opportunities != 'zone_separate':
                worker.handwash()
            if worker.colonized or (other.isGoo and other.colonized):
                raise RuntimeError('Clean p=0 background became contaminated')
            if other.isPatient:
                weight = model.isolation_factor if other.isolated else 1.0
                model.bridge_event('patient', worker.unique_id, other.unique_id,
                                   bool(other.colonized), weight)
            else:
                model.bridge_event('environment', worker.unique_id, other.unique_id)
            # Frozen p=0 clean-hand control flow still consumes the pickup draw.
            if other.colonized:
                model.rng.random()
            if other.isPatient and other.colonized:
                worker.testCRE(other)
            if model.hygiene_opportunities == 'separate':
                worker.handwash()

    agents.HCW.handwash = wash
    agents.HCW.spread = spread

    def reset_wrapper(original, role):
        def step(worker):
            model = worker.model
            if role == 'Nurse':
                reset = model.tick % (model.ticks_in_hour * worker.workHours) == 0
            elif role == 'Dr':
                reset = model.tick in (worker.startTime, worker.endTime % model.ticks_in_day)
            else:
                # Both Xray resets precede all movement and pending zone washes.
                reset = model.tick % (model.ticks_in_day / worker.shiftsPerDay) == 0
                finishing = (worker.activated and worker._route_entry_completed
                             and worker.moveTick == len(worker.path) - 1)
                # A scheduled restart sets moveTick=0, entry_completed=False,
                # so it suppresses the second condition in the original method.
                reset = reset or finishing
            if reset:
                model.bridge_event('reset', worker.unique_id)
            original(worker)
        return step

    for role in ('Nurse', 'Dr', 'XrayDr'):
        cls = getattr(agents, role)
        cls.step = reset_wrapper(cls.step, role)

    class BackgroundModel(core.CPE_Model):
        def __init__(self, *, warmup_ticks, pairs, **kwargs):
            if kwargs.get('prob_transmission', 0) != 0 or kwargs.get('initial_environment_contaminated', False):
                raise ValueError('Bridge theorem requires p=0 and clean initial environment')
            # Epidemiological replay law is always the ORIGINAL first profile.
            kwargs.update(prob_transmission=0.0, initial_environment_contaminated=False,
                          icu_hcw_wash_rate=.90, outside_hcw_wash_rate=.90,
                          hygiene_autocorrelation=0., hygiene_efficacy=1.)
            self.bridge_warmup_ticks = warmup_ticks
            self.bridge_pairs = pairs
            self.bridge_arms = [BridgeCoefficient(pair['h0'], h, pair['rho'], pair['e'],
                                                 2 if pair.get('legacy_id') else None)
                               for pair in pairs for h in (pair['h0'], pair['h1'])]
            self.bridge_counts = Counter()
            self.observation_counts_by_role = {role: Counter() for role in ('Nurse', 'Dr', 'XrayDr')}
            self.raw_bridges = RawBridgeHistogram(12)
            self.bridge_digest = hashlib.sha256()
            super().__init__(**kwargs)
            self.bridge_roles = {a.unique_id: a.hygiene_role for a in self.schedule.agents if a.isHCW}

        def bridge_event(self, kind, uid, *details):
            observed = self.total_ticks > self.bridge_warmup_ticks
            event = (self.total_ticks, kind, uid, *details)
            # Keep original digest bytes exactly: role accounting consumes no RNG.
            self.bridge_digest.update(json.dumps(event, separators=(',', ':'), allow_nan=False).encode() + b'\n')
            self.bridge_counts[kind] += 1
            if observed:
                role = self.bridge_roles[uid]
                label = 'surface' if kind == 'environment' else kind
                self.observation_counts_by_role[role][label] += 1
            if kind == 'opportunity':
                self.raw_bridges.opportunity(uid)
                for arm in self.bridge_arms:
                    arm.opportunity(uid, observed)
            elif kind == 'reset':
                self.raw_bridges.reset(uid)
                for arm in self.bridge_arms:
                    arm.reset(uid)
            elif kind == 'patient':
                patient_id, carrier, weight = details
                self.raw_bridges.patient_contact(uid, carrier, observed)
                for arm in self.bridge_arms:
                    arm.patient_contact(uid, patient_id, carrier, weight, observed)

        def bridge_result(self):
            self.assert_invariants()
            if self.acquired_colonizations or any(a.colonized for a in self.schedule.agents if a.isHCW or a.isGoo):
                raise RuntimeError('p=0 reservoir invariant failed')
            result = {'pairs': {pair['id']: {'baseline': self.bridge_arms[2*i].result(),
                                             'policy': self.bridge_arms[2*i+1].result()}
                               for i, pair in enumerate(self.bridge_pairs)},
                      'event_counts': dict(self.bridge_counts),
                      'observation_event_counts_by_role': {
                          role: {kind: counts[kind] for kind in ('opportunity', 'patient', 'surface', 'reset')}
                          for role, counts in self.observation_counts_by_role.items()},
                      'ordered_bridge_events_sha256': self.bridge_digest.hexdigest()}
            result.update(self.raw_bridges.result())
            for pair in self.bridge_pairs:
                for arm in result['pairs'][pair['id']].values():
                    if arm['raw_pairs_without_reset'] != result['raw_pairs_without_reset']:
                        raise RuntimeError('Raw bridge arm reconciliation failed')
            return result

    return BackgroundModel

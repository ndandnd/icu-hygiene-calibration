"""Passive proximate-lineage observations over isolated frozen ICU classes.

No extra random draws. Environment multiplier !=1 is an explicit intervention.
Origin tags describe the last successful pickup, not causal attribution.
"""
from collections import Counter
import hashlib
import importlib.util
from pathlib import Path
import sys
HERE = Path(__file__).resolve().parent
ABM = HERE.parents[1]
SRC = ABM / 'src'
CORE_SHA256 = {'CPEmodel.py': '5659a207818cdc8f1dcd40dd41b2900fa9d02a302a92b9dc6c5f6b376213ff53', 'agents.py': 'bb4d52252f2c4d4fba599245fa1abc07254afe12fc2165e692ad1279504fc656', 'routes.py': '48c1ce8440cf71a82996afd92ede83b1237744db56b439bd9163ce7b60f0dcd5'}

def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module
for _name, _hash in CORE_SHA256.items():
    if hashlib.sha256((SRC / _name).read_bytes()).hexdigest() != _hash:
        raise RuntimeError('Frozen scientific source changed: ' + _name)
agents = _load('_response_agents', SRC / 'agents.py')
_previous = sys.modules.get('agents')
sys.modules['agents'] = agents
sys.path.insert(0, str(SRC))
try:
    core = _load('_response_core', SRC / 'CPEmodel.py')
finally:
    sys.path.pop(0)
    if _previous is None:
        sys.modules.pop('agents', None)
    else:
        sys.modules['agents'] = _previous
_original_handwash = agents.HCW.handwash
_original_clean = agents.Goo.handwash

# Exact frozen spread control flow with deterministic success callbacks.
# The core hash guard refuses future revisions rather than silently drifting.
def spread(self):
    self.model._observe_implicit_reset(self)
    if self.model.hygiene_opportunities == 'zone_separate':
        self._update_zone_hygiene()
    self.reached = False
    eligible = [a for a in self.model.grid.get_cell_list_contents([self.pos]) if a.isPatient or a.isGoo]
    if self._contact_position != self.pos:
        self._contacted_objects.clear()
    self._contact_position = self.pos
    present = {a.unique_id for a in eligible}
    self._contacted_objects.intersection_update(present)
    for other in eligible:
        if other.unique_id in self._contacted_objects:
            continue
        self._contacted_objects.add(other.unique_id)
        self.reached = True
        counter = 'patient_contacts' if other.isPatient else 'environment_contacts'
        setattr(self, counter, getattr(self, counter) + 1)
        setattr(self.model, counter, getattr(self.model, counter) + 1)
        self.last_contact_tick = self.model.total_ticks
        self.last_contact_pos = self.pos
        if self.model.hygiene_opportunities != 'zone_separate':
            self.handwash()
        prob = self.model.prob_transmission
        if other.isolated:
            prob *= self.model.isolation_factor
        if other.isGoo:
            prob *= self.model.environment_multiplier
        if self.colonized and (not other.colonized):
            if self.model.rng.random() < prob:
                other.colonized = True
                if other.isPatient:
                    other.state = 'colonized'
                    other.stay += self.model.stay_extension_ticks
                    self.model.cumul_sick_patients += 1
                    self.model.num_infecByHCW += 1
                    self.model.acquired_colonizations += 1
                    self.model._observe_acquisition(self, other)
                else:
                    self.model._observe_deposition(self, other)
        elif not self.colonized and other.colonized:
            if self.model.rng.random() < prob:
                self.colonized = True
                self.model._observe_pickup(self, other)
        if other.isPatient and other.colonized:
            self.testCRE(other)
        if self.model.hygiene_opportunities == 'separate':
            self.handwash()

def handwash(self):
    model = self.model
    origin = getattr(self, '_response_origin', None)
    if not self.colonized:
        origin = None
        self._response_origin = None
    _original_handwash(self)
    if origin is not None:
        n = self.hygiene_opportunity_count - origin['opportunities_at_pickup']
        model.response_counts['hand_at_risk_opportunities'][str(n)] += 1
        key = 'hand_survived_opportunities' if self.colonized else 'hand_cleared_opportunities'
        model.response_counts[key][str(n)] += 1
        if not self.colonized:
            self._response_origin = None

def clean(self):
    was_dirty = self.colonized
    _original_clean(self)
    if hasattr(self.model, 'response_environment_count'):
        if was_dirty:
            self.model.response_environment_count -= 1
        self._response_environment_origin = None
        self.model.response_counts['environment_events']['scheduled_clear'] += 1
        if was_dirty:
            self.model.response_counts['environment_events']['dirty_scheduled_clear'] += 1
agents.HCW.spread = spread
agents.HCW.handwash = handwash
agents.Goo.handwash = clean

class InstrumentedModel(core.CPE_Model):

    def __init__(self, *args, environment_multiplier=1.0, **kwargs):
        if isinstance(environment_multiplier, bool) or not isinstance(environment_multiplier, (int, float)) or (not 0 <= environment_multiplier <= 1):
            raise ValueError('environment_multiplier must be finite and in [0,1]')
        self.environment_multiplier = float(environment_multiplier)
        self.response_counts = {k: Counter() for k in ('pickups', 'pickups_by_completion_state', 'acquisitions_by_patient_source_status', 'acquisitions_by_route', 'acquisitions_by_opportunities', 'acquisitions_by_route_and_opportunities', 'hand_at_risk_opportunities', 'hand_survived_opportunities', 'hand_cleared_opportunities', 'hand_implicit_reset_opportunities', 'environment_events')}
        self.response_admission_generation = Counter()
        self.response_last_admission = {}
        super().__init__(*args, **kwargs)
        self.response_sites = [a for a in self.schedule.agents if a.isGoo]
        self.response_workers = [a for a in self.schedule.agents if a.isHCW]
        self.response_environment_count = sum((a.colonized for a in self.response_sites))
        self.response_environment_ticks = 0
        for worker in self.response_workers:
            worker._response_origin = None
        for site in self.response_sites:
            zone = self.care_zone_by_position[site.pos]
            site._response_environment_origin = {'kind': 'initial_unknown', 'zone': zone, 'admission_generation': self.response_admission_generation[zone], 'donor_patient_id': None} if site.colonized else None

    def admit_patient(self, bed):
        patient = super().admit_patient(bed)
        zone = bed.unique_id
        self.response_admission_generation[zone] += 1
        self.response_last_admission[zone] = patient.unique_id
        patient._response_patient_origin = 'imported' if patient.colonized else 'susceptible'
        return patient

    def _observe_pickup(self, worker, source):
        zone = self.care_zone_by_position.get(source.pos)
        lineage = None if source.isPatient else getattr(source, '_response_environment_origin', None)
        worker._response_origin = {'kind': 'patient' if source.isPatient else 'environment', 'source_id': source.unique_id, 'zone': zone, 'tick': self.total_ticks, 'opportunities_at_pickup': worker.hygiene_opportunity_count, 'environment_lineage': None if lineage is None else dict(lineage),
            'completion_at_pickup': worker._hygiene_previous_completion,
            'patient_source_status': getattr(source, '_response_patient_origin', None) if source.isPatient else None}
        self.response_counts['pickups'][worker._response_origin['kind']] += 1
        state = worker._hygiene_previous_completion
        self.response_counts['pickups_by_completion_state']['none' if state is None else ('completed' if state else 'missed')] += 1

    def _observe_deposition(self, worker, site):
        self.response_environment_count += 1
        source = getattr(worker, '_response_origin', None)
        zone = self.care_zone_by_position[site.pos]
        donor = None
        if source is not None:
            if source['kind'] == 'patient':
                donor = source['source_id']
            elif source.get('environment_lineage'):
                donor = source['environment_lineage'].get('donor_patient_id')
        site._response_environment_origin = {'kind': 'hand_deposition', 'zone': zone, 'donor_patient_id': donor, 'admission_generation': self.response_admission_generation[zone]}
        self.response_counts['environment_events']['deposition'] += 1

    def _observe_acquisition(self, worker, patient):
        patient._response_patient_origin = 'acquired'
        origin = getattr(worker, '_response_origin', None)
        zone = self.care_zone_by_position[patient.pos]
        route = 'untracked_hand'
        n = None
        if origin is not None:
            n = worker.hygiene_opportunity_count - origin['opportunities_at_pickup']
            if origin['kind'] == 'patient':
                route = 'patient_pickup'
                self.response_counts['acquisitions_by_patient_source_status'][origin.get('patient_source_status') or 'unknown'] += 1
            elif origin['zone'] != zone:
                route = 'other_zone_environment_pickup'
            else:
                route = 'own_zone_environment_pickup'
                lineage = origin.get('environment_lineage')
                if lineage is not None and lineage['admission_generation'] < self.response_admission_generation[zone] and (self.response_last_admission.get(zone) == patient.unique_id):
                    route = 'own_zone_environment_pickup_after_readmission'
        gap = 'unknown' if n is None else str(n)
        self.response_counts['acquisitions_by_route'][route] += 1
        self.response_counts['acquisitions_by_opportunities'][gap] += 1
        self.response_counts['acquisitions_by_route_and_opportunities'][route + '|' + gap] += 1

    def _observe_implicit_reset(self, worker):
        origin = getattr(worker, '_response_origin', None)
        if not worker.colonized and origin is not None:
            n = worker.hygiene_opportunity_count - origin['opportunities_at_pickup']
            self.response_counts['hand_implicit_reset_opportunities'][str(n)] += 1
            worker._response_origin = None

    def step(self):
        self.response_environment_ticks += self.response_environment_count
        super().step()
        for worker in self.response_workers:
            self._observe_implicit_reset(worker)

    def assert_response_invariants(self):
        self.assert_invariants()
        assert self.response_environment_count == sum((a.colonized for a in self.response_sites))
        assert sum(self.response_counts['acquisitions_by_route'].values()) == self.acquired_colonizations
        assert sum(self.response_counts['acquisitions_by_opportunities'].values()) == self.acquired_colonizations
        assert sum(self.response_counts['acquisitions_by_route_and_opportunities'].values()) == self.acquired_colonizations
        assert self.response_counts['acquisitions_by_route']['untracked_hand'] == 0

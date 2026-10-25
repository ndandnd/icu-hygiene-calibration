"""p=0 bridge background under the original or the alternative ICU workflow.

Study: research/workflow_coefficient_2026-10-21 (see ../PROTOCOL.md).

Composition, without editing or copying any frozen scientific file:

* ``research/full_icu_confirmation_2026-10-08/bridge/background.py`` supplies the
  reviewed p=0 ``BackgroundModel`` (hand/opportunity/reset hooks, event digest,
  exact ``BridgeCoefficient`` arms). It is loaded by path and hash-checked.
* ``research/robustness_uq_2026-10-17/workflow_variant.py`` supplies the
  alternative workflow rules. Its four rule methods
  (``_clear_vacated_site``, ``_remove_discharged``, ``_move_patient``,
  ``_relocate_detected_patients``) are reused as the *same code objects*: the
  frozen class subclasses a different base (``InstrumentedModel``), so the two
  methods that use zero-argument ``super()`` get a fresh function object whose
  ``__class__`` cell names this adapter's class. Bytecode, constants and
  globals are the frozen ones; only the ``super()`` anchor changes. The small
  amount of state those methods expect from ``WorkflowRobustnessModel.__init__``
  and ``InstrumentedModel.__init__`` is re-expressed below with line
  citations.

``alternative=False`` returns a class that adds only passive bookkeeping (extra
coefficient arms fed the frozen event stream, a relocation audit, a
warmup-boundary counter snapshot, a patient-contact tally by carrier status and
bed type). It consumes no random draws and never changes positions, so it must
reproduce the archived original-workflow records exactly (Gate 1).
``alternative=True`` adds the four frozen rule methods on top of that class.
Those rules make no random draws either; they change only where detected
patients lie (and, at p=0, clear sites that are already clean), which can
change later contacts, admissions and therefore the RNG path.

Never writes bytecode: ``sys.dont_write_bytecode`` is set before any frozen
module is executed, so frozen ``__pycache__`` folders are untouched even when a
caller forgets ``python -B``.
"""
import sys

sys.dont_write_bytecode = True

from collections import Counter  # noqa: E402
import hashlib  # noqa: E402
import importlib.util  # noqa: E402
import json  # noqa: E402
from pathlib import Path  # noqa: E402
import types  # noqa: E402

HERE = Path(__file__).resolve().parent
STUDY = HERE.parent
ABM = HERE.parents[2]
BRIDGE = ABM / 'research/full_icu_confirmation_2026-10-08/bridge'
ROBUSTNESS = ABM / 'research/robustness_uq_2026-10-17'

# Every frozen file this adapter executes or relies on, ABM-root relative.
FROZEN_SHA256 = {
    'src/CPEmodel.py': '5659a207818cdc8f1dcd40dd41b2900fa9d02a302a92b9dc6c5f6b376213ff53',
    'src/agents.py': 'bb4d52252f2c4d4fba599245fa1abc07254afe12fc2165e692ad1279504fc656',
    'src/routes.py': '48c1ce8440cf71a82996afd92ede83b1237744db56b439bd9163ce7b60f0dcd5',
    'research/full_icu_confirmation_2026-10-08/bridge/background.py': '4613669a9ffd9fbf749bcc5d27bb1e62844590a49ddfb23eda579b11c3e0ca04',
    'research/full_icu_confirmation_2026-10-08/bridge/coefficient.py': 'da9b21c5c95ea208b08ef577e2305a151ee7173c139b0e17cefcdf20a4c52898',
    'research/full_icu_confirmation_2026-10-08/bridge/run.py': '9a0959a91007425abf9cf270d3be2142209d7e2f92dbe27f5ebc6bbb2025526b',
    'research/full_icu_confirmation_2026-10-08/bridge/seedbook.json': '7f38db102597f653b70f35560a5aeebe1b01231be49b4802bd7bb209991a518c',
    'research/full_icu_confirmation_2026-10-08/bridge/configs/bridge_all.json': '8e9f0f431ded41e2e2d35a60414851798572ce6a3a6784cc65edad14d582015c',
    'research/robustness_uq_2026-10-17/workflow_variant.py': 'a1037d7cfa83889e6da010671cfa54ade6c599a8de4c45f34113744aa693ef74',
    'research/full_icu_response_2026-10-08/instrumented_model.py': '5dcbe23b392694c1baffa7cd2475588ab12d8a3333778c17e73e0b5e6131888b',
}

# Archived extended run behind supplement Table S5 (outside the ABM root).
EXTENDED_RESULTS = (ABM.parent / 'research/author_notes_2026-10-10/collections/'
                    'bridge_20261011_061144_utc/results/bridge')
EXTENDED_MANIFEST = EXTENDED_RESULTS / 'part_00/manifest.json'
EXTENDED_MANIFEST_SHA256 = 'e50eaefb3f8aac31600abcd4fda1cf3b5d357f6e79e701bdcf59560d54bc0ed9'

# The 22 Table S5 comparisons, literally as in
# EXTENDED_MANIFEST['resolved_config']['science']['pairs'] (verified by
# ``verify_surface_pairs``). ``legacy_id`` marks the three primary profiles,
# which are evaluated by the frozen BackgroundModel's own profile arms.
SURFACE_PAIRS = (
    dict(id='h0.9_to_0.95_rho0_e1', h0=.9, h1=.95, rho=0., e=1., legacy_id='iid_perfect'),
    dict(id='h0.9_to_0.95_rho0_e0.9', h0=.9, h1=.95, rho=0., e=.9),
    dict(id='h0.9_to_0.95_rho0_e0.8', h0=.9, h1=.95, rho=0., e=.8, legacy_id='iid_imperfect'),
    dict(id='h0.9_to_0.95_rho0_e0.6', h0=.9, h1=.95, rho=0., e=.6),
    dict(id='h0.9_to_0.95_rho0.25_e1', h0=.9, h1=.95, rho=.25, e=1.),
    dict(id='h0.9_to_0.95_rho0.25_e0.9', h0=.9, h1=.95, rho=.25, e=.9),
    dict(id='h0.9_to_0.95_rho0.25_e0.8', h0=.9, h1=.95, rho=.25, e=.8),
    dict(id='h0.9_to_0.95_rho0.25_e0.6', h0=.9, h1=.95, rho=.25, e=.6),
    dict(id='h0.9_to_0.95_rho0.5_e1', h0=.9, h1=.95, rho=.5, e=1., legacy_id='persistent_perfect'),
    dict(id='h0.9_to_0.95_rho0.5_e0.9', h0=.9, h1=.95, rho=.5, e=.9),
    dict(id='h0.9_to_0.95_rho0.5_e0.8', h0=.9, h1=.95, rho=.5, e=.8),
    dict(id='h0.9_to_0.95_rho0.5_e0.6', h0=.9, h1=.95, rho=.5, e=.6),
    dict(id='h0.9_to_0.95_rho0.75_e1', h0=.9, h1=.95, rho=.75, e=1.),
    dict(id='h0.9_to_0.95_rho0.75_e0.9', h0=.9, h1=.95, rho=.75, e=.9),
    dict(id='h0.9_to_0.95_rho0.75_e0.8', h0=.9, h1=.95, rho=.75, e=.8),
    dict(id='h0.9_to_0.95_rho0.75_e0.6', h0=.9, h1=.95, rho=.75, e=.6),
    dict(id='h0.8_to_0.9_rho0_e1', h0=.8, h1=.9, rho=0., e=1.),
    dict(id='h0.8_to_0.9_rho0.5_e1', h0=.8, h1=.9, rho=.5, e=1.),
    dict(id='h0.8_to_0.9_rho0_e0.8', h0=.8, h1=.9, rho=0., e=.8),
    dict(id='h0.95_to_0.99_rho0_e1', h0=.95, h1=.99, rho=0., e=1.),
    dict(id='h0.95_to_0.99_rho0.5_e1', h0=.95, h1=.99, rho=.5, e=1.),
    dict(id='h0.95_to_0.99_rho0_e0.8', h0=.95, h1=.99, rho=0., e=.8),
)
PRIMARY_PAIR_IDS = {'iid_perfect': 'h0.9_to_0.95_rho0_e1',
                    'persistent_perfect': 'h0.9_to_0.95_rho0.5_e1',
                    'iid_imperfect': 'h0.9_to_0.95_rho0_e0.8'}
# Histogram cap for the 19 non-primary arms. The archived extended run used no
# histogram for them (cap None); K is computed independently of the cap.
SURFACE_HISTOGRAM_CAP = None


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def verify_frozen(abm_root=ABM):
    """Return the frozen hash table after checking every file byte for byte."""
    abm_root = Path(abm_root)
    for relative, expected in FROZEN_SHA256.items():
        actual = sha256(abm_root / relative)
        if actual != expected:
            raise RuntimeError('Frozen input changed: ' + relative)
    return dict(FROZEN_SHA256)


def adapter_sha256():
    return sha256(Path(__file__))


def verify_surface_pairs():
    """The hard-coded 22 pairs equal the archived extended manifest exactly."""
    if sha256(EXTENDED_MANIFEST) != EXTENDED_MANIFEST_SHA256:
        raise RuntimeError('Archived extended manifest changed')
    archived = json.loads(EXTENDED_MANIFEST.read_text())['resolved_config']['science']['pairs']
    if [dict(sorted(p.items())) for p in archived] != [dict(sorted(p.items())) for p in SURFACE_PAIRS]:
        raise RuntimeError('Surface pairs differ from the archived extended run')
    return archived


def _load(name, path, aliases=None):
    """Execute a frozen file under a private module name.

    ``aliases`` temporarily maps the bare names the frozen file imports
    (``coefficient``, ``background``) to the already loaded frozen modules, then
    restores whatever was registered before.
    """
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    aliases = aliases or {}
    previous = {key: sys.modules.get(key) for key in aliases}
    sys.modules.update(aliases)
    try:
        spec.loader.exec_module(module)
    finally:
        for key, value in previous.items():
            if value is None:
                sys.modules.pop(key, None)
            else:
                sys.modules[key] = value
    return module


verify_frozen()
bridge_coefficient = _load('_wfc_bridge_coefficient', BRIDGE / 'coefficient.py')
bridge_background = _load('_wfc_bridge_background', BRIDGE / 'background.py',
                          dict(coefficient=bridge_coefficient))
bridge_run = _load('_wfc_bridge_run', BRIDGE / 'run.py',
                   dict(coefficient=bridge_coefficient, background=bridge_background))
workflow_variant = _load('_wfc_workflow_variant', ROBUSTNESS / 'workflow_variant.py')
for _module, _relative in ((bridge_coefficient, 'research/full_icu_confirmation_2026-10-08/bridge/coefficient.py'),
                           (bridge_background, 'research/full_icu_confirmation_2026-10-08/bridge/background.py'),
                           (bridge_run, 'research/full_icu_confirmation_2026-10-08/bridge/run.py'),
                           (workflow_variant, 'research/robustness_uq_2026-10-17/workflow_variant.py')):
    if Path(_module.__file__).resolve() != (ABM / _relative).resolve():
        raise RuntimeError('Frozen module identity mismatch: ' + _relative)
if bridge_background.BridgeCoefficient is not bridge_coefficient.BridgeCoefficient:
    raise RuntimeError('Frozen background bound an unexpected coefficient module')

BridgeCoefficient = bridge_coefficient.BridgeCoefficient
FrozenVariant = workflow_variant.WorkflowRobustnessModel
VARIANT_METHODS = ('_clear_vacated_site', '_remove_discharged', '_move_patient',
                   '_relocate_detected_patients')
WORKFLOW_NAME = FrozenVariant.workflow_name


def rebind_frozen_method(function, owner):
    """Frozen function with its zero-argument ``super()`` anchored at ``owner``.

    Same code object, globals, defaults and docstring. Functions without a
    ``__class__`` free variable are returned unchanged.
    """
    freevars = function.__code__.co_freevars
    if not freevars:
        return function
    if freevars != ('__class__',):
        raise RuntimeError('Unexpected closure in frozen workflow method ' + function.__name__)
    clone = types.FunctionType(function.__code__, function.__globals__, function.__name__,
                               function.__defaults__, (types.CellType(owner),))
    clone.__kwdefaults__ = function.__kwdefaults__
    clone.__doc__ = function.__doc__
    clone.__qualname__ = owner.__qualname__ + '.' + function.__name__
    clone.__wrapped_frozen__ = function
    return clone


AUDIT_KEYS = ('moves', 'moves_from_isolation_bed', 'moves_into_isolation_bed',
              'moves_into_occupied_bed', 'moves_shared_to_isolation', 'swaps')
ENVIRONMENT_KEYS = ('discharge_clear', 'dirty_discharge_clear', 'relocation_clear', 'dirty_relocation_clear')
WORKFLOW_KEYS = ('isolation_relocation', 'isolation_deferred_no_empty_bed', 'occupied_destination_rejected')


def _validate_surface(surface_pairs, profiles, h0, h1):
    seen = set()
    extra, legacy = [], {}
    by_profile = {p['id']: p for p in profiles}
    for pair in surface_pairs:
        if pair['id'] in seen:
            raise ValueError('Duplicate surface pair id ' + pair['id'])
        seen.add(pair['id'])
        if pair.get('legacy_id'):
            profile = by_profile.get(pair['legacy_id'])
            if profile is None or (pair['h0'], pair['h1'], pair['rho'], pair['e']) != (h0, h1, profile['rho'], profile['e']):
                raise ValueError('Legacy surface pair does not match its profile arm: ' + pair['id'])
            legacy[pair['id']] = pair['legacy_id']
        else:
            extra.append(dict(pair, histogram_cap=pair.get('histogram_cap', SURFACE_HISTOGRAM_CAP)))
    return extra, legacy


def make_model_class(*, alternative, abm_root=ABM):
    """Return a BackgroundModel subclass for one workflow.

    Each call loads fresh private frozen classes (frozen ``load_core``), so an
    original-mode and an alternative-mode class can coexist in one process.
    Constructor: the frozen keywords (warmup_ticks, h0, h1, histogram_cap,
    profiles, model kwargs) plus ``surface_pairs`` (default: none).
    """
    if type(alternative) is not bool:
        raise TypeError('alternative must be a bool')
    verify_frozen(abm_root)
    Base = bridge_background.make_background_model(abm_root)

    class OriginalWorkflowBackground(Base):
        """Frozen BackgroundModel plus passive bookkeeping only."""
        workflow_alternative = False
        workflow_name = 'original_swap_isolation'

        def __init__(self, *, surface_pairs=(), **kwargs):
            # Passive tally of patient-contact events by period, carrier status
            # and bed type (isolation weight). Set before the frozen __init__.
            self.patient_contact_mix = Counter()
            super().__init__(**kwargs)
            if sum(self.bridge_counts.values()):
                raise RuntimeError('Bridge events before construction finished')
            first_policy = self.bridge_arms[1]
            extra, legacy = _validate_surface(surface_pairs, self.bridge_profiles,
                                              first_policy.h0, first_policy.h1)
            for i in range(len(self.bridge_profiles)):
                base_arm, policy_arm = self.bridge_arms[2*i], self.bridge_arms[2*i+1]
                if (base_arm.h0, base_arm.h1, policy_arm.h0, policy_arm.h1) != (first_policy.h0, first_policy.h0, first_policy.h0, first_policy.h1):
                    raise RuntimeError('Unexpected frozen profile arm layout')
            self.surface_pairs = [dict(p) for p in surface_pairs]
            self.surface_legacy = legacy
            self.surface_extra = extra
            self.surface_extra_offset = {p['id']: len(self.bridge_arms) + 2*i for i, p in enumerate(extra)}
            # Extra arms receive the identical event stream inside the frozen
            # bridge_event loop; they consume no RNG and do not touch the digest.
            self.bridge_arms.extend(BridgeCoefficient(p['h0'], h, p['rho'], p['e'], p['histogram_cap'])
                                    for p in extra for h in (p['h0'], p['h1']))
            self.bed_isolation_positions = frozenset(bed.pos for bed in self.isol_beds)
            self.bed_positions = frozenset(bed.pos for bed in self.shared_beds + self.isol_beds)
            self.relocation_audit = Counter(dict.fromkeys(AUDIT_KEYS, 0))
            self._audit_step_moves = []
            self._workflow_boundary = self.workflow_counters() if self.bridge_warmup_ticks == 0 else None

        # Frozen bridge_event (digest, counts, every arm) runs first and
        # unchanged; the tally below only reads its arguments.
        def bridge_event(self, kind, uid, *details):
            super().bridge_event(kind, uid, *details)
            if kind == 'patient':
                _, carrier, weight = details
                period = 'observation' if self.total_ticks > self.bridge_warmup_ticks else 'warmup'
                self.patient_contact_mix[(period, 'carrier' if carrier else 'susceptible',
                                          'isolation_bed' if weight != 1.0 else 'shared_bed')] += 1

        # Passive audit of every patient move made by any workflow rule.
        def _move_patient(self, patient, destination):
            origin = patient.pos
            occupants = tuple(a.unique_id for a in self.grid.get_cell_list_contents([destination])
                              if a.isPatient and a is not patient)
            super()._move_patient(patient, destination)
            audit = self.relocation_audit
            audit['moves'] += 1
            audit['moves_into_occupied_bed'] += bool(occupants)
            audit['moves_from_isolation_bed'] += origin in self.bed_isolation_positions
            audit['moves_into_isolation_bed'] += tuple(destination) in self.bed_isolation_positions
            audit['moves_shared_to_isolation'] += (origin not in self.bed_isolation_positions
                                                   and tuple(destination) in self.bed_isolation_positions)
            self._audit_step_moves.append((patient.unique_id, tuple(origin), tuple(destination)))

        def step(self):
            self._audit_step_moves = []
            super().step()
            moves = self._audit_step_moves
            if len(moves) > 1:
                # Occupant swap: within one tick, a moves X->Y and b moves Y->X.
                for i, (a, x, y) in enumerate(moves):
                    for b, x2, y2 in moves[i+1:]:
                        if a != b and x2 == y and y2 == x:
                            self.relocation_audit['swaps'] += 1
            if self.total_ticks == self.bridge_warmup_ticks:
                self._workflow_boundary = self.workflow_counters()

        def workflow_counters(self):
            if set(self.relocation_audit) != set(AUDIT_KEYS):
                raise RuntimeError('Unexpected relocation audit key')
            return {'audit.' + k: self.relocation_audit[k] for k in AUDIT_KEYS}

        def check_workflow_invariants(self):
            return None

        def bridge_result(self):
            result = super().bridge_result()
            self.check_workflow_invariants()
            pairs = {}
            for pair in self.surface_pairs:
                if pair['id'] in self.surface_legacy:
                    arms = result['profiles'][self.surface_legacy[pair['id']]]
                else:
                    offset = self.surface_extra_offset[pair['id']]
                    arms = {'baseline': self.bridge_arms[offset].result(),
                            'policy': self.bridge_arms[offset+1].result()}
                parameters = {k: pair[k] for k in ('h0', 'h1', 'rho', 'e')}
                if pair.get('legacy_id'):
                    parameters['legacy_id'] = pair['legacy_id']
                pairs[pair['id']] = dict(parameters=parameters, **arms)
            total = self.workflow_counters()
            boundary = self._workflow_boundary
            observed = None if boundary is None else {k: v - boundary.get(k, 0) for k, v in total.items()}
            mix ={'|'.join(key): count for key, count in sorted(self.patient_contact_mix.items())}
            if sum(mix.values()) != self.bridge_counts['patient']:
                raise RuntimeError('Patient-contact tally does not reconcile with bridge event counts')
            result['pairs'] = pairs
            result['workflow'] = dict(name=self.workflow_name, alternative=self.workflow_alternative,
                                      counters_total=total, counters_observation=observed,
                                      patient_contacts_by_period_status_bed=mix)
            return result

    if not alternative:
        return OriginalWorkflowBackground

    class AlternativeWorkflowBackground(OriginalWorkflowBackground):
        """Vacated-site clearing and empty-isolation-bed-only relocation."""
        workflow_alternative = True
        workflow_name = WORKFLOW_NAME

        def __init__(self, **kwargs):
            super().__init__(**kwargs)
            # State used by the frozen rule methods. Correspondence:
            #   instrumented_model.py:126  self.response_counts = {k: Counter() ...}; the rules
            #     read/write only its 'environment_events' Counter (plus workflow_events, below).
            #   instrumented_model.py:130  self.response_sites = [a ... if a.isGoo]
            #   instrumented_model.py:132  self.response_environment_count = sum(colonized sites)
            #   instrumented_model.py:136-138  site._response_environment_origin (None for a
            #     clean site; every site is clean at p=0 with a clean initial environment)
            #   workflow_variant.py:59-62  sites_by_zone = {care_zone_by_position[site.pos]: site}
            #   workflow_variant.py:63-66  response_site_by_bed_position = {bed.pos: sites_by_zone[bed.unique_id]}
            #   workflow_variant.py:67-68  exactly 30 fixed sites and 30 bed positions
            #   workflow_variant.py:69     response_counts['workflow_events'] = Counter()
            self.response_counts = {'environment_events': Counter(dict.fromkeys(ENVIRONMENT_KEYS, 0)),
                                    'workflow_events': Counter(dict.fromkeys(WORKFLOW_KEYS, 0))}
            self.response_sites = [a for a in self.schedule.agents if a.isGoo]
            self.response_environment_count = sum(a.colonized for a in self.response_sites)
            sites_by_zone = {self.care_zone_by_position[site.pos]: site for site in self.response_sites}
            self.response_site_by_bed_position = {bed.pos: sites_by_zone[bed.unique_id]
                                                  for bed in self.shared_beds + self.isol_beds}
            if len(sites_by_zone) != 30 or len(self.response_site_by_bed_position) != 30:
                raise RuntimeError('Workflow requires one fixed environmental site per bed')
            for site in self.response_sites:
                site._response_environment_origin = None
            if self.bridge_warmup_ticks == 0:
                self._workflow_boundary = self.workflow_counters()

        def workflow_counters(self):
            counters = super().workflow_counters()
            if hasattr(self, 'response_counts'):
                for group, keys in (('environment_events', ENVIRONMENT_KEYS), ('workflow_events', WORKFLOW_KEYS)):
                    if set(self.response_counts[group]) != set(keys):
                        raise RuntimeError('Unexpected workflow counter in ' + group)
                    counters.update({group + '.' + k: self.response_counts[group][k] for k in keys})
            return counters

        def check_workflow_invariants(self):
            """Frozen workflow invariants (workflow_variant.py:125-132) plus p=0 and no-swap checks."""
            environment = self.response_counts['environment_events']
            workflow = self.response_counts['workflow_events']
            audit = self.relocation_audit
            failures = []
            # workflow_variant.py:129
            if environment['discharge_clear'] != self.discharged_patients:
                failures.append('discharge clears != discharged patients')
            # workflow_variant.py:130
            if environment['relocation_clear'] != workflow['isolation_relocation']:
                failures.append('relocation clears != relocations')
            # workflow_variant.py:131-132
            for reason in ('discharge', 'relocation'):
                if environment['dirty_' + reason + '_clear'] > environment[reason + '_clear']:
                    failures.append('dirty clears exceed clears: ' + reason)
            # p=0 clean reservoir: no site is ever dirty, so no clear removes mass.
            if environment['dirty_discharge_clear'] or environment['dirty_relocation_clear']:
                failures.append('dirty clear in clean p=0 background')
            if self.response_environment_count != sum(a.colonized for a in self.response_sites) or self.response_environment_count:
                failures.append('site contamination count inconsistent or nonzero')
            # Rule semantics: every move is a detected patient into an empty isolation bed.
            if workflow['occupied_destination_rejected']:
                failures.append('occupied destination attempted')
            if audit['moves'] != workflow['isolation_relocation']:
                failures.append('audited moves != relocations')
            if audit['moves_into_occupied_bed'] or audit['moves_from_isolation_bed'] or audit['swaps']:
                failures.append('occupant swap or displacement observed')
            if audit['moves_shared_to_isolation'] != audit['moves']:
                failures.append('move not from a shared bed into an isolation bed')
            if failures:
                raise RuntimeError('Alternative workflow invariant failed: ' + '; '.join(failures))

    # Reuse the frozen rule methods as the same code objects.
    for name in VARIANT_METHODS:
        setattr(AlternativeWorkflowBackground, name,
                rebind_frozen_method(FrozenVariant.__dict__[name], AlternativeWorkflowBackground))
    return AlternativeWorkflowBackground


def archived_scientific_config():
    """The archived original batch config (bridge_all.json) as resolved by frozen run.normalize."""
    raw = json.loads((BRIDGE / 'configs/bridge_all.json').read_text())
    seedbook = json.loads((BRIDGE / 'seedbook.json').read_text())
    science = {k: v for k, v in raw.items() if k != 'seeds'}
    if science != seedbook['scientific_config'] or raw['seeds'] != seedbook['all_seeds']:
        raise RuntimeError('Archived config and seedbook disagree')
    return raw, bridge_run.normalize(raw)


def seed_banks():
    seedbook = json.loads((BRIDGE / 'seedbook.json').read_text())
    banks = {}
    for bank, seeds in seedbook['seed_blocks'].items():
        for seed in seeds:
            banks[seed] = bank
    if len(banks) != 256:
        raise RuntimeError('Expected 256 unique seeds in banks A and B')
    return banks


def make_factory(*, alternative, surface_pairs=(), abm_root=ABM):
    """Callable with the frozen ``simulate_seed`` factory signature."""
    cls = make_model_class(alternative=alternative, abm_root=abm_root)
    surface_pairs = tuple(dict(p) for p in surface_pairs)

    def factory(**kwargs):
        return cls(surface_pairs=surface_pairs, **kwargs)
    factory.model_class = cls
    return factory


def simulate_seed(seed, config, factory, deadline=None):
    """Frozen ``run.simulate_seed`` (horizon loop, boundary statistics) unchanged."""
    return bridge_run.simulate_seed(seed, config, factory, deadline)


__all__ = ['ABM', 'BRIDGE', 'FROZEN_SHA256', 'SURFACE_PAIRS', 'PRIMARY_PAIR_IDS', 'BridgeCoefficient',
           'FrozenVariant', 'VARIANT_METHODS', 'archived_scientific_config', 'adapter_sha256',
           'bridge_background', 'bridge_run', 'make_factory', 'make_model_class',
           'rebind_frozen_method', 'seed_banks', 'simulate_seed', 'verify_frozen', 'verify_surface_pairs']

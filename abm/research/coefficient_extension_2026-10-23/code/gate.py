"""Per-seed gate (PROTOCOL.md "Gate (every seed)").

Against the archived extended record of the same seed (the run behind
supplement Table S5; payload digest verified before use):

1. ordered bridge-event hash equal;
2. the 3 primary profiles equal (record ``profiles[legacy]`` versus the
   archived legacy pair, whose arms are the frozen cap-2 profile arms);
3. all 22 archived pairs equal, arm by arm;
4. every fine-grid row that repeats an archived pair equal to the archived
   arms (the two cap-dependent histogram keys are excluded only for the three
   primary duplicates; compared keys are stored);
5. (internal) every crossing row equal to the fine-grid row it repeats, when the
   record carries the fine grid.

"Bitwise" is implemented as equality of canonical JSON (sorted keys, ``repr``
floats, no NaN): two values compare equal only if every float has the same
binary64 value and the same sign of zero, and every container has the same
structure. A failed gate does not stop the record from being written; the
record stores the full comparison and ``passed=False``, and the seed is then
excluded and reported by summarize.py.
"""
import sys

sys.dont_write_bytecode = True

import hashlib  # noqa: E402
import json  # noqa: E402
from pathlib import Path  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from adapter import wb  # noqa: E402
import pairs  # noqa: E402

ARMS = ('baseline', 'policy')
PRIMARY_PAIR_IDS = dict(wb.PRIMARY_PAIR_IDS)
EXPOSURE_FIELDS = ('event_counts', 'executed_ticks', 'model_kwargs', 'observation_patient_days',
                   'observation_susceptible_patient_days', 'observation_admissions',
                   'observation_imported_colonizations', 'total_admissions', 'total_imported_colonizations')


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False)


def sha_text(text):
    return hashlib.sha256(text.encode()).hexdigest()


def bitwise_equal(a, b):
    return canonical(a) == canonical(b)


def extended_record_path(seed):
    paths = sorted(wb.EXTENDED_RESULTS.glob(f'part_*/completed/{seed}.json'))
    if len(paths) != 1:
        raise RuntimeError(f'Expected one archived extended record for seed {seed}, found {len(paths)}')
    return paths[0]


def load_extended(seed):
    """Archived extended record of ``seed`` with its payload digest verified."""
    if wb.sha256(wb.EXTENDED_MANIFEST) != wb.EXTENDED_MANIFEST_SHA256:
        raise RuntimeError('Archived extended manifest changed')
    path = extended_record_path(seed)
    raw = path.read_bytes()
    row = json.loads(raw)
    payload = {k: v for k, v in row.items() if k != 'payload_sha256'}
    if sha_text(canonical(payload)) != row.get('payload_sha256'):
        raise RuntimeError('Archived extended record payload digest mismatch: ' + str(path))
    if row.get('seed') != seed or row.get('status') != 'complete':
        raise RuntimeError('Archived extended record has wrong seed or status: ' + str(path))
    return dict(path=str(path), file_sha256=hashlib.sha256(raw).hexdigest(), row=row)


def _arm_detail(new, old, keys=None):
    keys = sorted(set(new) | set(old)) if keys is None else list(keys)
    equal = all(k in new and k in old and bitwise_equal(new[k], old[k]) for k in keys)
    out = dict(equal=equal)
    if not equal:
        out.update(K_new=new.get('K'), K_archived=old.get('K'),
                   differing_keys=sorted(k for k in keys if not (k in new and k in old and bitwise_equal(new[k], old[k]))))
    return out


def compare_arms(new_arms, old_arms, keys=None):
    detail = {arm: _arm_detail(new_arms[arm], old_arms[arm], keys) for arm in ARMS}
    return dict(equal=all(d['equal'] for d in detail.values()), arms=detail)


def run_gate(row, extended):
    """Gate comparison of one simulated record ``row`` with the archived ``extended`` row."""
    out = {}
    new_hash, old_hash = row['ordered_bridge_events_sha256'], extended['ordered_bridge_events_sha256']
    out['event_hash'] = dict(equal=new_hash == old_hash, new=new_hash, archived=old_hash)
    out['primary_profiles'] = {legacy: compare_arms(row['profiles'][legacy], extended['pairs'][pid])
                               for legacy, pid in sorted(PRIMARY_PAIR_IDS.items())}
    out['archived_pairs'] = {pid: compare_arms(row['pairs'][pid], extended['pairs'][pid]) for pid in pairs.ARCHIVED_IDS}
    duplicates = {}
    for new_id, archived_id in pairs.ARCHIVED_DUPLICATES:
        if new_id not in row['pairs']:
            continue
        keys = pairs.compared_keys_for_archived_duplicate(archived_id, row['pairs'][new_id]['baseline'])
        result = compare_arms(row['pairs'][new_id], extended['pairs'][archived_id], keys)
        result.update(archived_pair=archived_id, compared_keys=keys,
                      excluded_keys=sorted(set(row['pairs'][new_id]['baseline']) - set(keys)))
        if result['excluded_keys']:
            # Excluded only because the new arm has no histogram (cap None) by construction.
            result['excluded_keys_are_none_in_new_arm'] = all(
                row['pairs'][new_id][arm][k] is None for arm in ARMS for k in result['excluded_keys'])
            result['equal'] = result['equal'] and result['excluded_keys_are_none_in_new_arm']
        duplicates[new_id] = result
    out['archived_duplicates'] = duplicates
    internal = {}
    for cross_id, fine_id in pairs.NEW_DUPLICATES:
        if cross_id in row['pairs'] and fine_id in row['pairs']:
            internal[cross_id] = dict(compare_arms(row['pairs'][cross_id], row['pairs'][fine_id]), fine_pair=fine_id)
    out['internal_duplicates'] = internal
    out['informational_exposure_fields_equal'] = {k: bitwise_equal(row.get(k), extended.get(k)) for k in EXPOSURE_FIELDS}
    counts = dict(primary_profiles_equal=sum(v['equal'] for v in out['primary_profiles'].values()),
                  archived_pairs_equal=sum(v['equal'] for v in out['archived_pairs'].values()),
                  archived_duplicates_checked=len(duplicates),
                  archived_duplicates_equal=sum(v['equal'] for v in duplicates.values()),
                  internal_duplicates_checked=len(internal),
                  internal_duplicates_equal=sum(v['equal'] for v in internal.values()))
    out['counts'] = counts
    out['passed'] = bool(out['event_hash']['equal'] and counts['primary_profiles_equal'] == 3
                         and counts['archived_pairs_equal'] == 22
                         and counts['archived_duplicates_equal'] == counts['archived_duplicates_checked']
                         and counts['internal_duplicates_equal'] == counts['internal_duplicates_checked'])
    out['definition'] = ('passed = event hash equal AND 3 primary profiles AND 22 archived pairs AND every archived '
                         'duplicate AND every internal duplicate bitwise equal (canonical JSON). '
                         'informational_exposure_fields_equal is recorded, not part of passed.')
    return out


def gate_summary(gate):
    """Compact view of a stored gate result."""
    return dict(passed=gate['passed'], event_hash_equal=gate['event_hash']['equal'], **gate['counts'])


__all__ = ['ARMS', 'EXPOSURE_FIELDS', 'PRIMARY_PAIR_IDS', 'bitwise_equal', 'canonical', 'compare_arms',
           'extended_record_path', 'gate_summary', 'load_extended', 'run_gate', 'sha_text']

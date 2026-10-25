"""Every coefficient pair of the extension study (PROTOCOL.md "New pairs").

Order (fixed): the 22 archived Table S5 pairs first (``wb.SURFACE_PAIRS``,
literally, including the three ``legacy_id`` primaries that are read from the
frozen profile arms), then the new pairs, all ids prefixed ``x_``:

1. crossing (O-2)      h .90->.95, rho in {0,.5} x e in {.91,.92,.93,.94,.95,.96,.98}   14
2. completion (O-10)   h .60->.70 and .70->.80 at (rho,e) in {(0,1),(.5,1),(0,.8),(.5,.8)};
                       h .80->.90 and .95->.99 at (.5,.8)                                10
3. transition (new)    per-boundary h0, h1 giving WHO transition compliance 90% then 95%
                       (exact roots below), rho in {0,.5} x e in {1,.8}                   4
4. fine grid (O-6)     h .90->.95, rho in {0,.1,...,.9,.95} x e in {.5,.6,.7,.75,.78,.8,
                       .82,.85,.9,.95,.98,1}                                           132

160 new pairs, 182 in all, 364 coefficient arms (6 frozen profile arms + 358
extra arms). Every new pair is an extra ``BridgeCoefficient`` arm with the
adapter's default histogram cap (None), the same convention as the 19
archived non-primary pairs; K does not depend on the cap.

Duplicates are computed, never skipped:
* 8 fine-grid rows repeat an archived pair (rho in {0,.5} x e in {1,.9,.8,.6});
  each must equal the archived extended record bitwise. Three of them repeat a
  primary profile, whose archived arms carry the cap-2 histogram fields; for
  those the two histogram keys are excluded from the comparison (they are None
  in the new arm by construction) and the compared keys are recorded.
* 4 crossing rows repeat a fine-grid row (rho in {0,.5} x e in {.95,.98}); they
  must equal each other bitwise (an internal check, only on seeds that carry the
  fine grid).

Only ``id, h0, h1, rho, e`` reach the adapter; groups and duplicate maps live
in ``PAIR_INFO``.

Kinds (PROTOCOL "Budget and stopping rule"):
  all            archived 22 + all 160 new pairs
  no_fine_grid   archived 22 + crossing + completion + transition (50 pairs)
"""
import sys

sys.dont_write_bytecode = True

from pathlib import Path  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from adapter import wb  # noqa: E402
import closed_forms  # noqa: E402

H_CROSS = (.9, .95)
CROSSING_RHO = (0., .5)
CROSSING_E = (.91, .92, .93, .94, .95, .96, .98)
COMPLETION_FOUR = (((.6, .7), ((0., 1.), (.5, 1.), (0., .8), (.5, .8))),
                   ((.7, .8), ((0., 1.), (.5, 1.), (0., .8), (.5, .8))),
                   ((.8, .9), ((.5, .8),)),
                   ((.95, .99), ((.5, .8),)))
TRANSITION_TARGETS = ('0.90', '0.95')   # transition compliance before and after
TRANSITION_RHO = (0., .5)
TRANSITION_E = (1., .8)
FINE_RHO = (0., .1, .2, .3, .4, .5, .6, .7, .8, .9, .95)
FINE_E = (.5, .6, .7, .75, .78, .8, .82, .85, .9, .95, .98, 1.)
# Protocol text (rounded to 6 digits) for the transition roots; checked below.
PROTOCOL_TRANSITION_H6 = {0.: (.683772, .776393), .5: (.829180, .908392)}
HISTOGRAM_SENSITIVE_KEYS = ('K_by_opportunities', 'raw_pairs_without_reset_by_opportunities')


def pair_id(h0, h1, rho, e):
    """Archived id format (e.g. h0.9_to_0.95_rho0.5_e0.8)."""
    return f'h{h0:g}_to_{h1:g}_rho{rho:g}_e{e:g}'


def _pair(prefix, h0, h1, rho, e):
    return dict(id=prefix + pair_id(h0, h1, rho, e), h0=float(h0), h1=float(h1), rho=float(rho), e=float(e))


ARCHIVED = tuple(dict(p) for p in wb.SURFACE_PAIRS)
ARCHIVED_IDS = tuple(p['id'] for p in ARCHIVED)
LEGACY = {p['id']: p['legacy_id'] for p in ARCHIVED if p.get('legacy_id')}

CROSSING = tuple(_pair('x_cross_', *H_CROSS, rho, e) for rho in CROSSING_RHO for e in CROSSING_E)
COMPLETION = tuple(_pair('x_range_', h0, h1, rho, e) for (h0, h1), cells in COMPLETION_FOUR for rho, e in cells)

TRANSITION_ROOTS = {}
_transition = []
for _rho in TRANSITION_RHO:
    _records = [closed_forms.transition_root_record(_rho, target) for target in TRANSITION_TARGETS]
    _h0, _h1 = _records[0]['h_double'], _records[1]['h_double']
    if (round(_h0, 6), round(_h1, 6)) != PROTOCOL_TRANSITION_H6[_rho]:
        raise RuntimeError('Transition roots disagree with the protocol text')
    TRANSITION_ROOTS[f'rho{_rho:g}'] = dict(h0=_records[0], h1=_records[1])
    for _e in TRANSITION_E:
        _transition.append(dict(id=f'x_trans_c0.9_to_0.95_rho{_rho:g}_e{_e:g}', h0=_h0, h1=_h1, rho=_rho, e=_e))
TRANSITION = tuple(_transition)

FINE_GRID = tuple(_pair('x_fine_', *H_CROSS, rho, e) for rho in FINE_RHO for e in FINE_E)

GROUPS = (('crossing', CROSSING), ('completion', COMPLETION), ('transition', TRANSITION), ('fine_grid', FINE_GRID))
NEW = tuple(p for _, group in GROUPS for p in group)
ALL = ARCHIVED + NEW
KINDS = {'all': ALL,
         'no_fine_grid': ARCHIVED + CROSSING + COMPLETION + TRANSITION}


def _key(p):
    return (p['h0'], p['h1'], p['rho'], p['e'])


def _build_info():
    info = {}
    archived_by_key = {_key(p): p['id'] for p in ARCHIVED}
    for p in ARCHIVED:
        info[p['id']] = dict(group='archived', legacy_id=p.get('legacy_id'))
    for group, members in GROUPS:
        for p in members:
            info[p['id']] = dict(group=group)
    by_key = {}
    for p in NEW:
        by_key.setdefault(_key(p), []).append(p['id'])
    for p in NEW:
        archived = archived_by_key.get(_key(p))
        if archived is not None:
            info[p['id']]['duplicates_archived'] = archived
            info[p['id']]['archived_is_primary_profile'] = archived in LEGACY
        others = [q for q in by_key[_key(p)] if q != p['id']]
        if others:
            info[p['id']]['duplicates_new'] = others
    return info


PAIR_INFO = _build_info()
# Ordered (new_id, archived_id) and (crossing_id, fine_id) duplicate pairs.
ARCHIVED_DUPLICATES = tuple((p['id'], PAIR_INFO[p['id']]['duplicates_archived'])
                            for p in NEW if 'duplicates_archived' in PAIR_INFO[p['id']])
NEW_DUPLICATES = tuple((p['id'], q) for p in CROSSING for q in PAIR_INFO[p['id']].get('duplicates_new', ()))


def compared_keys_for_archived_duplicate(archived_id, arm_keys):
    """Keys compared between a new duplicate arm and its archived counterpart."""
    if archived_id in LEGACY:
        return sorted(k for k in arm_keys if k not in HISTOGRAM_SENSITIVE_KEYS)
    return sorted(arm_keys)


def pair_table():
    """Science table for the manifest: every pair with its group and duplicate map."""
    return [dict(p, **PAIR_INFO[p['id']]) for p in ALL]


def _self_check():
    ids = [p['id'] for p in ALL]
    if len(ids) != len(set(ids)):
        raise RuntimeError('Duplicate pair ids')
    if ALL[:22] != ARCHIVED or len(ARCHIVED) != 22:
        raise RuntimeError('Archived pairs must come first, all 22')
    if not all(p['id'].startswith('x_') and set(p) == {'id', 'h0', 'h1', 'rho', 'e'} for p in NEW):
        raise RuntimeError('New pair ids must be prefixed x_ and carry only id, h0, h1, rho, e')
    sizes = {name: len(members) for name, members in GROUPS}
    if sizes != {'crossing': 14, 'completion': 10, 'transition': 4, 'fine_grid': 132}:
        raise RuntimeError('Unexpected group sizes: ' + str(sizes))
    if len(ARCHIVED_DUPLICATES) != 8 or len(NEW_DUPLICATES) != 4:
        raise RuntimeError('Unexpected duplicate structure')
    if not all(a.startswith('x_fine_') for a, _ in ARCHIVED_DUPLICATES):
        raise RuntimeError('Only fine-grid rows may repeat an archived pair')
    if len(KINDS['all']) != 182 or len(KINDS['no_fine_grid']) != 50:
        raise RuntimeError('Unexpected kind sizes')


_self_check()

__all__ = ['ALL', 'ARCHIVED', 'ARCHIVED_DUPLICATES', 'ARCHIVED_IDS', 'COMPLETION', 'CROSSING', 'FINE_E', 'FINE_GRID',
           'FINE_RHO', 'GROUPS', 'HISTOGRAM_SENSITIVE_KEYS', 'KINDS', 'LEGACY', 'NEW', 'NEW_DUPLICATES', 'PAIR_INFO',
           'TRANSITION', 'TRANSITION_ROOTS', 'compared_keys_for_archived_duplicate', 'pair_id', 'pair_table']

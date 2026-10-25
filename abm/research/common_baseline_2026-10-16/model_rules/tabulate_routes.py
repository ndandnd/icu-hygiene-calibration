#!/usr/bin/env python3
"""Tabulate passive route|opportunity labels for observed acquisitions.

Reads the original full_icu_confirmation_2026-10-08 holdout records (448 cluster
+ 64 local seeds per profile) and writes route_opportunity_table.csv and
route_opportunity_table.md next to this script. Read-only on all inputs.

Observation-window quantities: arms.<arm>.observation.response_counts is
already final_snapshot minus warmup_snapshot (full_icu_response_2026-10-08/
run.py:107-111). This script re-derives that difference from the two
cumulative snapshots for every record and asserts equality, then asserts that
the joint route|opportunity histogram sums to observation.acquired_colonizations.

The labels are passive proximate-pickup labels (the pickup that began the
current uncleared hand-contamination episode); they are not causal shares.
"""
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
ABM = HERE.parents[2]
ROOT = ABM / 'results/full_icu_confirmation_2026-10-08'
FOLDERS = ('holdout_cluster', 'holdout_local')
PROFILE = {
    'independent_target_candidate': ('Independent', 0.0, 1.0),
    'persistent_target_candidate': ('Clustered misses', 0.5, 1.0),
    'imperfect_target_candidate': ('Imperfect clearing', 0.0, 0.8),
}
ORDER = list(PROFILE)
ARMS = ('baseline', 'policy')
ENV_ROUTES = ('own_zone_environment_pickup', 'other_zone_environment_pickup',
              'own_zone_environment_pickup_after_readmission')
S13 = {  # Supplementary Table S13 (supplement.tex:444-446): site/total
    ('independent_target_candidate', 'baseline'): (1717, 2375),
    ('independent_target_candidate', 'policy'): (627, 770),
    ('persistent_target_candidate', 'baseline'): (894, 2255),
    ('persistent_target_candidate', 'policy'): (405, 995),
    ('imperfect_target_candidate', 'baseline'): (885, 2352),
    ('imperfect_target_candidate', 'policy'): (671, 1718),
}
DAYS = 56
SITES = 30


def load():
    records = []
    for folder in FOLDERS:
        manifest = json.loads((ROOT / folder / 'manifest.json').read_text())
        jobs = {job['job_id']: job for job in manifest['jobs']}
        files = sorted((ROOT / folder / 'completed').glob('job_*.json'))
        assert len(files) == len(jobs), (folder, len(files), len(jobs))
        for path in files:
            record = json.loads(path.read_text())
            assert record['status'] == 'complete', path
            assert record['identity'] == manifest['identity'], path
            job = jobs[record['job_id']]
            assert record['seed'] == job['seed'] and record['scenario'] == job['scenario'], path
            record['_folder'] = folder
            records.append(record)
    return records


def observation_counts(record, arm):
    """Return the observation-window response_counts after verifying them."""
    start = record['warmup_snapshot']['response_counts']
    end = record['arms'][arm]['final_snapshot']['response_counts']
    stored = record['arms'][arm]['observation']['response_counts']
    derived = {kind: {k: v - start.get(kind, {}).get(k, 0) for k, v in counts.items()}
               for kind, counts in end.items()}
    assert derived == stored, (record['job_id'], arm)
    return stored


def main():
    records = load()
    per_scenario = Counter(r['scenario']['id'] for r in records)
    seeds = defaultdict(set)
    for r in records:
        seeds[r['scenario']['id']].add(r['seed'])
    assert all(per_scenario[s] == 512 and len(seeds[s]) == 512 for s in ORDER), per_scenario

    joint = defaultdict(Counter)
    totals = Counter()
    clears = Counter()
    checks = Counter()
    for r in records:
        sid = r['scenario']['id']
        for arm in ARMS:
            obs = r['arms'][arm]['observation']
            counts = observation_counts(r, arm)
            acquired = obs['acquired_colonizations']
            j = counts['acquisitions_by_route_and_opportunities']
            assert sum(j.values()) == acquired, (r['job_id'], arm)
            assert sum(counts['acquisitions_by_route'].values()) == acquired
            assert sum(counts['acquisitions_by_opportunities'].values()) == acquired
            assert counts['acquisitions_by_route'].get('untracked_hand', 0) == 0
            ev = counts['environment_events']
            assert ev.get('scheduled_clear', 0) == SITES * DAYS, (r['job_id'], arm, ev)
            assert obs['environment_contaminated_sites_endpoint_start'] == 0
            assert obs['environment_contaminated_sites_endpoint_end'] == 0
            clears[(sid, arm, 'scheduled_clear')] += ev.get('scheduled_clear', 0)
            clears[(sid, arm, 'dirty_scheduled_clear')] += ev.get('dirty_scheduled_clear', 0)
            clears[(sid, arm, 'deposition')] += ev.get('deposition', 0)
            for key, value in j.items():
                joint[(sid, arm)][key] += value
            totals[(sid, arm)] += acquired
            checks['arm_records'] += 1

    rows = []
    summary = []
    for sid in ORDER:
        name, rho, e = PROFILE[sid]
        for arm in ARMS:
            total = totals[(sid, arm)]
            cells = joint[(sid, arm)]
            assert sum(cells.values()) == total
            route_tot = Counter()
            for key, value in cells.items():
                route_tot[key.split('|')[0]] += value
            env = sum(route_tot[k] for k in ENV_ROUTES)
            exp_site, exp_total = S13[(sid, arm)]
            assert (env, total) == (exp_site, exp_total), (sid, arm, env, total)

            def gap_key(k):
                g = k.split('|')[1]
                return (k.split('|')[0], int(g) if g.isdigit() else 10**9)
            for key in sorted(cells, key=gap_key):
                if cells[key] == 0:
                    continue  # key present only because it occurred during warmup
                route, gap = key.split('|')
                rows.append(dict(profile=name, scenario_id=sid, rho=rho, e=e,
                                 arm='control (h=.90)' if arm == 'baseline' else 'policy (h=.95)',
                                 route=route, opportunities=gap, count=cells[key],
                                 arm_total_acquisitions=total,
                                 share_of_arm_total=round(cells[key] / total, 6)))
            env0 = cells.get('own_zone_environment_pickup|0', 0) + \
                cells.get('other_zone_environment_pickup|0', 0) + \
                cells.get('own_zone_environment_pickup_after_readmission|0', 0)
            summary.append(dict(
                profile=name, rho=rho, e=e, arm=arm, total=total,
                patient=route_tot['patient_pickup'],
                own=route_tot['own_zone_environment_pickup'],
                other=route_tot['other_zone_environment_pickup'],
                readm=route_tot['own_zone_environment_pickup_after_readmission'],
                env=env, env0=env0,
                env0_share=env0 / total, env_share=env / total,
                env0_of_env=env0 / env if env else float('nan'),
                own0=cells.get('own_zone_environment_pickup|0', 0),
                patient0=cells.get('patient_pickup|0', 0),
                other0=cells.get('other_zone_environment_pickup|0', 0),
                readm0=cells.get('own_zone_environment_pickup_after_readmission|0', 0),
                clears=clears[(sid, arm, 'scheduled_clear')],
                dirty_clears=clears[(sid, arm, 'dirty_scheduled_clear')],
                depositions=clears[(sid, arm, 'deposition')]))

    with (HERE / 'route_opportunity_table.csv').open('w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    with (HERE / 'route_opportunity_summary.csv').open('w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(summary[0]))
        w.writeheader()
        w.writerows(summary)

    pct = lambda x: f'{100 * x:.1f}'
    lines = [
        '# Route|opportunity labels of observed acquisitions, original holdout',
        '',
        'Source: `results/full_icu_confirmation_2026-10-08/holdout_cluster/completed` (448 seeds per',
        'profile) and `holdout_local/completed` (64 seeds per profile); 512 paired seeds per profile,',
        '1,536 records, all `status=complete`, record identities equal to their folder manifests.',
        'Generated by `tabulate_routes.py` in this folder.',
        '',
        '**What is counted.** `arms.<arm>.observation.response_counts` is the 56-day observation-window',
        'histogram: the runner stores `final_snapshot - warmup_snapshot`',
        '(`research/full_icu_response_2026-10-08/run.py:107-111`, used at :161 and :168). The script',
        're-derived that difference from the two cumulative snapshots for every record and arm and',
        'found exact agreement; it then checked, for every record and arm, that the joint',
        '`acquisitions_by_route_and_opportunities` histogram, the `acquisitions_by_route` histogram and',
        'the `acquisitions_by_opportunities` histogram each sum to `observation.acquired_colonizations`,',
        'and that `untracked_hand` is 0. Pooled surface counts reproduce Supplementary Table S13 exactly.',
        'Every arm of every record also has `environment_events.scheduled_clear = 1680` (30 sites x 56',
        'days) and zero contaminated sites at both daily-boundary endpoints, consistent with one',
        'synchronized complete cleaning of all sites per day and no other site-clearing event.',
        '',
        '**What the labels mean.** The route is the kind of the successful pickup that began the',
        'acquiring worker\'s current uncleared hand-contamination episode (patient, own-zone site,',
        'other-zone site, or own-zone site contaminated before the current admission to that bed).',
        'The opportunity number is how many hygiene opportunities (zone entries/exits under the',
        '`zone_separate` clock used here) that worker passed between that pickup and the transfer,',
        'without the hands being cleared. The pickup may have occurred during warmup. These are',
        'passive, proximate labels of realized paths. They are not causal or attributable shares:',
        'removing a route would also change later carriers, screening/isolation and contacts.',
        '',
        '## Summary by profile and arm',
        '',
        '| Profile (rho, e) | Arm | Acquisitions | Patient pickup | Own-zone site | Other-zone site | Own-zone site, earlier admission | Any site (S13) | Site pickup with 0 opportunities (own-zone + earlier-admission) | % of all | % of site-labelled |',
        '|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|',
    ]
    for s in summary:
        arm = 'control (.90)' if s['arm'] == 'baseline' else 'policy (.95)'
        lines.append(f"| {s['profile']} ({s['rho']:g}, {s['e']:g}) | {arm} | {s['total']} | {s['patient']} | "
                     f"{s['own']} | {s['other']} | {s['readm']} | {s['env']} ({pct(s['env_share'])}%) | "
                     f"{s['env0']} ({s['own0']} + {s['readm0']}) | {pct(s['env0_share'])} | {pct(s['env0_of_env'])} |")
    lines += ['',
              'Notes on reading the cells:',
              '',
              '- Under the entry/exit (`zone_separate`) clock every zone change costs an exit and an entry',
              '  opportunity, so a 0-opportunity acquisition needs pickup and transfer inside the same',
              '  care zone during one visit. In these records every such acquisition follows a site pickup',
              '  (wash on entry, touch the bedside site on the way to the bed, touch the patient). No',
              '  opportunity lies between that pickup and transfer, so completion changes this cell only',
              '  indirectly: through how often sites are contaminated, and through whether the entry wash',
              '  cleared earlier hand contamination. A pickup is registered only on clean hands',
              '  (`instrumented_model.py:82-85`); if the entry wash fails on already dirty hands, the same',
              '  site-to-patient visit is labelled with the earlier origin and at least 2 opportunities.',
              '  All observed gaps are even (exit+entry pairs).',
              '- `own_zone_environment_pickup_after_readmission` means the site contamination was deposited',
              '  before the acquiring patient\'s admission to that bed (a previous occupant\'s stay or an',
              '  empty-bed interval) and the acquiring patient is the latest admission to that bed',
              '  (`instrumented_model.py:186-188`). With no cleaning at discharge or admission, this is the',
              '  carry-over path that cleaning between occupants would largely interrupt; almost all such',
              '  acquisitions have 0 opportunities.',
              '- A relocation swap is not an admission (`response_admission_generation` changes only in',
              '  `admit_patient`, `instrumented_model.py:140-146`), so a susceptible patient swapped into a',
              '  detected carrier\'s shared bed who acquires from that bed\'s site is labelled plain',
              '  `own_zone_environment_pickup`, not `..._after_readmission`. The histograms cannot separate',
              '  swap-inherited contamination from contamination deposited during the current stay.',
              '- Keys whose observation-window count is 0 (present only because they occurred in warmup)',
              '  are omitted from the cell table below.',
              '',
              'Checked cells and site-event counts (pooled over 512 seeds):',
              '']
    for s in summary:
        arm = 'control' if s['arm'] == 'baseline' else 'policy'
        lines.append(f"- {s['profile']}, {arm}: patient_pickup|0 = {s['patient0']}, "
                     f"other_zone_environment_pickup|0 = {s['other0']}, "
                     f"own_zone_environment_pickup_after_readmission|0 = {s['readm0']}; "
                     f"site clears = {s['clears']} (dirty at clearing: {s['dirty_clears']}); "
                     f"hand-to-site depositions = {s['depositions']}.")
    lines += ['', 'Depositions equal dirty scheduled clears in every profile and arm, and both daily-boundary',
              'endpoints have zero dirty sites: every site contamination in the window ends at a scheduled',
              'daily clear, with no other clearing path.']
    lines += ['', '## Full route|opportunity cells', '',
              '| Profile | Arm | Route | Opportunities | Count | Share of arm total (%) |',
              '|---|---|---|---:|---:|---:|']
    for r in rows:
        lines.append(f"| {r['profile']} | {r['arm']} | {r['route']} | {r['opportunities']} | "
                     f"{r['count']} | {100 * r['share_of_arm_total']:.2f} |")
    (HERE / 'route_opportunity_table.md').write_text('\n'.join(lines) + '\n')
    print(json.dumps(dict(checks, records=len(records)), indent=1))
    for s in summary:
        print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in s.items()})


if __name__ == '__main__':
    main()

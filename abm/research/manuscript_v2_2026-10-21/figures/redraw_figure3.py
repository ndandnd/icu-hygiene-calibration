"""Redraw Figure 3 only, verifying its existing numerical content is unchanged."""
from pathlib import Path
import json
import math
import sys

import make_figures as figures

HERE = Path(__file__).resolve().parent
REVISION = HERE.parent / 'figure3_revision_2026-10-24'


def equal_values(old, new, path='values'):
    if isinstance(old, dict):
        if old.keys() != new.keys():
            raise AssertionError('Different keys: ' + path)
        for key in old:
            equal_values(old[key], new[key], path + '.' + key)
    elif isinstance(old, list):
        assert len(old) == len(new), path
        for i, (a, b) in enumerate(zip(old, new)):
            equal_values(a, b, path + f'[{i}]')
    elif isinstance(old, (float, int)) and not isinstance(old, bool):
        assert math.isclose(old, new, rel_tol=1e-12, abs_tol=1e-12), path
    else:
        assert old == new, path


def main():
    previous = json.loads((REVISION / 'before/figures/data_sources.json').read_text())
    alt, reason = figures.load_alt_coefficients(figures.ALT_COEF_SUMMARY, False)
    assert alt is not None, reason
    outputs = figures.fig_surface(HERE, alt, save_workflow_difference=False)
    current = figures.SOURCES['fig_surface']
    old_values = previous['figures']['fig_surface']['values']
    equal_values(old_values, {k: v for k, v in current['values'].items()
                              if k != 'display_semantics'})
    before = json.loads((REVISION / 'before/MANIFEST.json').read_text())
    for name in ('figures/fig_surface_diff.pdf', 'figures/fig_surface_diff.png'):
        assert figures.sha256(HERE.parent / name) == before[name]['sha256'], name
    receipt = dict(
        status='PASS', scope='Figure 3 presentation only; no new scientific result',
        scientific_values_unchanged=True, tolerance=1e-12,
        original_data_sources_sha256=figures.sha256(REVISION / 'before/figures/data_sources.json'),
        generator_sha256=figures.sha256(HERE / 'make_figures.py'),
        entrypoint_sha256=figures.sha256(Path(__file__)),
        python=sys.version.split()[0], numpy=figures.np.__version__,
        matplotlib=figures.matplotlib.__version__, figure=current,
        outputs={str(Path(p).name): figures.sha256(Path(p)) for p in outputs},
        workflow_difference_figures_unchanged=True,
    )
    (HERE / 'fig_surface_provenance.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({k: receipt[k] for k in ('status', 'scientific_values_unchanged', 'outputs')}, indent=2))


if __name__ == '__main__':
    main()

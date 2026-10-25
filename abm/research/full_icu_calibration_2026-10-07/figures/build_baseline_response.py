#!/usr/bin/env python3
"""Render the completed baseline grid as shared native/vector primitives.

No model execution and no policy inputs. PNG, PDF, SVG and portable inline
LaTeX use the same data and coordinates. Requires matplotlib for rendering.
"""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import tempfile

os.environ.setdefault('MPLCONFIGDIR', str(Path(tempfile.gettempdir())/'abm-baseline-figure-mpl'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

HERE = Path(__file__).resolve().parent
ABM = HERE.parents[2]
DEFAULT_REPORT = ABM/'results/full_icu_calibration_2026-10-07/baseline_grid_report/report.json'
PROFILES = (
    ('reference', 'Zone / independent / full clearing', '2463A6', 'o'),
    ('persistent', 'Zone / persistent / full clearing', 'C66A13', 's'),
    ('imperfect', 'Zone / independent / 80% clearing', '8657A6', 'd'),
    ('inherited', 'Per contact / independent / full clearing', '607080', 'open'),
)


def pale(color, white_fraction=.68):
    channels = [int(color[i:i+2], 16) for i in (0, 2, 4)]
    return ''.join(f'{round(value*(1-white_fraction)+255*white_fraction):02X}' for value in channels)


def tex_escape(value):
    replacements = {'&': r'\&', '%': r'\%', '_': r'\_', '#': r'\#'}
    return ''.join(replacements.get(c, c) for c in value)


class Scene:
    def __init__(self, width=540, height=432):
        self.width, self.height = width, height
        self.figure, self.ax = plt.subplots(figsize=(width/72, height/72))
        self.figure.subplots_adjust(0, 0, 1, 1)
        self.ax.set(xlim=(0, width), ylim=(0, height))
        self.ax.axis('off')
        self.layer = 0
        self.tex = [r'\begingroup\setlength{\unitlength}{1pt}',
                    r'\begin{picture}(%g,%g)' % (width, height)]

    def next_layer(self):
        self.layer += 1
        return self.layer

    def line(self, x0, y0, x1, y1, color='333333', width=.65):
        self.ax.plot([x0, x1], [y0, y1], color='#'+color, lw=width,
                     solid_capstyle='butt', zorder=self.next_layer())
        if abs(y1-y0) < 1e-8:
            self.tex.append(r'{\color[HTML]{%s}\put(%.4f,%.4f){\rule{%.4fpt}{%.4fpt}}}' %
                            (color, min(x0, x1), y0-width/2, max(abs(x1-x0), .001), width))
        elif abs(x1-x0) < 1e-8:
            self.tex.append(r'{\color[HTML]{%s}\put(%.4f,%.4f){\rule{%.4fpt}{%.4fpt}}}' %
                            (color, x0-width/2, min(y0, y1), width, max(abs(y1-y0), .001)))
        else:
            self.tex.append(r'{\color[HTML]{%s}\linethickness{%gpt}\qbezier(%.4f,%.4f)(%.4f,%.4f)(%.4f,%.4f)}' %
                            (color, width, x0, y0, (x0+x1)/2, (y0+y1)/2, x1, y1))

    def dashed(self, x0, x1, y, color='7C838A', width=.8):
        for x in range(int(x0), int(x1), 9):
            self.line(x, y, min(x+5, x1), y, color, width)

    def text(self, x, y, value, size=8, align='left', color='222222', rotation=0, weight='normal'):
        self.ax.text(x, y, value, fontsize=size, ha=align, va='center', color='#'+color,
                     rotation=rotation, weight=weight, family='DejaVu Sans', zorder=self.next_layer())
        anchor = {'left':'l', 'center':'c', 'right':'r'}[align]
        content = r'{\color[HTML]{%s}\fontsize{%g}{%g}\selectfont %s%s}' % (
            color, size, size*1.15, r'\bfseries ' if weight=='bold' else '', tex_escape(value))
        if rotation:
            content = r'\rotatebox{%g}{%s}' % (rotation, content)
        self.tex.append(r'\put(%.4f,%.4f){\makebox(0,0)[%s]{%s}}' % (x, y, anchor, content))

    def rectangle(self, x, y, width, height, color='FFFFFF', border=None):
        self.ax.add_patch(Rectangle((x, y), width, height, facecolor='#'+color,
                                   edgecolor='#'+border if border else 'none', linewidth=.5,
                                   zorder=self.next_layer()))
        self.tex.append(r'{\color[HTML]{%s}\put(%.4f,%.4f){\rule{%.4fpt}{%.4fpt}}}' % (color, x, y, width, height))
        if border:
            for a, b, c, d in ((x,y,x+width,y), (x+width,y,x+width,y+height),
                               (x+width,y+height,x,y+height), (x,y+height,x,y)):
                self.line(a,b,c,d,border,.5)

    def point(self, x, y, color, size=3.7, marker='o'):
        mpl_marker = 'o' if marker == 'open' else marker
        self.ax.plot(x, y, marker=mpl_marker, ms=size, mec='#'+color,
                     mfc='white' if marker=='open' else '#'+color, mew=.6, zorder=self.next_layer())
        if marker == 's':
            self.tex.append(r'{\color[HTML]{%s}\put(%.4f,%.4f){\rule{%.4fpt}{%.4fpt}}}' %
                            (color, x-size/2, y-size/2, size, size))
        elif marker == 'd':
            half = size/2
            for a,b,c,d in ((x-half,y,x,y+half), (x,y+half,x+half,y),
                           (x+half,y,x,y-half), (x,y-half,x-half,y)):
                self.line(a,b,c,d,color,.9)
        else:
            self.tex.append(r'{\color[HTML]{%s}\put(%.4f,%.4f){\circle%s{%g}}}' %
                            (color, x, y, '' if marker=='open' else '*', size))

    def save(self, directory):
        self.tex += [r'\end{picture}', r'\endgroup']
        inline = '\n'.join(self.tex)+'\n'
        portable = '% Requires xcolor and graphicx; no external images or data.\n'+r'\resizebox{\linewidth}{!}{%'+'\n'+inline+'}\n'
        (directory/'baseline_response.tex').write_text(portable)
        standalone = (r'\documentclass{article}'+'\n'+
                      r'\usepackage{xcolor,graphicx}'+'\n'+
                      r'\usepackage[paperwidth=580pt,paperheight=480pt,margin=20pt]{geometry}'+'\n'+
                      r'\pagestyle{empty}\setlength{\parindent}{0pt}'+'\n'+
                      r'\begin{document}'+'\n'+portable+r'\end{document}'+'\n')
        (directory/'baseline_response_standalone.tex').write_text(standalone)
        for suffix in ('png', 'pdf', 'svg'):
            self.figure.savefig(directory/f'baseline_response.{suffix}', dpi=300, facecolor='white',
                                metadata={'Title':'Full ICU baseline response',
                                          'Creator':'Completed-grid figure builder'} if suffix in ('pdf','svg') else None)
        plt.close(self.figure)


def checked_report(path):
    report = json.loads(path.read_text())
    if report['status'] != 'complete_baseline_grid' or report['complete_replicates'] != 224:
        raise ValueError('Only the complete declared 224-replicate baseline report is permitted')
    if report['seeds_per_cell'] != 8 or report['scenario_count'] != 28 or report['policy_outcomes_read']:
        raise ValueError('Figure requires the declared baseline-only eight-seed grid')
    if report['reference_scenario']['p'] != .3:
        raise ValueError('Prospectively selected training reference must be p=.3')
    for cell in report['cells']:
        if cell['scenario']['h0'] != .9 or cell['scenario']['h1'] != .9:
            raise ValueError('Policy input refused')
        values = list(cell['counts_by_seed'].values())
        mean = sum(values)/len(values)
        if len(values) != 8 or not math.isclose(mean, cell['observation']['acquired_colonizations']['mean']):
            raise ValueError('Replicate values do not reproduce reported cell mean')
    return report


def draw_curves(scene, groups, xx, yy, *, inset=False):
    for profile, _, color, marker in PROFILES:
        points = groups[profile]
        if inset:
            points = [point for point in points if point['scenario']['p'] <= .1]
        # Pale replicate dots, deterministic horizontal spread; means remain at p.
        for cell in points:
            p = cell['scenario']['p']
            for index, (_, value) in enumerate(sorted(cell['counts_by_seed'].items())):
                jitter = (index-3.5)*(.0010 if inset else .0016)
                adjusted = max(0., min(.1 if inset else .8, p+jitter))
                scene.point(xx(adjusted), yy(value), pale(color), size=1.6 if inset else 2.0)
        for left, right in zip(points, points[1:]):
            scene.line(xx(left['scenario']['p']), yy(left['observation']['acquired_colonizations']['mean']),
                       xx(right['scenario']['p']), yy(right['observation']['acquired_colonizations']['mean']),
                       color, 1.1 if inset else 1.45)
        for cell in points:
            stats = cell['observation']['acquired_colonizations']
            x = xx(cell['scenario']['p'])
            lower, upper = yy(stats['mean']-stats['se']), yy(stats['mean']+stats['se'])
            scene.line(x, lower, x, upper, color, .8)
            for y in (lower, upper):
                scene.line(x-2, y, x+2, y, color, .8)
            scene.point(x, yy(stats['mean']), color, size=3.4 if inset else 4.4, marker=marker)


def build(report, output):
    output.mkdir(parents=True, exist_ok=True)
    groups = {profile: sorted([c for c in report['cells'] if c['profile']==profile], key=lambda c:c['scenario']['p'])
              for profile, *_ in PROFILES}
    if any(len(group)!=7 for group in groups.values()):
        raise ValueError('Every hygiene profile needs its complete seven-point grid')
    scene = Scene()
    scene.text(12,417,'Full ICU baseline response to transmission probability',12,weight='bold')
    scene.text(12,399,'30 beds; 19 staff roles; 7 warmup days followed by 14 observation days',9)
    scene.text(12,383,'Completion is 90% in every cell. Pale dots: 8 seeds. Large markers: mean, with +/- 1 Monte Carlo SE.',8.3)
    x0,x1,y0,y1 = 72,522,136,343
    xx = lambda p:x0+p/.82*(x1-x0)
    yy = lambda count:y0+count/100*(y1-y0)
    for count in (0,20,40,60,80,100):
        scene.line(x0,yy(count),x1,yy(count),'E4E8EC',.5)
        scene.text(x0-8,yy(count),str(count),8.5,align='right',color='4B5157')
    for p in (0,.1,.2,.3,.5,.8):
        scene.line(xx(p),y0-3,xx(p),y0,'4B5157',.7)
        scene.text(xx(p),y0-13,f'{p:g}',8.5,align='center',color='4B5157')
    scene.line(x0,y0,x1,y0,'4B5157',.75)
    scene.line(x0,y0,x0,y1,'4B5157',.75)
    scene.text(19,(y0+y1)/2,'ICU-acquired colonizations per 14 days',9,align='center',rotation=90)
    scene.text((x0+x1)/2,107,'Common transmission probability p',9,align='center')
    target = report['reference_acquired_target']['mean']
    scene.dashed(x0,x1,yy(target),'828990',.8)
    draw_curves(scene,groups,xx,yy)
    scene.point(xx(.3),yy(target),'1A2734',size=8.5,marker='open')
    scene.text(517,yy(target)+11,f'Provisional eight-seed training target: {target:g}',8.2,align='right',color='535B63')
    inherited_max = max(c['observation']['acquired_colonizations']['mean'] for c in groups['inherited'])

    # Low-p inset in empty space below the high-burden curves; no data obscured.
    scene.rectangle(345,146,173,69,'FFFFFF','CDD3D9')
    scene.text(353,207,'Low-p detail',8,color='535B63')
    a,b,c,d = 363,507,162,198
    ix = lambda p:a+p/.1*(b-a)
    iy = lambda count:c+count/12*(d-c)
    for count in (0,6,12):
        scene.line(a,iy(count),b,iy(count),'EEF0F3',.45)
        scene.text(a-5,iy(count),str(count),6.6,align='right',color='69727B')
    for p in (0,.05,.1):
        scene.text(ix(p),153,f'{p:g}',6.6,align='center',color='69727B')
    scene.line(a,c,b,c,'69727B',.5)
    scene.line(a,c,a,d,'69727B',.5)
    draw_curves(scene,groups,ix,iy,inset=True)

    for (profile,label,color,marker),(x,y) in zip(PROFILES,((14,84),(282,84),(14,63),(282,63))):
        scene.line(x,y,x+14,y,color,1.4)
        scene.point(x+7,y,color,size=4.2,marker=marker)
        scene.text(x+21,y,label,8.3)
    scene.text(14,42,f'Persistent correlation: 0.5. Per-contact maximum mean: {inherited_max:g}; no baseline fit resolved on the grid.',8.1)
    scene.text(14,22,'Training only: a later 32-seed bank refined the reference target before validation.',8.5,color='535B63')
    scene.save(output)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--report', type=Path, default=DEFAULT_REPORT)
    parser.add_argument('--output', type=Path, default=HERE)
    args = parser.parse_args(argv)
    report = checked_report(args.report)
    build(report,args.output)
    caption = (
        'Baseline-only response curves in the full synthetic ICU model: 30 beds and 19 staff roles, '
        'with environment, admissions/discharges, screening and isolation retained. Each cell uses '
        'eight seeds, seven warmup days and 14 observation days, with washing completion fixed at 90%. '
        'Pale dots are individual seed outcomes, horizontally spread for visibility; large markers '
        'are means and bars are +/- one Monte Carlo standard error. Connecting lines aid reading '
        'and do not assert monotonicity or a continuous fitted response. Zone profiles place washing '
        'at the entry/exit of a patient-and-associated-environment care zone. Persistent completion '
        'has correlation 0.5; imperfect clearing has conditional efficacy 0.8. The inherited per-contact '
        'profile washes before/after each patient or environmental contact. The ring and dashed line '
        'mark the prospectively selected p=0.3 reference and provisional eight-seed acquisition target '
        '40.25, later refined with a separate 32-seed baseline bank. The inherited grid does '
        'not resolve a match to that target within the tested p values; this is not proof of mathematical '
        'infeasibility. The inset shows the low-transmission cells. These are training results; no '
        'policy outcomes or fresh-seed baseline-equivalence conclusions are shown.\n')
    (args.output/'baseline_response_caption.txt').write_text(caption)
    outputs = [args.output/f'baseline_response.{suffix}' for suffix in ('png','pdf','svg','tex')]
    outputs += [args.output/'baseline_response_standalone.tex',args.output/'baseline_response_caption.txt']
    try:
        source_label = str(Path(__file__).relative_to(ABM))
        data_label = str(args.report.resolve().relative_to(ABM))
    except ValueError:
        source_label, data_label = Path(__file__).name, args.report.name
    receipt = {'source': source_label, 'data': data_label,
               'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
               'data_sha256': hashlib.sha256(args.report.read_bytes()).hexdigest(),
               'outputs': {path.name:hashlib.sha256(path.read_bytes()).hexdigest() for path in outputs},
               'replicates':224,'seeds_per_cell':8,'complete_baseline_only':True,
               'uncertainty':'mean +/- one Monte Carlo SE, descriptive training summaries',
               'policy_outcomes_used':False,'matplotlib_version':matplotlib.__version__,
               'inline_dependencies':['xcolor','graphicx'],
               'note':'Common primitives generate native inline LaTeX and standalone figures. No model evaluations.'}
    (args.output/'baseline_response_manifest.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'outputs':[path.name for path in outputs],'replicates':224,'policy_outcomes_used':False},sort_keys=True))


if __name__ == '__main__':
    main()

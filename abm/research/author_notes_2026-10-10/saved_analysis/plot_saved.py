#!/usr/bin/env python3
"""Plot frozen saved weekly diagnostics; verify plot source/input before writing."""
from pathlib import Path
import hashlib
import json
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def canonical(obj):return json.dumps(obj,sort_keys=True,separators=(',',':'))
def write(path,value):
    with path.open('x') as stream:stream.write(json.dumps(value,indent=2,sort_keys=True)+'\n')

def main():
    source=HERE/'plot_saved.py';input_path=HERE/'outputs/weekly.json';freeze_path=HERE/'plot_freeze.json'
    if sys.argv[1:] == ['prepare']:
        payload=dict(schema_version=1,source_sha256=sha(source),weekly_input_sha256=sha(input_path),scope='prospective_plotting_of_existing_derived_results')
        payload['identity_sha256']=hashlib.sha256(canonical(payload).encode()).hexdigest();write(freeze_path,payload);return
    frozen=json.loads(freeze_path.read_text());payload={k:v for k,v in frozen.items() if k!='identity_sha256'}
    if hashlib.sha256(canonical(payload).encode()).hexdigest()!=frozen['identity_sha256'] or sha(source)!=frozen['source_sha256'] or sha(input_path)!=frozen['weekly_input_sha256']:raise ValueError('Plot freeze mismatch')
    rows=json.loads(input_path.read_text())['rows'];target=HERE/'outputs/weekly_control.png'
    if target.exists():raise FileExistsError(target)
    fig,axes=plt.subplots(3,1,figsize=(9.5,8.8),sharex=True)
    fig.subplots_adjust(left=.10,right=.985,top=.92,bottom=.14,hspace=.32)
    names=['independent','persistent','imperfect']; colors=['#19647e','#7e4380','#368052']
    for ax,name,color in zip(axes,names,colors):
        cells=[r for r in rows if r['cohort']=='pooled' and r['profile']==name+'_target_candidate']
        midpoint=[(r['day_start']+r['day_end'])/2 for r in cells];values=[r['rate_per_1000_patient_days']['estimate'] for r in cells]
        intervals=[r['rate_per_1000_patient_days']['uncertainty']['interval95'] for r in cells]
        ax.axvspan(0,90,color='#f2f0e8');ax.axvline(90,color='#777777',linestyle='--',linewidth=1)
        ax.errorbar(midpoint,values,yerr=[[v-i[0] for v,i in zip(values,intervals)],[i[1]-v for v,i in zip(values,intervals)]],color=color,fmt='o-',capsize=2,linewidth=1.2,markersize=4)
        ax.axhline(3,color='#aaaaaa',linestyle=':',linewidth=1)
        ax.set(title=name.capitalize(),ylabel='Acquisitions / 1,000 patient-days',xlim=(0,146),ylim=(0,5));ax.grid(axis='y',alpha=.18)
        ax.spines[['top','right']].set_visible(False)
    axes[0].text(40,4.7,'Warmup (90 days)',ha='center',color='#666666');axes[0].text(118,4.7,'Observation (56 days)',ha='center',color='#666666')
    axes[-1].set_xlabel('Elapsed simulation days (week midpoints)')
    fig.suptitle('Weekly baseline acquisition rates across the saved 146-day trajectories',fontsize=13)
    fig.text(.5,.025,'Pooled 512 seeds per profile; descriptive whole-seed bootstrap 95% intervals (50,000 draws).\nWeek 13 [84,91) spans the phase boundary; final week [140,146) contains six days.\nFinite-horizon diagnostics do not establish equilibrium or stationarity.',ha='center',fontsize=9)
    fig.savefig(target,dpi=160,bbox_inches='tight');plt.close(fig)
    write(HERE/'outputs/plot_receipt.json',dict(schema_version=1,status='complete',plot_freeze_sha256=sha(freeze_path),plot_freeze_identity=frozen['identity_sha256'],source_sha256=sha(source),input_hashes={'outputs/weekly.json':sha(input_path)},output_hashes={'weekly_control.png':sha(target)},software=dict(matplotlib=matplotlib.__version__)))

if __name__=='__main__':main()

"""Deterministic prospective design; never reads new simulation outcomes."""
import hashlib
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
PROFILES={"independent":dict(p=.14327558866543832,rho=0.,e=1.),
          "clustered":dict(p=.07602158239372622,rho=.5,e=1.),
          "imperfect":dict(p=.07483995991409355,rho=0.,e=.8)}

def write_new(path,value):
    text=json.dumps(value,indent=2,sort_keys=True,allow_nan=False)+"\n"
    path.parent.mkdir(parents=True,exist_ok=True)
    if path.exists():
        if path.read_text()!=text: raise ValueError("Refusing to change existing design: "+str(path))
    else: path.write_text(text)

def seeds(bank,n):
    return [int.from_bytes(hashlib.sha256(f"icu_robustness_uq_2026-10-17/{bank}/{i}".encode()).digest()[:4],"big") for i in range(n)]

def main():
    book={"grid":seeds("grid",128),"policy":seeds("policy",512),"probe":seeds("probe",3)}
    all_seeds=sum(book.values(),[])
    assert len(set(all_seeds))==len(all_seeds)
    cfg=dict(warmup_days=90,observation_days=56,block_days=1,
             model_overrides=dict(prob_patient_sick=.03,initial_environment_contaminated=False),
             description="Combined vacancy cleaning and empty-bed-only isolation; separately frozen robustness study.",
             seeds=book["grid"],scenarios=[])
    for name,profile in PROFILES.items():
        for label,m in (("m0p5",.5),("m1",1),("m2",2),("m4",4)):
            cfg["scenarios"].append(dict(id=name+"_"+label,p=profile["p"]*m,
                rho=profile["rho"],e=profile["e"],h0=.9,h1=.9,
                opportunities="zone_separate",environment_multiplier=1.))
    design=dict(schema_version=1,study=HERE.name,profiles=PROFILES,
        grid_multipliers=[.5,1,2,4],grid_seed_count=128,grid_records=1536,
        policy_seed_count=512,maximum_policy_pairs=1536,target_rate=3.,
        target_band=[2.7,3.3],target_checks="descriptive; not equivalence",
        bootstrap_draws=20000,bootstrap_seed=20261017,
        primary_contrasts=["independent-clustered","independent-imperfect","clustered-imperfect"],
        contrast_interval_level=1-.05/3,
        selection="lowest exact target hit, else lowest-p adjacent positive increasing bracket; log rate against log p; no fallback or extrapolation",
        stopping="fixed workload only; preserve missing bracket, undefined draws and failed target checks",
        workflow="clean bedside site on vacancy by discharge/relocation; isolate only into empty isolation bed; no occupant swaps",
        policy_uncertainty="conditional on selected p; calibration bootstrap reported separately; interpolation-model error excluded")
    write_new(HERE/"design.json",design);write_new(HERE/"seedbook.json",book)
    write_new(HERE/"configs/grid.json",cfg)
    probe=dict(cfg,warmup_days=1,observation_days=1,seeds=book["probe"],
               scenarios=[dict(cfg["scenarios"][i],id=name,h1=.95) for i,name in ((1,"independent"),(5,"clustered"),(9,"imperfect"))])
    write_new(HERE/"configs/probe.json",probe)
    print(json.dumps({"grid_records":1536,"maximum_policy_pairs":1536,"probe_pairs":9}))

if __name__=="__main__":main()

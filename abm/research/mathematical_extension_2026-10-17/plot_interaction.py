"""Plot the stationary two-opportunity result, not full-ICU estimates."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
rho = np.linspace(0, 0.9, 401)
h0, h1 = .9, .95

def prevention(r, e):
    q0 = (1-h0*e)**2+r*e*e*h0*(1-h0)
    q1 = (1-h1*e)**2+r*e*e*h1*(1-h1)
    return 100*(1-q1/q0)

threshold = json.loads((ROOT/"verification.json").read_text())["clearing_threshold_e"]
plt.rcParams.update({"font.size":11,"axes.spines.top":False,"axes.spines.right":False})
fig,ax=plt.subplots(figsize=(8.4,5.3),layout="constrained")
for e,color,label in ((1.,"#1c5679","Perfect clearing: e = 1"),
                       (.8,"#a45024","Imperfect clearing: e = 0.8"),
                       (threshold,"#74797c","Direction-change threshold: e ≈ 0.924")):
    ax.plot(rho,prevention(rho,e),color=color,lw=2.5,
            ls="--" if e==threshold else "-",label=label)
for e in (1.,.8):
    ax.scatter([0,.5],[prevention(0,e),prevention(.5,e)],color="#18272e",s=26,zorder=5)
ax.set(xlim=(-.01,.91),ylim=(20,80),xlabel="Persistence of completion/misses (ρ)",
       ylabel="Predicted prevention (%)",
       title="Clustering can change hygiene response in opposite directions")
ax.grid(axis="y",alpha=.2)
ax.legend(loc="upper right",frameon=False,fontsize=10)
fig.suptitle("Exact stationary two-opportunity bridge calculation; completion 90% → 95%\n"
             "Conditional model result, not a full-ICU estimate or clinical observation",fontsize=10,y=1.07)
fig.savefig(ROOT/"interaction.png",dpi=180,bbox_inches="tight")
fig.savefig(ROOT/"interaction.svg",bbox_inches="tight")
print(ROOT/"interaction.png")

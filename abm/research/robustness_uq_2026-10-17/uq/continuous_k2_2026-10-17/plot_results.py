"""Presentation of saved continuous-k2 results; adds no measurements."""
import hashlib
import json
import os
from pathlib import Path
import tempfile

os.environ.setdefault("MPLCONFIGDIR",tempfile.mkdtemp(prefix="abm-uq-cont-mpl-"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


root=Path(__file__).resolve().parent
results=json.loads((root/"run_001/results.json").read_text())
cases=[case for case in results["cases"] if case["n_per_assay"]==8000]
fig,axis=plt.subplots(figsize=(9.8,5.7))
x=np.arange(4)
series=(
    (-.25,"q1_only","continuous_width","#A4B0BD","One-opportunity assay, continuous"),
    (0,"q1_and_q2","grid_width","#89B5C7","Both assays, finite grid"),
    (.25,"q1_and_q2","continuous_width","#176A8A","Both assays, continuous"),
)
for offset,information,field,color,label in series:
    widths=[100*case["information_sets"][information][field]["quantiles_10_50_90"][1] for case in cases]
    bars=axis.bar(x+offset,widths,width=.23,color=color,label=label)
    axis.bar_label(bars,labels=[f"{value:.2f}" for value in widths],fontsize=9,padding=3)
axis.set_xticks(x,["Independent misses","Clustered misses","Imperfect clearing","Combined"])
axis.set_ylabel("Median prediction interval width (percentage points)")
axis.set_ylim(0,29)
axis.spines[["top","right"]].set_visible(False)
axis.set_axisbelow(True);axis.grid(axis="y",color="#E4E8EC",linewidth=.7)
fig.suptitle("Information gain persists over continuous parameters",fontsize=15,y=.97)
fig.text(.5,.895,"Two-opportunity bridges · known completion 90% → 95% · 8,000 measurements per assay",ha="center",fontsize=10)
fig.legend(*axis.get_legend_handles_labels(),loc="lower center",bbox_to_anchor=(.5,.08),ncol=1,frameon=False,fontsize=9)
fig.text(.5,.027,"Post hoc sensitivity of synthetic controlled assays. Widths condition on nonempty sets; adding the second assay doubles measurement cost.",ha="center",fontsize=7.8,color="#454D56")
fig.subplots_adjust(left=.09,right=.985,top=.82,bottom=.285)
for extension in ("png","svg"):
    fig.savefig(root/f"continuous_uncertainty.{extension}",dpi=190)
files=["plot_results.py","continuous_uncertainty.png","continuous_uncertainty.svg","run_001/results.json"]
(root/"figure_receipt.json").write_text(json.dumps({"scope":"Presentation of saved results only",
    "files":{name:hashlib.sha256((root/name).read_bytes()).hexdigest() for name in files}},indent=2)+"\n")

"""Plot the saved UQ benchmark; does not regenerate measurements or estimates."""
import json
import os
from pathlib import Path
import tempfile

os.environ.setdefault("MPLCONFIGDIR", tempfile.mkdtemp(prefix="abm-uq-mpl-"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


root = Path(__file__).resolve().parent
results = json.loads((root/"run_001/results.json").read_text())
cases = [case for case in results["cases"] if case["n_per_assay"] == 8000]
labels = ["Independent\nmisses", "Clustered\nmisses", "Imperfect\nclearing", "Combined"]
fig, axes = plt.subplots(1, 2, figsize=(11, 5.3), sharey=True)
for axis, index, title in zip(axes, (2, 3), ("Two-opportunity bridge", "Six-opportunity bridge")):
    x = np.arange(4)
    for offset, information, color, label in (
        (-.18, "q1_only", "#A4B0BD", "One-opportunity assay"),
        (.18, "q1_and_q2", "#176A8A", "One + two-opportunity assays"),
    ):
        widths = [100*case["information_sets"][information]["targets"][index]["width"]["quantiles_10_50_90"][1]
                  for case in cases]
        bars = axis.bar(x+offset, widths, width=.34, color=color, label=label)
        axis.bar_label(bars, labels=[f"{value:.1f}" for value in widths], padding=3, fontsize=9)
    axis.set_xticks(x, labels)
    axis.set_title(title, fontsize=12)
    axis.set_ylim(0, 55)
    axis.spines[["top", "right"]].set_visible(False)
    axis.set_axisbelow(True)
    axis.grid(axis="y", color="#E5E8EC", linewidth=.7)
axes[0].set_ylabel("Median prediction interval width (percentage points)")
fig.suptitle("Measuring longer survival narrows forecast uncertainty", fontsize=15, x=.5, y=.97)
fig.text(.5, .895, "Known completion: 90% → 95% · 8,000 controlled measurements per assay", ha="center", fontsize=10)
handles, legend_labels = axes[0].get_legend_handles_labels()
fig.legend(handles, legend_labels, loc="lower center", bbox_to_anchor=(.5,.09), ncol=2, frameon=False)
fig.text(.5,.035,"Synthetic stationary Markov grid; 1,000 experiments per truth. Widths exclude empty sets; coverage counts them as failures.",
         ha="center", fontsize=8.4, color="#454D56")
fig.subplots_adjust(top=.81, bottom=.25, left=.08, right=.985, wspace=.12)
fig.savefig(root/"uncertainty_contraction.png", dpi=190)
fig.savefig(root/"uncertainty_contraction.svg")

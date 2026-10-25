#!/usr/bin/env python3
"""Publication figures for manuscript v2 (ICU hand-hygiene bridge model).

Run with the project venv:

    $PROJECT_ROOT/abm/.venv/bin/python make_figures.py

Writes PDF + PNG for four figures and ``data_sources.json`` (input paths,
SHA256 and every plotted value) into this folder. Every number is read from a
machine-readable archive; nothing is digitized from an image. The script only
reads its inputs and writes into its own output directory.

Figures (drawn at their final printed width, so every font size is the printed size;
all text is at least 8 pt; round-1 revision 22 Oct 2026)
  fig_bridge_schematic  165 mm  vector schematic (no data)
  fig_surface           165 mm  leading-order F*(rho, e), original workflow, with the isolated
                                two-opportunity and renewal (F_geo) comparators
  fig_surface_diff      124 mm  supplementary: alternative minus original F*, same backgrounds
  fig_uq                165 mm  k=2 continuous interval widths
  fig_uq_k6             165 mm  supplementary: k=6 grid-conditional interval widths
  fig_icu               165 mm  finite-p F and acquisitions prevented, both workflows,
                                each workflow's own p->0 limit (hollow diamonds)

Alternative-workflow p->0 coefficients
  Read from ``workflow_coefficient_2026-10-21/summary.json`` (prespecified study,
  256 shared backgrounds). It is used only when status == "complete",
  self_checks_pass is true and its SHA256 matches RESULTS_SHA256.txt. Otherwise
  the figures are still drawn, but carry a visible "PLACEHOLDER" note where the
  alternative-workflow limits belong, and the status line says so.
  ``--alt-summary`` and ``--allow-incomplete-alt`` exist only to test that path
  into a scratch ``--outdir``.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import sys
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import patches as mpatches  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

# --------------------------------------------------------------------------
# Paths (read-only inputs)
# --------------------------------------------------------------------------
HERE = Path(__file__).resolve().parent
PROJECT = Path(__file__).resolve().parents[4]
ABM = PROJECT / "abm"

SUPPLEMENT = ABM / "manuscript" / "supplement.tex"
COEF_CSV = (PROJECT / "research" / "author_notes_2026-10-10" / "analysis_2026-10-11"
            / "bridge_pooled" / "coefficient_pairs.csv")
MATH_VERIFICATION = ABM / "research" / "mathematical_extension_2026-10-17" / "verification.json"
UQ_GRID = ABM / "research" / "robustness_uq_2026-10-17" / "uq" / "run_001" / "results.json"
UQ_CONT = (ABM / "research" / "robustness_uq_2026-10-17" / "uq" / "continuous_k2_2026-10-17"
           / "run_001" / "results.json")
ORIG_HOLDOUT = (ABM / "research" / "full_icu_confirmation_2026-10-08" / "holdout_report"
                / "report.json")
ALT_ANALYSIS = (PROJECT / "research" / "robustness_uq_2026-10-17" / "collection_2026-10-21"
                / "execution01" / "outputs" / "analysis.json")
ALT_VERIFIED_RESULTS = (ABM / "research" / "robustness_uq_2026-10-17"
                        / "verification_2026-10-21" / "RESULTS.md")
ALT_COEF_SUMMARY = ABM / "research" / "workflow_coefficient_2026-10-21" / "summary.json"
ALT_COEF_SHA_FILE = ABM / "research" / "workflow_coefficient_2026-10-21" / "RESULTS_SHA256.txt"
ALT_COEF_RESULTS = ABM / "research" / "workflow_coefficient_2026-10-21" / "RESULTS.md"
MATH_PLOT = ABM / "research" / "mathematical_extension_2026-10-17" / "plot_interaction.py"
# SHA256 of analysis.json quoted (abbreviated) in the verified report
ALT_ANALYSIS_SHA_PREFIX, ALT_ANALYSIS_SHA_SUFFIX = "063c0d95f376", "579f"

# --------------------------------------------------------------------------
# Style
# --------------------------------------------------------------------------
MM = 1 / 25.4
# Final printed widths: main figures use the 165.1-mm manuscript linewidth;
# compact supplementary panels retain the 124-mm width.
W1, W2 = 124 * MM, 165 * MM
ROUND1 = HERE.parent / "reviews_round1" / "analysis" / "round1_numbers.json"

# Okabe-Ito colour-blind-safe hues, assigned to entities in a fixed order.
OI = dict(black="#000000", orange="#E69F00", sky="#56B4E9", green="#009E73",
          yellow="#F0E442", blue="#0072B2", vermillion="#D55E00", purple="#CC79A7")
INK = "#1A1A1A"        # primary text / axes
INK2 = "#555555"       # secondary text
GRID = "#E3E3E3"
WORKFLOW_COLOR = {"original": OI["blue"], "alternative": OI["vermillion"]}
INFO_COLOR = {"q1_only": "#555555", "q1_and_q2": OI["green"]}
ACCENT = OI["purple"]  # contaminated hand in the schematic
# Ordinal clearing effectiveness e: one hue (violet), light -> dark, darkest = most
# effective. Passes the dataviz ordinal checks (monotone OKLCH L, adjacent dL >= .06,
# hue spread 7 deg, light end 3.2:1 on white). Violet keeps blue/vermillion free for
# the two workflows. Every level also has its own marker and a direct label.
E_LEVELS = (1.0, 0.9, 0.8, 0.6)
E_COLOR = {1.0: "#36096E", 0.9: "#523596", 0.8: "#6E5DAE", 0.6: "#8F8AC4"}
E_MARKER = {1.0: "o", 0.9: "s", 0.8: "^", 0.6: "D"}

PROFILES = (  # (key, display name, rho, e, legacy id used by the coefficient study)
    ("independent", "Independent misses", 0.0, 1.0, "iid_perfect"),
    ("clustered", "Clustered misses", 0.5, 1.0, "persistent_perfect"),
    ("imperfect", "Imperfect clearing", 0.0, 0.8, "iid_imperfect"),
)

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 8.5,
    "axes.titlesize": 8.5,
    "axes.labelsize": 8.5,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 8,
    "axes.linewidth": 0.6,
    "axes.edgecolor": INK,
    "axes.labelcolor": INK,
    "text.color": INK,
    "xtick.color": INK,
    "ytick.color": INK,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.major.size": 2.5,
    "ytick.major.size": 2.5,
    "xtick.minor.size": 1.5,
    "ytick.minor.size": 1.5,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "lines.linewidth": 1.0,
    "legend.frameon": False,
    "legend.handlelength": 1.6,
    "mathtext.fontset": "custom",
    "mathtext.rm": "Arial",
    "mathtext.it": "Arial:italic",
    "mathtext.bf": "Arial:bold",
    "mathtext.sf": "Arial",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
    "savefig.dpi": 600,
    "figure.dpi": 150,
})


# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------
def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


SOURCES: dict = {}


def register(fig: str, path: Path, role: str) -> None:
    SOURCES.setdefault(fig, {"inputs": [], "values": {}})
    SOURCES[fig]["inputs"].append({"path": str(path), "sha256": sha256(path), "role": role})


def record(fig: str, key: str, value) -> None:
    SOURCES.setdefault(fig, {"inputs": [], "values": {}})
    SOURCES[fig]["values"][key] = value


def check(cond: bool, msg: str) -> None:
    if not cond:
        raise SystemExit("CHECK FAILED: " + msg)


def panel_label(ax, letter: str, x: float = -0.02, y: float = 1.02) -> None:
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=10, fontweight="bold",
            ha="right", va="bottom")


def save(fig, outdir: Path, name: str) -> list[str]:
    paths = []
    for ext in ("pdf", "png"):
        p = outdir / f"{name}.{ext}"
        fig.savefig(p, dpi=600 if ext == "png" else None,
                    metadata=None if ext == "png" else {"Creator": "make_figures.py", "Producer": None})
        paths.append(str(p))
    plt.close(fig)
    return paths


def q2_response(rho, e, h0=0.9, h1=0.95):
    """Isolated stationary two-opportunity prevention 1 - q2(h1)/q2(h0)."""
    q0 = (1 - h0 * e) ** 2 + rho * e * e * h0 * (1 - h0)
    q1 = (1 - h1 * e) ** 2 + rho * e * e * h1 * (1 - h1)
    return 1 - q1 / q0


def geo_response(rho, e, h0=0.9, h1=0.95):
    """Renewal comparator F_geo = 1 - S(h1)/S(h0), S(h) = sum_j q_{2j}(h) = pi D P D (I-(PD)^2)^{-1} 1:
    one susceptible contact every two opportunities, unbounded stationary contamination episode."""
    def S(h):
        P = np.array([[rho + (1 - rho) * (1 - h), (1 - rho) * h],
                      [(1 - rho) * (1 - h), rho + (1 - rho) * h]])
        D = np.diag([1.0, 1.0 - e])
        M = P @ D
        pi = np.array([1 - h, h])
        return float((pi @ D @ M @ np.linalg.inv(np.eye(2) - M @ M)).sum())
    return 1 - S(h1) / S(h0)


# --------------------------------------------------------------------------
# Optional alternative-workflow coefficient summary (hook)
# --------------------------------------------------------------------------
def load_alt_coefficients(path: Path, allow_incomplete: bool):
    """Return (summary or None, reason). Never raises for an absent file."""
    if not path.exists():
        return None, "absent"
    try:
        summary = json.loads(path.read_text())
    except (OSError, ValueError):
        return None, "unreadable"
    if allow_incomplete:
        summary["_source_path"] = str(path)
        return summary, "test mode (completeness not enforced)"
    if summary.get("status") != "complete" or not summary.get("self_checks_pass"):
        return None, f"incomplete (status={summary.get('status')!r})"
    if ALT_COEF_SHA_FILE.exists() and path.resolve() == ALT_COEF_SUMMARY.resolve():
        sealed = dict(reversed(line.split()) for line in ALT_COEF_SHA_FILE.read_text().splitlines()
                      if line.strip())
        if sealed.get("summary.json") != sha256(path):
            return None, "SHA256 does not match RESULTS_SHA256.txt"
    summary["_source_path"] = str(path)
    return summary, "complete; SHA256 matches RESULTS_SHA256.txt"


def alt_pair_id(h0, h1, rho, e) -> str:
    return f"h{h0:g}_to_{h1:g}_rho{rho:g}_e{e:g}"


def placeholder(ax, text="PLACEHOLDER: alternative-workflow limits pending"):
    ax.text(0.5, 0.5, text, transform=ax.transAxes, ha="center", va="center", fontsize=7,
            color=OI["vermillion"], fontweight="bold", zorder=10,
            bbox=dict(facecolor="white", edgecolor=OI["vermillion"], lw=0.6, pad=2))


# --------------------------------------------------------------------------
# Table S5 parsing (cross-check of the CSV against the published table)
# --------------------------------------------------------------------------
def parse_table_s5() -> dict:
    tex = SUPPLEMENT.read_text()
    start = tex.index(r"\label{tab:author-coefficients}")
    end = tex.index(r"\end{longtable}", start)
    rows = {}
    for line in tex[start:end].splitlines():
        cells = [c.strip().rstrip("\\").strip() for c in line.split("&")]
        if len(cells) != 8:
            continue
        try:
            h0, h1, rho, e = (float(c) for c in cells[:4])
        except ValueError:
            continue
        m = re.match(r"([\d.]+)\s*\[([\d.]+),\s*([\d.]+)\]", cells[6])
        rows[(h0, h1, rho, e)] = dict(F=float(m.group(1)), lo=float(m.group(2)),
                                      hi=float(m.group(3)), q2=float(cells[7]))
    return rows


def load_coefficients() -> dict:
    import csv
    out = {}
    with open(COEF_CSV, newline="") as f:
        for r in csv.DictReader(f):
            if r["bank"] != "pooled":
                continue
            key = (float(r["h0"]), float(r["h1"]), float(r["rho"]), float(r["e"]))
            out[key] = dict(F=float(r["limiting_response"]), lo=float(r["ci95_low"]),
                            hi=float(r["ci95_high"]), q2=float(r["stationary_q2_response"]),
                            n=int(r["n"]), pair=r["pair"])
    return out


# ==========================================================================
# Figure 1: bridge schematic
# ==========================================================================
def token(ax, x, y, kind, r=1.45, lw=0.7):
    """Opportunity outcome token: 'clear' filled, 'miss' open, 'fail' half-filled."""
    if kind == "miss":
        ax.add_patch(mpatches.Circle((x, y), r, facecolor="white", edgecolor=INK, lw=lw, zorder=4))
    elif kind == "clear":
        ax.add_patch(mpatches.Circle((x, y), r, facecolor=INK, edgecolor=INK, lw=lw, zorder=4))
    elif kind == "fail":
        ax.add_patch(mpatches.Circle((x, y), r, facecolor="white", edgecolor=INK, lw=lw, zorder=4))
        ax.add_patch(mpatches.Wedge((x, y), r, 90, 270, facecolor=INK, edgecolor="none", zorder=5))
        ax.add_patch(mpatches.Circle((x, y), r, facecolor="none", edgecolor=INK, lw=lw, zorder=6))
    else:
        raise ValueError(kind)


def bracket(ax, x0, x1, y, label, h=1.4, lw=0.6):
    ax.plot([x0, x0, x1, x1], [y - h, y, y, y - h], color=INK, lw=lw, solid_capstyle="butt")
    ax.text((x0 + x1) / 2, y + 0.9, label, ha="center", va="bottom", fontsize=9)


def arrow(ax, p0, p1, color=INK, lw=0.7, ms=6, style="-|>", shrinkA=0, shrinkB=0, z=3):
    ax.add_patch(mpatches.FancyArrowPatch(p0, p1, arrowstyle=style, mutation_scale=ms, color=color,
                                          lw=lw, shrinkA=shrinkA, shrinkB=shrinkB, zorder=z))


def fig_bridge_schematic(outdir: Path):
    # Data coordinates span 175 x 70 units; the figure is printed 165 mm wide, so one unit is
    # 0.943 mm. All font sizes are printed sizes (>= 8 pt).
    W_units, H_units = 175, 74
    W_mm = 165
    fig = plt.figure(figsize=(W_mm * MM, W_mm * H_units / W_units * MM))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, W_units)
    ax.set_ylim(0, H_units)
    ax.set_aspect("equal")
    ax.axis("off")

    # ---------------- a: one bridge -------------------------------------
    yb = 52.0
    xc, xs = 11.0, 101.0           # carrier and susceptible contacts
    xo = (38.0, 66.0)              # two hygiene opportunities
    rp = 5.0
    ax.text(1.0, 73.0, "a", fontsize=10, fontweight="bold", va="top")
    # contaminated hand carried by one worker
    ax.plot([xc + rp + 0.8, xs - rp - 0.8], [yb, yb], color=ACCENT, lw=3.2,
            solid_capstyle="butt", zorder=1, alpha=0.9)
    arrow(ax, (xs - rp - 4.5, yb), (xs - rp - 0.4, yb), color=ACCENT, lw=0.01, ms=9, z=2)
    # patients
    ax.add_patch(mpatches.Circle((xc, yb), rp, facecolor=INK, edgecolor=INK, lw=0.8, zorder=3))
    ax.add_patch(mpatches.Circle((xs, yb), rp, facecolor="white", edgecolor=INK, lw=0.8, zorder=3))
    ax.text(xc, yb - rp - 2.0, "Carrier\ncontact", ha="center", va="top", fontsize=8)
    ax.text(xs, yb - rp - 2.0, "Susceptible\ncontact", ha="center", va="top", fontsize=8)
    ax.text((xc + xo[0]) / 2 + 1.0, yb + 1.6, "pickup", ha="center", va="bottom", fontsize=8,
            color=INK2, style="italic")
    ax.text((xo[1] + xs) / 2 - 1.0, yb + 1.6, "transfer", ha="center", va="bottom", fontsize=8,
            color=INK2, style="italic")
    # opportunities: exit from the carrier's zone, entry to the next zone (same worker)
    s = 6.4
    for i, x in enumerate(xo, 1):
        ax.add_patch(mpatches.FancyBboxPatch((x - s / 2, yb - s / 2), s, s,
                                             boxstyle="round,pad=0,rounding_size=1.2",
                                             facecolor="white", edgecolor=INK, lw=0.8, zorder=3))
        ax.text(x, yb, str(i), ha="center", va="center", fontsize=8.5, zorder=4)
    ax.text(xo[0], yb - s / 2 - 2.0, "exit from\ncarrier zone", ha="center", va="top", fontsize=8)
    ax.text(xo[1], yb - s / 2 - 2.0, "entry to\nnext zone", ha="center", va="top", fontsize=8)
    # survival brackets
    bracket(ax, xc + rp + 0.8, xo[0] + s / 2 + 1.5, yb + 7.5, r"$q_1$")
    bracket(ax, xc + rp + 0.8, xo[1] + s / 2 + 1.5, yb + 14.5, r"$q_2$")

    # ---------------- b: one opportunity --------------------------------
    ax.text(112.0, 73.0, "b", fontsize=10, fontweight="bold", va="top")
    x0, y0 = 119.0, 52.0           # opportunity box centre
    xb, yu, yd = 131.0, 61.0, 43.0  # completion junction, upper/lower branch heights
    xl = 142.0                      # leaf position
    yc, yn = yu + 6.0, yu - 4.0     # cleared / not-cleared leaves
    ax.plot([x0 - s / 2 - 4.0, x0 - s / 2], [y0, y0], color=ACCENT, lw=3.2, solid_capstyle="butt",
            alpha=0.9, zorder=1)
    ax.add_patch(mpatches.FancyBboxPatch((x0 - s / 2, y0 - s / 2), s, s,
                                         boxstyle="round,pad=0,rounding_size=1.2",
                                         facecolor="white", edgecolor=INK, lw=0.8, zorder=3))
    pink = dict(color=ACCENT, lw=2.0, alpha=0.9, solid_capstyle="butt", zorder=1)
    ax.plot([x0 + s / 2, xb], [y0, yu], **pink)            # completed (clearing trial pending)
    ax.plot([x0 + s / 2, xb, xl], [y0, yd, yd], **pink)    # missed: survives
    ax.plot([xb, xl], [yu, yn], **pink)                    # completed, not cleared: survives
    ax.plot([xb, xl], [yu, yc], color=INK, lw=0.6, zorder=1)  # completed and cleared: ends
    token(ax, xl + 2.2, yc, "clear")
    token(ax, xl + 2.2, yn, "fail")
    token(ax, xl + 2.2, yd, "miss")
    ax.text(xl + 5.0, yc, "cleared: ends", va="center", fontsize=8)
    ax.text(xl + 5.0, yn, "not cleared", va="center", fontsize=8)
    ax.text(xl + 5.0, yd, "missed", va="center", fontsize=8)
    # branch probabilities sit inside the fork, clear of tokens and labels
    ax.text(x0 + s / 2 + 4.6, y0 + 3.4, r"$h$", ha="left", va="center", fontsize=8.5, color=INK2)
    ax.text(x0 + s / 2 + 4.0, y0 - 4.4, r"$1-h$", ha="left", va="center", fontsize=8.5, color=INK2)
    ax.text((xb + xl) / 2, (yu + yc) / 2 + 1.0, r"$e$", ha="right", va="bottom", fontsize=8.5, color=INK2)
    ax.text((xb + xl) / 2 + 1.0, (yu + yn) / 2 - 1.4, r"$1-e$", ha="center", va="top", fontsize=8.5,
            color=INK2)
    ax.text(x0, y0 - s / 2 - 2.0, "Opportunity", ha="center", va="top", fontsize=8)

    # ---------------- c: three survival laws ----------------------------
    ax.text(1.0, 32.0, "c", fontsize=10, fontweight="bold", va="top")
    laws = (
        ("Independent misses", ["clear", "clear", "miss", "clear", "clear", "clear", "clear",
                                 "clear", "miss", "clear"],
         r"$q_k=(1-h)^k$", r"$h=0.9$: $q_1=0.10$, $q_2=0.010$"),
        ("Clustered misses",
         ["clear", "clear", "clear", "miss", "miss", "miss", "clear", "clear", "clear", "clear"],
         r"$q_k=(1-h)\,[\rho+(1-\rho)(1-h)]^{k-1}$",
         r"$\rho=0.5$: $q_1=0.10$, $q_2=0.055$"),
        ("Imperfect clearing",
         ["clear", "fail", "clear", "clear", "miss", "clear", "fail", "clear", "clear", "clear"],
         r"$q_k=(1-he)^k$", r"$e=0.8$: $q_1=0.28$, $q_2=0.078$"),
    )
    cols = (6.0, 64.0, 122.0)
    for (title, seq, formula, values), cx in zip(laws, cols):
        ax.text(cx, 29.0, title, fontsize=8.5, fontweight="bold", va="center")
        for j, kind in enumerate(seq):
            token(ax, cx + 1.6 + j * 4.3, 22.0, kind, r=1.4)
        ax.annotate("", xy=(cx + 1.6 + 9 * 4.3 + 3.6, 22.0), xytext=(cx + 1.6 + 9 * 4.3 + 1.9, 22.0),
                    arrowprops=dict(arrowstyle="-|>", color=INK2, lw=0.5, mutation_scale=5))
        ax.text(cx, 14.5, formula, fontsize=8.5, va="center")
        ax.text(cx, 6.5, values, fontsize=8, va="center", color=INK2)
    return save(fig, outdir, "fig_bridge_schematic")


# ==========================================================================
# Figure 3: leading-order surface
# ==========================================================================
def fig_surface(outdir: Path, alt, *, save_workflow_difference=True):
    name = "fig_surface"
    register(name, COEF_CSV, "pooled 256-background p=0 coefficients, original workflow (Table S5 source)")
    register(name, SUPPLEMENT, "Table S5 text, parsed for cross-check only")
    register(name, MATH_VERIFICATION, "math extension: exact isolated two-opportunity values and "
                                      "clearing_threshold_e (sign change of the rho effect)")
    register(name, MATH_PLOT, "math extension: isolated two-opportunity prevention formula")
    coef = load_coefficients()
    s5 = parse_table_s5()
    check(len(s5) == 22, f"Table S5 parsed {len(s5)} rows, expected 22")
    for key, row in s5.items():
        c = coef[key]
        check(abs(round(100 * c["F"], 2) - row["F"]) < 0.0051, f"F* mismatch at {key}")
        check(abs(round(100 * c["lo"], 2) - row["lo"]) < 0.0051, f"F* lower mismatch at {key}")
        check(abs(round(100 * c["hi"], 2) - row["hi"]) < 0.0051, f"F* upper mismatch at {key}")
        check(abs(round(100 * c["q2"], 2) - row["q2"]) < 0.0051, f"q2 mismatch at {key}")
        check(abs(c["q2"] - q2_response(key[2], key[3], key[0], key[1])) < 1e-12,
              f"q2 formula mismatch at {key}")
        check(c["n"] == 256, "pooled n != 256")
    mv = json.loads(MATH_VERIFICATION.read_text())
    check(mv.get("status") == "all_checks_passed", "math extension verification status")
    threshold = mv["clearing_threshold_e"]
    # isolated curve formula against the math extension's exact rational values
    n_exact = 0
    for blk in mv["interaction"]:
        e = blk["efficacy"]["decimal"]
        for r, F in zip(blk["rho"], blk["F"]):
            r = r["decimal"] if isinstance(r, dict) else float(r)
            check(abs(q2_response(r, e) - F["decimal"]) < 1e-12, f"isolated curve vs exact at rho={r}, e={e}")
            n_exact += 1
    check(abs(q2_response(0.0, threshold) - q2_response(0.75, threshold)) < 1e-9, "threshold flatness")
    check("(1-h0*e)**2+r*e*e*h0*(1-h0)" in MATH_PLOT.read_text().replace(" ", ""),
          "math extension formula text changed")

    rhos = (0.0, 0.25, 0.5, 0.75)
    # renewal comparator: cross-check against the round-1 analysis output
    r1 = json.loads(ROUND1.read_text())
    register(name, ROUND1, "round-1 analysis: renewal comparator F_geo for the 22 surface rows "
                           "(renewal_comparator.rows), cross-checked against the closed form here")
    n_geo = 0
    for key in s5:
        row = r1["renewal_comparator"]["rows"][alt_pair_id(*key)]
        check(abs(100 * geo_response(key[2], key[3], key[0], key[1]) - row["F_geo_pct"]) < 1e-9,
              f"F_geo mismatch at {key}")
        n_geo += 1

    # Separate the orthogonal symbol and line encodings from the data region.
    # All dimensions/fonts below refer to the final full-width manuscript figure.
    fig = plt.figure(figsize=(W2, 80 * MM))
    ax = fig.add_axes([0.085, 0.16, 0.585, 0.80])
    key_ax = fig.add_axes([0.725, 0.16, 0.27, 0.80])
    key_ax.set_xlim(0, 1)
    key_ax.set_ylim(0, 1)
    key_ax.set_axis_off()
    rr = np.linspace(0, 0.75, 301)
    plotted = {}
    halfwidth_max = 0.0
    for e in E_LEVELS:
        col = E_COLOR[e]
        ax.plot(rr, 100 * q2_response(rr, e), color=col, lw=0.7, ls=(0, (3.2, 1.6)), zorder=2)
        ax.plot(rr, [100 * geo_response(r, e) for r in rr], color=col, lw=0.9, ls=(0, (0.8, 1.4)),
                zorder=2)
        ys = [100 * coef[(0.9, 0.95, r, e)]["F"] for r in rhos]
        for r in rhos:
            c = coef[(0.9, 0.95, r, e)]
            halfwidth_max = max(halfwidth_max, 100 * (c["hi"] - c["lo"]) / 2)
        ax.plot(rhos, ys, color=col, lw=1.3, zorder=3)
        ax.plot(rhos, ys, ls="none", marker=E_MARKER[e], ms=5.0, mfc=col, mec="white", mew=0.5,
                zorder=4)
        plotted[f"e={e}"] = {f"rho={r}": dict(F_star_pct=100 * coef[(0.9, 0.95, r, e)]["F"],
                                              ci95_pct=[100 * coef[(0.9, 0.95, r, e)]["lo"],
                                                        100 * coef[(0.9, 0.95, r, e)]["hi"]],
                                              q2_isolated_pct=100 * coef[(0.9, 0.95, r, e)]["q2"],
                                              F_geo_pct=100 * geo_response(r, e))
                             for r in rhos}
    # isolated sign-change threshold: no rho effect at e = threshold
    ax.plot(rr, 100 * q2_response(rr, threshold), color="#9A9A9A", lw=0.7, ls=(0, (4, 1.5, 1, 1.5)),
            zorder=1)

    # alternative-workflow limits: checked here, reported in the supplementary difference figure
    alt_plotted, diffs = None, None
    if alt is not None and alt.get("surface"):
        register(name, Path(alt["_source_path"]), "alternative-workflow p=0 coefficients, same 256 "
                 "backgrounds (checked; plotted in fig_surface_diff)")
        surf = alt["surface"]
        check(surf.get("n_seeds") == 256, "alternative surface n_seeds != 256")
        check(alt.get("planned_surface") == "full", "alternative surface not full")
        check(len(surf["pairs"]) == 22, "alternative surface pair count != 22")
        alt_plotted, diffs, maxdiff_all = {}, {}, 0.0
        for key in s5:
            blk = surf["pairs"][alt_pair_id(*key)]
            check(blk["n"] == 256 and blk["alternative"]["n"] == 256, f"alt n at {key}")
            check(abs(blk["original"]["limiting_response"] - coef[key]["F"]) < 1e-9, f"orig mismatch {key}")
            maxdiff_all = max(maxdiff_all, abs(100 * (blk["alternative"]["limiting_response"]
                                                      - blk["original"]["limiting_response"])))
        check(maxdiff_all < 0.3, f"some |alt-orig| >= 0.3 pp ({maxdiff_all:.3f})")
        maxdiff_panel = 0.0
        for e in E_LEVELS:
            for r in rhos:
                pid = alt_pair_id(0.9, 0.95, r, e)
                blk = surf["pairs"][pid]
                a = blk["alternative"]
                y = 100 * a["limiting_response"]
                d = blk["delta_F_alt_minus_orig"]
                dest, div = 100 * d["estimate"], [100 * v for v in d["bootstrap_95"]["interval"]]
                check(abs(dest - (y - 100 * blk["original"]["limiting_response"])) < 1e-9, f"delta {pid}")
                check(d["bootstrap_95"]["undefined_draws"] == 0, f"undefined draws {pid}")
                maxdiff_panel = max(maxdiff_panel, abs(dest))
                iv = a["F_bootstrap_95"]["interval"]
                alt_plotted[pid] = dict(F_star_pct=y, ci95_pct=[100 * v for v in iv],
                                        original_F_star_pct=100 * blk["original"]["limiting_response"])
                diffs[pid] = dict(e=e, rho=r, delta_pp=dest, paired_bootstrap95_pp=div)
        check(maxdiff_panel < 0.06, f"panel |alt-orig| >= 0.06 pp ({maxdiff_panel:.3f})")
        record(name, "alternative_workflow_check", dict(
            path=alt["_source_path"], n_seeds=surf.get("n_seeds"), banks=surf.get("banks"),
            values=alt_plotted, max_abs_difference_pp_plotted_rows=maxdiff_panel,
            max_abs_difference_pp_all_22_rows=maxdiff_all))
    else:
        placeholder(ax, "PLACEHOLDER: alternative-workflow\ncoefficients pending")

    # Letters identify the three reference profiles without surrounding their
    # points with shapes that could be mistaken for uncertainty regions.
    prof_pts = {}
    label_positions = {"independent": (0.13, 80.0), "clustered": (0.50, 63.0),
                       "imperfect": (0.13, 35.0)}
    profile_letters = {"independent": "A", "clustered": "B", "imperfect": "C"}
    for key, label, rho, e, legacy in PROFILES:
        y = 100 * coef[(0.9, 0.95, rho, e)]["F"]
        prof_pts[key] = dict(rho=rho, e=e, F_star_pct_original=y)
        ax.annotate(profile_letters[key], xy=(rho, y), xytext=label_positions[key],
                    fontsize=9, fontweight="bold", color=INK, ha="center", va="center",
                    bbox=dict(facecolor="white", edgecolor="none", pad=0.8),
                    arrowprops=dict(arrowstyle="-", color=INK2, lw=0.6,
                                    shrinkA=3, shrinkB=4, connectionstyle="arc3,rad=0"),
                    zorder=6)

    ax.set_xlim(-0.05, 0.80)
    ax.set_ylim(8, 86)
    ax.set_xticks(rhos)
    ax.set_xticklabels(["0", "0.25", "0.5", "0.75"])
    ax.set_xlabel(r"Clustering of missed actions, $\rho$ (0 = independent)", fontsize=8.2)
    ax.set_ylabel(r"Relative prevented fraction $F^*$ (%)", fontsize=8.2)
    ax.yaxis.grid(True, color=GRID, lw=0.5)
    ax.set_axisbelow(True)
    for e, y in zip(E_LEVELS, (0.94, 0.86, 0.78, 0.70)):
        key_ax.plot([0.035], [y], marker=E_MARKER[e], ms=5,
                    color=E_COLOR[e], mec="white", mew=0.5, ls="none")
        key_ax.text(0.14, y, f"$e={e:g}$", fontsize=8.5,
                    ha="left", va="center")
    line_keys = [
        (0.55, "-", INK2, 1.3, "ICU calculation"),
        (0.45, (0, (3.2, 1.6)), INK2, 0.7, "Two-opportunity model"),
        (0.35, (0, (0.8, 1.4)), INK2, 0.9, "Renewal approximation"),
        (0.20, (0, (4, 1.5, 1, 1.5)), "#9A9A9A", 0.7,
         f"Two-opportunity threshold\n$e^*={threshold:.3f}$"),
    ]
    for y, style, colour, width, text in line_keys:
        key_ax.plot([0, 0.12], [y, y], ls=style, color=colour, lw=width)
        key_ax.text(0.17, y, text, fontsize=8, linespacing=1.15,
                    ha="left", va="center")
    record(name, "display_semantics", dict(
        marker_encodes="Probability e that a completed wash clears contamination",
        markers={str(e): E_MARKER[e] for e in E_LEVELS},
        profile_labels=profile_letters, profile_ovals=False,
        line_encodes="ICU calculation or approximation, not efficacy",
        context="Intervention, workflow, limit and A/B/C definitions stated in caption",
        printed_width_mm=165))
    record(name, "surface_h0.90_to_0.95_original", plotted)
    record(name, "profiles", prof_pts)
    record(name, "max_pointwise_ci95_halfwidth_pct_original", halfwidth_max)
    record(name, "isolated_sign_change_threshold_e", threshold)
    record(name, "isolated_curve_exact_checks", n_exact)
    record(name, "renewal_comparator_rows_checked", n_geo)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    for text in [*fig.texts, *ax.texts, *key_ax.texts, ax.xaxis.label, ax.yaxis.label]:
        box = text.get_window_extent(renderer)
        check(box.x0 >= 0 and box.y0 >= 0 and box.x1 <= fig.bbox.x1
              and box.y1 <= fig.bbox.y1, "Figure 3 label outside canvas: " + text.get_text())
    out = save(fig, outdir, name)

    # ---- supplementary: alternative minus original, same backgrounds ----
    if diffs and save_workflow_difference:
        nd = "fig_surface_diff"
        register(nd, Path(alt["_source_path"]), "surface.pairs.<pair>.delta_F_alt_minus_orig "
                 "(paired whole-seed bootstrap, 20,000 draws)")
        figd, axd = plt.subplots(1, 1, figsize=(W1, 62 * MM), layout="constrained")
        eoff = {1.0: -0.033, 0.9: -0.011, 0.8: 0.011, 0.6: 0.033}
        for pid, dd in diffs.items():
            e, r = dd["e"], dd["rho"]
            dest, div = dd["delta_pp"], dd["paired_bootstrap95_pp"]
            axd.errorbar([r + eoff[e]], [dest], yerr=[[dest - div[0]], [div[1] - dest]],
                         fmt=E_MARKER[e], ms=4.5, color=E_COLOR[e], mfc=E_COLOR[e], mec="white",
                         mew=0.4, elinewidth=0.8, capsize=0, zorder=3)
        axd.axhline(0, color=INK2, lw=0.5, zorder=1)
        axd.set_ylim(-0.04, 0.11)
        axd.set_yticks([0, 0.05, 0.1])
        axd.set_yticklabels(["0", "0.05", "0.10"])
        axd.set_xticks(rhos)
        axd.set_xticklabels(["0", "0.25", "0.5", "0.75"])
        axd.set_xlim(-0.06, 0.81)
        axd.set_xlabel(r"Lag-one correlation of misses, $\rho$ (markers dodged by $e$)")
        axd.set_ylabel(r"$\Delta F$*, alternative $-$ original (pp)")
        axd.yaxis.grid(True, color=GRID, lw=0.5)
        axd.set_axisbelow(True)
        hd = [Line2D([], [], ls="none", marker=E_MARKER[e], ms=4.5, mfc=E_COLOR[e], mec="white",
                     label=(r"$e=1$" if e == 1.0 else f"$e={e:g}$")) for e in E_LEVELS]
        axd.legend(handles=hd, loc="upper right", ncol=4, fontsize=8, handletextpad=0.2,
                   columnspacing=0.9, borderaxespad=0.1)
        record(nd, "paired_difference_alt_minus_orig_pp", diffs)
        out += save(figd, outdir, nd)
    return out


# ==========================================================================
# Figure 3: synthetic UQ widths
# ==========================================================================
def fig_uq(outdir: Path):
    name = "fig_uq"
    register(name, UQ_CONT, "continuous-domain k=2 widths (h 0.90->0.95)")
    register(name, UQ_GRID, "finite-grid k=6 widths at h=0.950 (grid-conditional)")
    cont = json.loads(UQ_CONT.read_text())
    grid = json.loads(UQ_GRID.read_text())
    truths = (("independent", "Independent misses", r"$\rho=0,\ e=1$"),
              ("clustered", "Clustered misses", r"$\rho=0.5,\ e=1$"),
              ("imperfect", "Imperfect clearing", r"$\rho=0,\ e=0.8$"),
              ("combined", "Clustered + imperfect", r"$\rho=0.5,\ e=0.8$"))
    ns = (500, 2000, 8000)
    sets = (("q1_only", r"$q_1$ only"), ("q1_and_q2", r"$q_1$ and $q_2$"))

    def cont_q(truth, n, s):
        c = next(c for c in cont["cases"] if c["truth"] == truth and c["n_per_assay"] == n)
        w = c["information_sets"][s]["continuous_width"]
        return [100 * v for v in w["quantiles_10_50_90"]], w["defined"], w["undefined"]

    def grid_q(truth, n, s, target="h=0.950,k=6"):
        c = next(c for c in grid["cases"] if c["truth"] == truth and c["n_per_assay"] == n)
        t = next(t for t in c["information_sets"][s]["targets"] if t["target"] == target)
        return [100 * v for v in t["width"]["quantiles_10_50_90"]], t["width"]["defined"], t["width"]["undefined"]

    # Cross-checks against the reported RESULTS tables (n = 8,000)
    exp_joint = dict(independent=6.27, clustered=1.13, imperfect=5.71, combined=6.69)
    exp_q1 = dict(independent=25.04, clustered=25.03, imperfect=11.45, combined=11.45)
    exp_k6 = dict(independent=(48.44, 2.39), clustered=(48.44, 3.68), imperfect=(16.21, 5.84),
                  combined=(16.21, 7.57))
    for t, *_ in truths:
        check(abs(round(cont_q(t, 8000, "q1_and_q2")[0][1], 2) - exp_joint[t]) < 0.006, f"k2 joint {t}")
        check(abs(round(cont_q(t, 8000, "q1_only")[0][1], 2) - exp_q1[t]) < 0.006, f"k2 q1 {t}")
        check(abs(round(grid_q(t, 8000, "q1_only")[0][1], 2) - exp_k6[t][0]) < 0.006, f"k6 q1 {t}")
        check(abs(round(grid_q(t, 8000, "q1_and_q2")[0][1], 2) - exp_k6[t][1]) < 0.006, f"k6 joint {t}")

    values = {"k2_continuous": {}, "k6_grid": {}}
    markers = {"q1_only": "o", "q1_and_q2": "s"}
    out = []
    for getter, rkey, nm, ylim, rowlab in ((cont_q, "k2_continuous", name, 30, None),
                                           (grid_q, "k6_grid", "fig_uq_k6", 55, None)):
        if nm != name:
            register(nm, UQ_GRID, "finite-grid k=6 widths at h=0.950 (grid-conditional)")
        fig, axes = plt.subplots(1, 4, figsize=(W2, 64 * MM), sharex=True, sharey=True,
                                 layout="constrained")
        fig.get_layout_engine().set(w_pad=0.02, wspace=0.03)
        for col, (t, title, params) in enumerate(truths):
            ax = axes[col]
            for s_, lab in sets:
                q = [getter(t, n, s_) for n in ns]
                lo = [v[0][0] for v in q]
                med = [v[0][1] for v in q]
                hi = [v[0][2] for v in q]
                colr = INFO_COLOR[s_]
                ax.fill_between(ns, lo, hi, color=colr, alpha=0.15, lw=0, zorder=1)
                ax.plot(ns, med, color=colr, lw=1.3, ls="-", zorder=3)
                ax.plot(ns, med, ls="none", marker=markers[s_], ms=4.5, mfc=colr, mec="white",
                        mew=0.5, zorder=4)
                values[rkey].setdefault(t, {})[s_] = {
                    str(n): dict(p10=v[0][0], median=v[0][1], p90=v[0][2], defined=v[1], undefined=v[2])
                    for n, v in zip(ns, q)}
            ax.set_xscale("log")
            ax.set_xticks(ns)
            ax.set_xticklabels(["500", "2,000", "8,000"])
            ax.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
            ax.set_xlim(380, 10500)
            ax.yaxis.grid(True, color=GRID, lw=0.5)
            ax.set_axisbelow(True)
            ax.set_title(f"{title}\n{params}", fontsize=8.5, pad=4)
        axes[0].set_ylabel("Interval width (pp)")
        axes[0].set_ylim(0, ylim)
        fig.supxlabel("Measurements per assay", fontsize=8.5)
        handles = [Line2D([], [], color=INFO_COLOR[s_], lw=1.3, ls="-", marker=markers[s_], ms=4.5,
                          mfc=INFO_COLOR[s_], mec="white", label=lab) for s_, lab in sets]
        handles.append(mpatches.Patch(color="#555555", alpha=0.15, lw=0, label="10th to 90th percentile"))
        fig.legend(handles=handles, loc="outside upper center", ncol=3, fontsize=8, columnspacing=1.6,
                   handlelength=2.4)
        out += save(fig, outdir, nm)
    record(name, "widths_pct_points", values["k2_continuous"])
    record("fig_uq_k6", "widths_pct_points", values["k6_grid"])
    record(name, "target", "prevention for completion 0.90 -> 0.95, isolated k=2 bridge, continuous (rho,e)")
    record("fig_uq_k6", "target", "prevention for completion 0.90 -> 0.95, isolated k=6 bridge, h=0.950 "
                                  "grid target (grid-conditional)")
    record(name, "summary_statistic", "median width among defined intervals; band = 10th to 90th percentile "
                                      "over 1,000 repeated synthetic experiments per cell; empty joint sets "
                                      "are undefined and excluded from widths")
    record(name, "scope", "synthetic controlled assays, known h, stationary first-order Markov completion "
                          "with constant clearing; q1+q2 uses twice the total measurements")
    return out


# ==========================================================================
# Figure 4: finite-p ICU results, both workflows
# ==========================================================================
def fig_icu(outdir: Path, alt):
    name = "fig_icu"
    register(name, ORIG_HOLDOUT, "original-workflow holdout report (512 pairs/profile)")
    register(name, ALT_ANALYSIS, "alternative-workflow cluster analysis (512 pairs/profile)")
    register(name, COEF_CSV, "original-workflow p->0 limits (pooled 256 backgrounds)")
    register(name, ALT_VERIFIED_RESULTS, "verified alternative-workflow report (values cross-checked)")
    alt_sha = sha256(ALT_ANALYSIS)
    check(alt_sha.startswith(ALT_ANALYSIS_SHA_PREFIX) and alt_sha.endswith(ALT_ANALYSIS_SHA_SUFFIX),
          "analysis.json SHA256 differs from the verified report")
    orig = json.loads(ORIG_HOLDOUT.read_text())
    altw = json.loads(ALT_ANALYSIS.read_text())
    coef = load_coefficients()
    orig_ids = {"independent": "independent_target_candidate", "clustered": "persistent_target_candidate",
                "imperfect": "imperfect_target_candidate"}

    rows = {}
    for key, label, rho, e, legacy in PROFILES:
        o = next(s for s in orig["scenarios"] if s["scenario"]["id"] == orig_ids[key])
        a = next(s for s in altw["scenarios"] if s["scenario"]["id"] == key)
        check(abs(o["scenario"]["rho"] - rho) < 1e-12 and abs(o["scenario"]["e"] - e) < 1e-12, f"orig {key}")
        check(abs(a["scenario"]["rho"] - rho) < 1e-12 and abs(a["scenario"]["e"] - e) < 1e-12, f"alt {key}")
        # both F are pooled-rate prevented fractions, 1 - policy rate / control rate
        check(abs(a["F"] - (1 - a["policy_rate"] / a["control_rate"])) < 1e-12, f"alt F definition {key}")
        lim = coef[(0.9, 0.95, rho, e)]
        rows[key] = dict(
            label=label,
            original=dict(p=o["scenario"]["p"], n=o["n"], F=100 * o["rate_prevented_fraction"],
                          F95=[100 * v for v in o["rate_prevented_fraction_uncertainty"]["interval95"]],
                          carriers=o["mean_prevented"], carriers95=list(o["mean_prevented_uncertainty"]["interval95"]),
                          control_rate=o["arms"]["baseline"]["rate_per_1000_patient_days"]),
            alternative=dict(p=a["scenario"]["p"], n=a["n"], F=100 * a["F"],
                             F95=[100 * v for v in a["F95"]["interval"]],
                             carriers=a["mean_carriers_prevented"],
                             carriers95=list(a["mean_carriers_prevented95"]["interval"]),
                             control_rate=a["control_rate"]),
            limit_original=dict(F=100 * lim["F"], ci95=[100 * lim["lo"], 100 * lim["hi"]],
                                q2_isolated=100 * lim["q2"]),
        )
    # cross-check against manuscript / verified report values
    exp = dict(independent=(67.6, 76.5, 3.13, 3.64, 75.13), clustered=(55.9, 53.7, 2.46, 2.31, 53.65),
               imperfect=(27.0, 27.3, 1.24, 1.17, 27.99))
    for k, (fo, fa, co, ca, fl) in exp.items():
        r = rows[k]
        check(round(r["original"]["F"], 1) == fo, f"orig F {k}")
        check(round(r["alternative"]["F"], 1) == fa, f"alt F {k}")
        check(round(r["original"]["carriers"], 2) == co, f"orig carriers {k}")
        check(round(r["alternative"]["carriers"], 2) == ca, f"alt carriers {k}")
        check(round(r["limit_original"]["F"], 2) == fl, f"limit {k}")
        check(r["original"]["n"] == 512 and r["alternative"]["n"] == 512, f"n {k}")

    alt_lim = None
    if alt is not None and alt.get("primary"):
        register(name, Path(alt["_source_path"]), "alternative-workflow p->0 limits "
                 "(primary.paired.<profile>.alternative.limiting_response, 256 shared backgrounds)")
        check(alt["primary"].get("n_seeds") == 256, "alternative primary n_seeds != 256")
        alt_lim = {}
        exp_alt = dict(independent=75.18, clustered=53.65, imperfect=27.99)
        posthoc = alt.get("post_hoc_finite_p", {}).get("by_profile", {})
        for key, label, rho, e, legacy in PROFILES:
            blk = alt["primary"]["paired"][legacy]
            a = blk["alternative"]
            check(blk["n"] == 256 and a["n"] == 256, f"alt primary n {key}")
            check(abs(blk["original"]["limiting_response"] - coef[(0.9, 0.95, rho, e)]["F"]) < 1e-9,
                  f"summary original limit vs CSV {key}")
            F = 100 * a["limiting_response"]
            check(round(F, 2) == exp_alt[key], f"alt limit {key}: {F:.4f}")
            check(abs(F - 100 * alt["surface"]["pairs"][alt_pair_id(0.9, 0.95, rho, e)]["alternative"]
                      ["limiting_response"]) < 1e-9, f"primary vs surface {key}")
            iv = a["F_bootstrap_95"]["interval"]
            alt_lim[key] = dict(F=F, ci95=[100 * v for v in iv],
                                finite_p_minus_limit_pp=rows[key]["alternative"]["F"] - F)
            ph = posthoc.get(legacy)
            if ph:   # consistency with the coefficient study's own post hoc block
                check(abs(100 * ph["finite_p_F"] - rows[key]["alternative"]["F"]) < 1e-9, f"post hoc F {key}")
                check(abs(ph["selected_p"] - rows[key]["alternative"]["p"]) < 1e-12, f"post hoc p {key}")
        record(name, "alternative_workflow_limits",
               dict(path=alt["_source_path"], n_seeds=alt["primary"].get("n_seeds"), values=alt_lim,
                    comparison="post hoc, descriptive (finite-p results were seen before the "
                               "coefficient protocol was written)"))

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(W2, 72 * MM), layout="constrained")
    fig.get_layout_engine().set(wspace=0.08)
    xs = np.arange(3)
    dx = 0.19
    lim_off = 0.085                    # limit marker sits just right of its own estimate
    mk = {"original": "o", "alternative": "s"}
    keys = [p[0] for p in PROFILES]
    for wf, sign in (("original", -1), ("alternative", 1)):
        col = WORKFLOW_COLOR[wf]
        for ax, val, ci in ((ax1, "F", "F95"), (ax2, "carriers", "carriers95")):
            y = np.array([rows[k][wf][val] for k in keys])
            lo = np.array([rows[k][wf][ci][0] for k in keys])
            hi = np.array([rows[k][wf][ci][1] for k in keys])
            ax.errorbar(xs + sign * dx, y, yerr=[y - lo, hi - y], fmt=mk[wf], ms=5.0, color=col,
                        mfc=col, mec="white", mew=0.5, elinewidth=1.0, capsize=2.2, capthick=0.9,
                        zorder=4)
    # each workflow's own p->0 limit: a hollow diamond just right of its own estimate
    for i, k in enumerate(keys):
        lims = {"original": rows[k]["limit_original"]["F"]}
        if alt_lim:
            lims["alternative"] = alt_lim[k]["F"]
        for wf, sign in (("original", -1), ("alternative", 1)):
            if wf not in lims:
                continue
            xc = i + sign * dx + lim_off
            ax1.plot([xc], [lims[wf]], ls="none", marker="D", ms=5.0, mfc="white",
                     mec=WORKFLOW_COLOR[wf], mew=1.1, zorder=5)
    if not alt_lim:
        placeholder(ax1)
    for ax in (ax1, ax2):
        ax.set_xticks(xs)
        ax.set_xticklabels([rows[k]["label"].replace(" ", "\n") for k in keys])
        ax.set_xlim(-0.55, 2.55)
        ax.yaxis.grid(True, color=GRID, lw=0.5)
        ax.set_axisbelow(True)
        ax.tick_params(axis="x", length=0)
    ax1.set_ylim(0, 90)
    ax1.set_ylabel("Prevented fraction (%)")
    ax2.set_ylim(0, 4.5)
    ax2.set_ylabel("Acquisitions prevented per run")
    panel_label(ax1, "a", x=-0.12)
    panel_label(ax2, "b", x=-0.12)
    def _h(wf, kind):
        col = WORKFLOW_COLOR[wf]
        if kind == "est":
            return Line2D([], [], color=col, marker=mk[wf], ms=5.0, mec="white", ls="none",
                          label=f"{wf.capitalize()}, finite $p$ [95% interval]")
        return Line2D([], [], ls="none", marker="D", ms=5.0, mfc="white", mec=col, mew=1.1,
                      label=f"{wf.capitalize()}, $p\\to0$ limit")
    handles = [_h("original", "est"), _h("original", "lim"), _h("alternative", "est")]
    if alt_lim:
        handles.append(_h("alternative", "lim"))
    fig.legend(handles=handles, loc="outside upper center", ncol=2, fontsize=8,
               handlelength=1.6, columnspacing=2.0)
    record(name, "values", rows)
    record(name, "definitions", dict(
        original_F="report.json scenarios[].rate_prevented_fraction (pooled-rate based) with "
                   "rate_prevented_fraction_uncertainty (95%, 200,000 whole-seed draws)",
        alternative_F="analysis.json scenarios[].F = 1 - policy_rate/control_rate with F95.interval "
                      "(95%, 20,000 synchronous whole-seed draws)",
        carriers="mean paired count reduction per 56-day observation window (mean_prevented / "
                 "mean_carriers_prevented) with 95% bootstrap intervals",
        limit_original="coefficient_pairs.csv bank=pooled limiting_response (256 original-workflow p=0 "
                       "backgrounds); equals summary.json primary.paired.<profile>.original",
        limit_alternative="workflow_coefficient summary.json primary.paired.<profile>.alternative."
                          "limiting_response (same 256 backgrounds, alternative workflow rules)",
        caveat="the two workflows' finite-p runs are separate studies with different selected p and "
               "fresh seeds; side by side is context only. Finite-p versus limit is descriptive and "
               "post hoc."))
    return save(fig, outdir, name)


# ==========================================================================
# Palette check (Machado 2009 CVD simulation, OKLab Delta E x100)
# ==========================================================================
def _lin(hexc):
    c = np.array([int(hexc[i:i + 2], 16) / 255 for i in (1, 3, 5)])
    return np.where(c <= 0.04045, c / 12.92, ((c + 0.055) / 1.055) ** 2.4)


def _oklab(lin):
    M1 = np.array([[0.4122214708, 0.5363325363, 0.0514459929],
                   [0.2119034982, 0.6806995451, 0.1073969566],
                   [0.0883024619, 0.2817188376, 0.6299787005]])
    M2 = np.array([[0.2104542553, 0.7936177850, -0.0040720468],
                   [1.9779984951, -2.4285922050, 0.4505937099],
                   [0.0259040371, 0.7827717662, -0.8086757660]])
    return M2 @ np.cbrt(M1 @ np.clip(lin, 0, 1))


CVD = {
    "protan": np.array([[0.152286, 1.052583, -0.204868], [0.114503, 0.786281, 0.099216],
                        [-0.003882, -0.048116, 1.051998]]),
    "deutan": np.array([[0.367322, 0.860646, -0.227968], [0.280085, 0.672501, 0.047413],
                        [-0.011820, 0.042940, 0.968881]]),
    "tritan": np.array([[1.255528, -0.076749, -0.178779], [-0.078411, 0.930809, 0.147602],
                        [0.004733, 0.691367, 0.303900]]),
}


def _contrast(hexc, bg="#FFFFFF"):
    def lum(h):
        return float(np.dot([0.2126, 0.7152, 0.0722], _lin(h)))
    a, b = sorted((lum(hexc), lum(bg)), reverse=True)
    return (a + 0.05) / (b + 0.05)


def _oklch(hexc):
    L, a, b = _oklab(_lin(hexc))
    return float(L), float(math.hypot(a, b)), float(math.degrees(math.atan2(b, a)) % 360)


def ordinal_report(cols_light_to_dark) -> dict:
    lch = [_oklch(c) for c in cols_light_to_dark]
    Ls = [x[0] for x in lch]
    hues = [x[2] for x in lch]
    gaps = [Ls[i] - Ls[i + 1] for i in range(len(Ls) - 1)]
    spread = max(hues) - min(hues)
    res = dict(colors_light_to_dark=cols_light_to_dark, oklch=[[round(v, 3) for v in x] for x in lch],
               lightness_monotone=all(g > 0 for g in gaps), min_adjacent_dL=round(min(gaps), 3),
               hue_spread_deg=round(spread, 1), light_end_contrast=round(_contrast(cols_light_to_dark[0]), 2))
    res["pass"] = (res["lightness_monotone"] and res["min_adjacent_dL"] >= 0.06 and spread <= 40
                   and res["light_end_contrast"] >= 2.0)
    return res


def palette_report() -> dict:
    sets = {"workflow": list(WORKFLOW_COLOR.values()), "information_set": list(INFO_COLOR.values())}
    out = {}
    for nm, cols in sets.items():
        pairs = []
        for i in range(len(cols)):
            for j in range(i + 1, len(cols)):
                a, b = _lin(cols[i]), _lin(cols[j])
                d = {"normal": 100 * float(np.linalg.norm(_oklab(a) - _oklab(b)))}
                for k, M in CVD.items():
                    d[k] = 100 * float(np.linalg.norm(_oklab(M @ a) - _oklab(M @ b)))
                pairs.append(dict(pair=[cols[i], cols[j]], deltaE_oklab_x100={k: round(v, 1) for k, v in d.items()},
                                  min_cvd=round(min(d[k] for k in CVD), 1)))
        adjacent = [p for p in pairs if abs(cols.index(p["pair"][0]) - cols.index(p["pair"][1])) == 1]
        out[nm] = dict(colors=cols, contrast_on_white={c: round(_contrast(c), 2) for c in cols}, pairs=pairs,
                       worst_adjacent_cvd_deltaE=min(p["min_cvd"] for p in adjacent),
                       worst_adjacent_normal_deltaE=round(min(p["deltaE_oklab_x100"]["normal"] for p in adjacent), 1),
                       worst_cvd_deltaE=min(p["min_cvd"] for p in pairs),
                       worst_normal_deltaE=round(min(p["deltaE_oklab_x100"]["normal"] for p in pairs), 1))
    out["clearing_e_ordinal"] = ordinal_report([E_COLOR[e] for e in sorted(E_LEVELS)])
    check(out["clearing_e_ordinal"]["pass"], "e ramp fails ordinal checks")
    check(out["workflow"]["worst_cvd_deltaE"] >= 8 and out["workflow"]["worst_normal_deltaE"] >= 15,
          "workflow colours fail CVD/normal-vision separation")
    out["note"] = ("Python re-implementation of the dataviz skill checks (Machado et al. 2009 full-severity "
                   "CVD matrices; OKLab Delta E x100; WCAG contrast); the skill's validate_palette.py gave the "
                   "same verdicts: workflow pair passes all categorical checks; e ramp passes the ordinal "
                   "checks (one hue, monotone lightness, adjacent dL >= .06, light end >= 2:1). The "
                   "information-set pair is deliberate emphasis (neutral grey reference for q1 only, green for "
                   "q1+q2), so the grey's chroma-floor result does not apply; the pair passes CVD and "
                   "normal-vision separation. Every series also has its own marker shape or line style.")
    return out


# ==========================================================================
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", type=Path, default=HERE)
    ap.add_argument("--alt-summary", type=Path, default=ALT_COEF_SUMMARY,
                    help="testing only: alternative-workflow coefficient summary.json")
    ap.add_argument("--allow-incomplete-alt", action="store_true",
                    help="testing only: accept a partial/test summary (refused for the canonical outdir)")
    args = ap.parse_args(argv)
    outdir = args.outdir.resolve()
    if args.allow_incomplete_alt and outdir == HERE:
        raise SystemExit("--allow-incomplete-alt may only write to a scratch --outdir")
    outdir.mkdir(parents=True, exist_ok=True)

    alt, reason = load_alt_coefficients(args.alt_summary, args.allow_incomplete_alt)
    status = f"USED ({reason})" if alt is not None else f"NOT USED ({reason}): PLACEHOLDER notes drawn"
    print(f"alternative-workflow coefficient summary ({args.alt_summary}): {status}")

    outputs = {}
    outputs["fig_bridge_schematic"] = fig_bridge_schematic(outdir)
    outputs["fig_surface"] = fig_surface(outdir, alt)
    outputs["fig_uq"] = fig_uq(outdir)
    outputs["fig_icu"] = fig_icu(outdir, alt)
    manuscript = ABM / "manuscript" / "manuscript.tex"
    register("fig_bridge_schematic", manuscript, "eq. (bridge-survival) and the two closed forms drawn in panel c")
    mtex = manuscript.read_text().replace(" ", "")
    check(r"$(1-he)^k$" in mtex and r"$(1-h)[\rho+(1-\rho)(1-h)]^{k-1}$" in mtex,
          "closed-form survival laws not found in manuscript")
    SOURCES["fig_bridge_schematic"]["values"]["note"] = (
        "Schematic only. Formulas from manuscript eq. (bridge-survival); q1/q2 values computed "
        "in closed form at h=0.90 (rho=0.5 for clustered, e=0.8 for imperfect).")
    SOURCES["fig_bridge_schematic"]["values"]["q_values_h0.90"] = {
        "independent": dict(q1=0.1, q2=0.1 ** 2),
        "clustered_rho0.5": dict(q1=0.1, q2=0.1 * (0.5 + 0.5 * 0.1)),
        "imperfect_e0.8": dict(q1=1 - 0.9 * 0.8, q2=(1 - 0.9 * 0.8) ** 2)}
    check(abs(q2_response(0.0, 1.0) - 0.75) < 1e-12, "q2 independent")

    manifest = dict(
        generated_by=str(Path(__file__).resolve()), script_sha256=sha256(Path(__file__).resolve()),
        python=sys.version.split()[0], numpy=np.__version__, matplotlib=matplotlib.__version__,
        alternative_coefficient_summary=dict(path=str(args.alt_summary), used=alt is not None, reason=reason,
                                             status=None if alt is None else alt.get("status"),
                                             sha256=sha256(args.alt_summary) if args.alt_summary.exists()
                                             else None),
        outputs=outputs, figures=SOURCES, palette=palette_report())
    if outdir == HERE:
        (outdir / "data_sources.json").write_text(json.dumps(manifest, indent=1, default=float) + "\n")
    else:
        (outdir / "data_sources_test.json").write_text(json.dumps(manifest, indent=1, default=float) + "\n")
    for k, v in outputs.items():
        print(k, *v)


if __name__ == "__main__":
    main()

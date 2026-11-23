"""Render Figs. 9 and 11 from the numerical summary tables."""

import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
import numpy as np

from figure_style import ALLOWED_COLOR, apply_style

ROOT = Path(__file__).resolve().parents[1]


def rows(relative):
    with (ROOT / relative).open(newline="") as stream:
        return list(csv.DictReader(stream))


def finish(figure, axes, output, name):
    for ax in axes:
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid()
    figure.savefig(output / f"{name}.pdf", bbox_inches="tight")
    plt.close(figure)


def area_robustness(output):
    target = next(row for row in rows("gaps/saving_methods.csv") if
                  row["shape"] == "triangle" and row["order"] == "3"
                  and row["filling_pct"] == "2.1937")
    variants = rows("gaps/paired_factors.csv")
    curves = rows("gaps/paired_curves.csv")
    factors = [float(r["quadratic_factor"]) for r in variants]
    intervals = [
        (float(target["method_low"]), float(target["method_high"])),
        (float(target["shifted_low"]), float(target["shifted_high"])),
        (min(factors), max(factors)),
    ]
    labels = ["Crossing method", r"Gap width $\pm0.2$ points", "Paired refinement"]
    fig, (left, right) = plt.subplots(1, 2, figsize=(9.5, 4.2), layout="constrained")
    for y, (lo, hi) in zip((2, 1, 0), intervals):
        left.hlines(y, lo, hi, colors=".25", linewidth=5)
        left.text(.98, y+.16, f"{lo:.2f}–{hi:.2f}", ha="right",
                  transform=left.get_yaxis_transform(),
                  bbox={"facecolor": "white", "edgecolor": "none", "pad": 1})
    left.axvline(float(target["saving"]), color=".4", linestyle=":", linewidth=.9)
    left.set_yticks([2,1,0], labels)
    left.set(xlim=(1.94, 2.20), ylim=(-.4, 2.57),
             xlabel=r"Saving factor $f_{\rm eq}/f$", title="(a) Sensitivity ranges")
    dense = [(float(r["filling_pct"]), float(r["widest_gap_pct"]))
             for r in rows("gaps/smooth_curves.csv") if r["shape"] == "triangle"
             and 3.4 <= float(r["filling_pct"]) <= 5.6]
    right.plot(*np.array(dense).T, color=".2", linestyle="--", label="Dense reference")
    palette = [ALLOWED_COLOR, "#b95824", "#238066", "#7656a4", "#a77515"]
    for row, color, marker, factor in zip(variants, palette, ("o","s","^","D","v"), factors):
        variant = row["variant"]
        curve = [(float(r["filling_pct"]), float(r["widest_gap_pct"])) for r in curves
                 if r["variant"] == variant and r["geometry"] == "smooth"]
        fk3 = next(float(r["widest_gap_pct"]) for r in curves
                   if r["variant"] == variant and r["geometry"] == "FK3")
        right.plot(*np.array(curve).T, marker=marker, markersize=4,
                   color=color, label=f"{variant}: {factor:.3f}")
        right.axhline(fk3, color=color, linestyle=":", linewidth=.75)
    right.set(xlabel=r"Smooth-triangle filling $f$ (%)",
              ylabel="Largest relative gap width (%)", title="(b) Paired refinement")
    right.legend(loc="upper left")
    finish(fig, (left, right), output, "area_estimate_robustness")


def energy_convergence(output):
    modes = rows("fields/recess_by_mode.csv")
    refinement = rows("fields/recess_refinement.csv")
    palette = {1: ".39", 4: "#238066"}
    fig, (left, right) = plt.subplots(1, 2, figsize=(8.8, 4.15), layout="constrained")
    for m in (1,4):
        pairs = [(float(r["volume"]), float(r["green"])) for r in modes
                 if int(r["best_mesh"]) == m]
        left.loglog(*np.array(pairs).T, "o", color=palette[m],
                    markersize=4.1, alpha=.68, label=rf"$m={m}$")
    bounds = [2e-5, 2.2e-2]
    left.loglog(bounds, bounds, "--", color=".45", linewidth=.9)
    left.set(xlim=bounds, ylim=bounds, xlabel=r"Volume estimate of $q_{\rm EM}$",
             ylabel=r"Green-identity estimate of $q_{\rm EM}$",
             title="(a) Recess energy across 214 modes")
    left.legend()
    spreads = []
    geometries = sorted({r["geometry"] for r in modes if r["best_mesh"] == "4"})
    for geometry in geometries:
        levels = [sorted((r for r in refinement if r["geometry"] == geometry
                          and int(r["boundary_factor"]) == m),
                         key=lambda r: float(r["omega_r"])) for m in (1,2,4)]
        assert all(len(level) == 7 for level in levels), geometry
        for mode in zip(*levels):
            volume = np.array([float(r["qEM_h8"]) for r in mode])
            green = np.array([float(r["qEM_green"]) for r in mode])
            difference = 100 * np.abs(volume-green)/(volume+green)
            right.plot([1,2,4], difference, "o-", color=palette[4],
                       markersize=3, linewidth=.8, alpha=.40)
            spreads.append(difference)
    right.plot([1,2,4], np.median(spreads, axis=0), "o-",
               color="#183a40", linewidth=1.9, label="Median")
    right.set(xticks=(1,2,4), xlabel=r"Boundary refinement $m$",
              ylabel="Relative half-difference (%)", title="(b) Refinement of 21 modes")
    right.legend()
    finish(fig, (left, right), output, "recess_energy_convergence")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    apply_style()
    area_robustness(args.output)
    energy_convergence(args.output)

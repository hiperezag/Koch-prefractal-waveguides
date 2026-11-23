"""Render Fig. 10 and Figs. S3--S5 from the archived intensity arrays."""

import argparse
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.path import Path as Polygon
from matplotlib.patheffects import withStroke
import numpy as np
from scipy.spatial import ConvexHull

from figure_style import ALLOWED_COLOR, FIELD_CMAP, apply_style

FIELDS = Path(__file__).resolve().parents[1] / "fields"


def read_csv(name):
    with (FIELDS / name).open(newline="") as source:
        return list(csv.DictReader(source))


def notation(value):
    if abs(value) < 1e-20:
        return "0"
    mantissa, power = f"{value:.1e}".split("e")
    return rf"{mantissa}\times10^{{{int(power)}}}"


def weights(record):
    return rf"$q_I={notation(float(record['qI']))},\quad q_{{\rm EM}}={notation(float(record['qEM']))}$"


def field_panel(ax, maps, key, grid, logarithmic=False):
    polygon = maps[f"{key}_polygon"]
    hull = polygon[ConvexHull(polygon).vertices]
    hull = np.vstack([hull, hull[0]])
    original = maps[f"{key}_I"]
    conductor = Polygon(polygon).contains_points(grid).reshape(original.shape)
    intensity = np.ma.array(original, mask=conductor)
    normalized = intensity / intensity.max()
    cmap = plt.get_cmap(FIELD_CMAP).with_extremes(bad="white")
    scale = {"norm": LogNorm(vmin=1e-7, vmax=1)} if logarithmic else {"vmin": 0, "vmax": 1}
    image = ax.imshow(normalized, origin="lower", extent=[0, 2*np.pi, 0, np.pi],
                      interpolation="nearest", cmap=cmap, **scale)
    ax.plot(*polygon.T, color=".25", linewidth=.45)
    line, = ax.plot(*hull.T, color="white", linestyle="--", linewidth=.85)
    line.set_path_effects([withStroke(linewidth=1.35, foreground=".20")])
    if logarithmic:
        ax.set_xlim(hull[:, 0].min() - .06, hull[:, 0].max() + .06)
        ax.set_ylim(hull[:, 1].min() - .06, hull[:, 1].max() + .06)
    ax.set(xticks=[], yticks=[])
    return image


def save(fig, output, name):
    fig.savefig(output / f"{name}.pdf", dpi=500, bbox_inches="tight")
    plt.close(fig)


def pocket_modes(maps, grid, index, records, widths, output, outward):
    shape = "externa_triangulo"
    inward = "interna_triangulo_FK3_igualn"
    if outward:
        cases = [(f"{shape}_FK0", "Smooth FK0"), (f"{shape}_FK1", "Outward FK1"),
                 (f"{shape}_FK3", "Outward FK3")]
        modes, columns = (1, 2, 4), 4
    else:
        cases = [(inward, "Inward FK3"), (f"{inward}_envolvente", "Convex hull"),
                 (f"{shape}_FK0", "Same-area smooth FK0")]
        modes, columns = (1, 2, 3, 4), 3
    fig, axes = plt.subplots(len(modes), columns, figsize=(11.1 if outward else 9.9, 2.08*len(modes)),
                             gridspec_kw={"hspace": .60, "wspace": .23})
    fig.subplots_adjust(left=.075, right=.90, top=.84, bottom=.12)
    profile_colors = [".28", ALLOWED_COLOR, "#b95824"]
    for i, mode in enumerate(modes):
        for j, (geometry, name) in enumerate(cases):
            key = index["pocket", mode, geometry]
            record = records[key]
            ax = axes[i,j]
            image = field_panel(ax, maps, key, grid)
            title = rf"$\omega_r={float(record['omega_r']):.4f},\ W={widths[geometry, mode]:.4f}$"
            ax.set_title((name + "\n" if i == 0 else "") + title, fontsize=12, pad=5)
            ax.text(.5, -.11, weights(record), transform=ax.transAxes, ha="center",
                    va="top", fontsize=11.5)
            if j == 0:
                ax.set_ylabel(f"Mode {mode}", labelpad=9)
            if outward:
                profile = maps[f"{key}_I"].sum(axis=0)
                axes[i,3].plot(maps["x"], profile/profile.max(), color=profile_colors[j],
                               label=name, linewidth=1.3)
        if outward:
            ax = axes[i,3]
            ax.set(xlim=(0, 2*np.pi), ylim=(-.02, 1.07), yticks=(0, 1))
            ax.set_xticks([0, np.pi, 2*np.pi], ["0", r"$\pi$", r"$2\pi$"])
            ax.grid()
            if i == 0:
                ax.set_title(r"$P(x)/\max P$")
            if i == len(modes)-1:
                ax.set_xlabel("$x$")
    if outward:
        handles, labels = axes[0,3].get_legend_handles_labels()
        fig.legend(handles, labels, ncol=3, loc="lower center",
                   bbox_to_anchor=(.49, .008))
    bar = fig.colorbar(image, cax=fig.add_axes([.925, .19, .017, .58]), ticks=[0,.5,1])
    bar.set_label(r"$I/\max I$")
    kind = "outward" if outward else "inward"
    save(fig, output, f"{kind}_triangle_pocket_modes")


def recess_modes(maps, grid, index, records, correlations, output, case):
    inward = case == "interna"
    headers = (["Smooth FK0", "Inward FK3", "Same-area twin", "Notches (log)"] if inward
               else ["Convex hull", "Outward FK3", "Same-area twin", "Valleys (log)"])
    fig, axes = plt.subplots(3, 4, figsize=(11.4, 6.7),
                             gridspec_kw={"hspace": .68, "wspace": .18})
    fig.subplots_adjust(left=.065, right=.91, top=.82, bottom=.19)
    for i, mode in enumerate((1,2,3)):
        row = correlations[case, mode]
        mismatches = [0., row["one_minus_rho_fractal"], row["one_minus_rho_twin"]]
        keys = [index[case, mode, role] for role in ("ref", "fractal", "contraste")]
        for j, (key, mismatch) in enumerate(zip(keys, mismatches)):
            ax = axes[i,j]
            linear_image = field_panel(ax, maps, key, grid)
            frequency = float(records[key]["omega_r"])
            title = rf"$\omega_r={frequency:.4f}$"
            ax.set_title((headers[j] + "\n" if i == 0 else "") + title, pad=5)
            ax.text(.5, -.12, rf"$1-\rho={notation(mismatch)}$", transform=ax.transAxes,
                    va="top", ha="center")
            if j == 0:
                ax.set_ylabel(f"Mode {mode}", labelpad=8)
        ax = axes[i,3]
        key = keys[1]
        log_image = field_panel(ax, maps, key, grid, logarithmic=True)
        record = records[key]
        ax.set_title((headers[3] + "\n" if i == 0 else "")
                     + rf"$q_I={notation(float(record['qI']))}$", pad=5)
        ax.text(.5, -.12, rf"$q_{{\rm EM}}={notation(float(record['qEM']))}$",
                transform=ax.transAxes, va="top", ha="center")
    bar = fig.colorbar(linear_image, cax=fig.add_axes([.13, .075, .45, .02]),
                       orientation="horizontal", ticks=[0, .5, 1])
    bar.set_label(r"$I/\max I$ (linear panels)")
    fig.colorbar(log_image, cax=fig.add_axes([.935, .23, .014, .57]),
                 ticks=[1e-7, 1e-5, 1e-3, 1e-1, 1])
    save(fig, output, "inward_triangle_recess_fields" if inward else "outward_triangle_hull_fields")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    apply_style()
    index = {(r["comparison"], int(r["mode"]), r["role"]): r["field_key"]
             for r in read_csv("manuscript_map_index.csv")}
    records = {r["key"]: r for r in read_csv("manuscript_maps.csv")}
    widths = {(r["geometria"], int(r["modo"])): float(r["W"])
              for r in read_csv("map_bandwidths_and_controls.csv")}
    correlations = {(r["case"], r["n"]): r
                    for r in json.loads((FIELDS/"intensity_correlations.json").read_text())}
    with np.load(FIELDS/"manuscript_maps.npz") as maps:
        x, y = np.meshgrid(maps["x"], maps["y"])
        grid = np.column_stack([x.ravel(), y.ravel()])
        pocket_modes(maps, grid, index, records, widths, args.output, outward=False)
        pocket_modes(maps, grid, index, records, widths, args.output, outward=True)
        for case in ("interna", "externa_env"):
            recess_modes(maps, grid, index, records, correlations, args.output, case)

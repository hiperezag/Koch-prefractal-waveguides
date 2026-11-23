"""Plot Figs. 6--8 and S2 using existing sweep, area, and FEM tables."""

import argparse
from collections import defaultdict
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

from figure_style import ALLOWED_COLOR, SHAPE_STYLE, apply_style


ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT/"gaps/shape_and_iteration_tables"
SHAPES = {
    "triangulo": "triangle",
    "cuadrado": "square",
    "pentagono": "pentagon",
    "hexagono": "hexagon",
    "rectangulo_horizontal": "horizontal_rectangle",
    "rectangulo_vertical": "vertical_rectangle",
}
ORDER_STYLE = {
    "FK0": ("#8a8a8a", "o", "--"),
    "FK1": ("#86b6ef", "s", "-"),
    "FK2": ("#3987e5", "^", "-"),
    "FK3": ("#184f95", "D", "-"),
}


def dat_rows(path):
    with path.open(encoding="utf-8") as stream:
        for line in stream:
            line = line.strip()
            if line and not line.startswith("#"):
                yield line.split()


def save(fig, output, file_name):
    fig.savefig(output/f"{file_name}.pdf", bbox_inches="tight")
    plt.close(fig)


def open_axes(ax):
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid()


def equal_filling(output):
    curves = defaultdict(list)
    for sense, shape, order, filling, _, relative_gap, _ in dat_rows(TABLES/"curvas_igual_n_TE.dat"):
        curves[sense, SHAPES[shape], order].append((float(filling), float(relative_gap)))
    fig, axes = plt.subplots(2, 6, figsize=(11.8, 5.6), sharex=True, sharey=True)
    for i, (sense, growth) in enumerate((("externa", "Outward"), ("interna", "Inward"))):
        for j, (shape, (label, _, _)) in enumerate(SHAPE_STYLE.items()):
            ax = axes[i,j]
            for order, (color, marker, line) in ORDER_STYLE.items():
                points = curves["externa" if order == "FK0" else sense, shape, order]
                if points:
                    ax.plot(*zip(*sorted(points)), color=color, marker=marker, linestyle=line,
                            linewidth=1.5, markersize=5.6, markeredgecolor="white",
                            markeredgewidth=.6, zorder=1 if order == "FK0" else 2)
            if i == 0:
                ax.set_title(label)
            else:
                ax.set_xlabel(r"$f$ (%)")
            if j == 0:
                ax.set_ylabel(growth + "\n" + r"$\Delta\omega/\omega_c$ (%)")
            open_axes(ax)
    handles = [Line2D([], [], color=color, marker=marker, linestyle=line,
                      label="FK0 (smooth)" if order == "FK0" else order)
               for order, (color, marker, line) in ORDER_STYLE.items()]
    fig.legend(handles=handles, ncol=4, loc="upper center", bbox_to_anchor=(.5, 1.))
    fig.tight_layout(rect=(0,0,1,.92))
    save(fig, output, "equal_filling_gap_comparison")


def aperture(output):
    curves = defaultdict(list)
    for shape, *values in dat_rows(TABLES/"barrido_FK0_TE.dat"):
        curves[SHAPES[shape]].append([float(v) for v in values])
    fig, axes = plt.subplots(2, 3, figsize=(10.1, 5.65), sharey="row")
    x_axes = [(1, r"$f$ (%)"), (2, r"Distance to plates $d_{\rm p}$"),
              (3, r"Gap between neighbors $d_{\rm n}$")]
    y_axes = [(4, r"$\Sigma$"), (5, r"$\Delta\omega/\omega_c$ (%)")]
    for i, (y, ylabel) in enumerate(y_axes):
        for j, (x, xlabel) in enumerate(x_axes):
            ax = axes[i,j]
            for shape, (_, color, marker) in SHAPE_STYLE.items():
                table = np.array(curves[shape])
                ax.plot(table[:,x], table[:,y], color=color, linewidth=1.3)
                ax.plot(table[::6,x], table[::6,y], marker, color=color, markersize=5.3,
                        markeredgecolor="white", markeredgewidth=.7)
            if i == 1:
                ax.set_xlabel(xlabel)
            if j == 0:
                ax.set_ylabel(ylabel)
            open_axes(ax)
    handles = [Line2D([], [], color=color, marker=marker, label=label)
               for label, color, marker in SHAPE_STYLE.values()]
    fig.legend(handles=handles, ncol=6, loc="upper center", bbox_to_anchor=(.5, 1.))
    fig.tight_layout(rect=(0,0,1,.91))
    save(fig, output, "filling_and_channel_aperture")


def area_factors(output):
    with (TABLES/"ahorro_consistente_TE.csv").open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    fig, axes = plt.subplots(1, 2, figsize=(9.7, 5.35), sharex=True, sharey=True)
    limits = (.95, 2.35)
    ticks = np.arange(1.,2.21,.2)
    markers = {"FK1": "o", "FK2": "s", "FK3": "D"}
    shapes = ("triangle", "horizontal_rectangle", "vertical_rectangle",
              "square", "pentagon", "hexagon")
    for ax, filling, letter in zip(axes, ("2.1937", "5.0661"), "ab"):
        ax.plot(limits, limits, color=".35", linewidth=1.2)
        for shape in shapes:
            _, color, _ = SHAPE_STYLE[shape]
            selected = sorted((r for r in rows if r["sense"] == "interna" and
                              SHAPES[r["shape"]] == shape and r["f"] == filling and
                              np.isfinite(float(r["reference"]))), key=lambda r: r["order"])
            ax.plot([float(r["area_ratio"]) for r in selected],
                    [float(r["reference"]) for r in selected],
                    color=color, linewidth=1, alpha=.65)
            for row in selected:
                x, y, lo, hi = (float(row[key]) for key in ("area_ratio", "reference", "lo", "hi"))
                ax.errorbar(x, y, yerr=[[y-lo],[hi-y]], fmt="none", color=color,
                            capsize=2.5, linewidth=.9, zorder=2)
                order = row["order"]
                ax.plot(x, y, markers[order], markersize=7.2 if order == "FK3" else 6.4,
                        markeredgecolor=color, markeredgewidth=1.2,
                        markerfacecolor=color if order == "FK3" else "white", zorder=3)
        ax.set(xlim=limits, ylim=limits, xticks=ticks, yticks=ticks,
               xlabel="Hull-to-conductor area ratio",
               title=rf"({letter}) $f={float(filling):.2f}\%$", aspect="equal")
        ax.set_axisbelow(True)
        open_axes(ax)
    axes[0].set_ylabel(r"Material-saving factor $f_{\rm eq}/f$")
    shape_keys = [Line2D([], [], color=SHAPE_STYLE[s][1], marker="D", label=SHAPE_STYLE[s][0])
                  for s in shapes]
    order_keys = [Line2D([], [], color=".35", label=r"Area estimate ($y=x$)")]
    order_keys.extend(Line2D([], [], linestyle="none", marker=mark, markeredgecolor=".25",
                            markerfacecolor=".25" if order == "FK3" else "white", label=order)
                      for order, mark in markers.items())
    fig.legend(handles=shape_keys, ncol=6, loc="lower center", bbox_to_anchor=(.5,.074))
    fig.legend(handles=order_keys, ncol=4, loc="lower center", bbox_to_anchor=(.5,.014))
    fig.tight_layout(rect=(0,.185,1,1), w_pad=1.3)
    save(fig, output, "equivalent_area_factors")


def paired_fem(output):
    rows = list(dat_rows(ROOT/"validation/fem_pares_TE.dat"))
    assert len(rows) == 6 and all(len(row) == 16 for row in rows)
    labels = [SHAPE_STYLE[SHAPES[row[0]]][0] for row in rows]
    y = np.arange(len(rows))[::-1]
    fig, axes = plt.subplots(1, 2, figsize=(8.8, 3.65), sharey=True, layout="constrained")
    for ax, first, reference in zip(axes, (10,13), ("convex hull", "smooth twin")):
        iem, fem = [[float(row[col]) for row in rows] for col in (first, first+2)]
        ax.axvline(0, color=".42", linewidth=.8, linestyle="--")
        ax.plot(fem, y, "o", markersize=7, markerfacecolor="none", markeredgecolor="black",
                label=r"FEM $h=0.04$", zorder=3)
        ax.plot(iem, y, ".", markersize=9, color=ALLOWED_COLOR, label="IEM", zorder=2)
        ax.set(xlabel="Relative difference (%)",
               title=rf"$\Sigma$ of FK3 relative to its {reference}", ylim=(-.5,5.5))
        ax.grid()
    axes[0].set_yticks(y, labels)
    axes[1].legend(loc="upper right")
    save(fig, output, "finite_element_hull_validation")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    directory = parser.parse_args().output
    directory.mkdir(parents=True, exist_ok=True)
    apply_style()
    for plot in (equal_filling, aperture, area_factors, paired_fem):
        plot(directory)

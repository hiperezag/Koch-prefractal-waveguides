"""Render Figs. 3--5, 13 and S1 using the archived band data."""

import argparse
import csv
import io
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.ticker import StrMethodFormatter
import numpy as np

from figure_style import ALLOWED_COLOR, FLAT_COLOR, GAP_ALPHA, GAP_COLOR, apply_style

ROOT = Path(__file__).resolve().parent.parent
BANDS = ROOT / "bands"
sys.path.insert(0, str(ROOT / "gaps"))
from iem_bandas import cargar_bandas, detectar_gaps_global, detectar_modos_discretos


def read_csv(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def save(figure, output, name):
    figure.savefig(output / f"{name}.pdf", bbox_inches="tight")
    plt.close(figure)


def fixed_side(output):
    gaps = [row.split() for row in (BANDS / "inventory/prohibidas_L1n2_TE.dat").read_text().splitlines()
            if row.strip() and not row.startswith("#")]
    bands = [row for row in read_csv(BANDS / "inventory/classified_bands.csv")
             if row["ensemble"] == "fixed_side" and row["flat"] == "True"]
    shapes = [("triangulo", "Triangle"), ("cuadrado", "Square"),
              ("pentagono", "Pentagon"), ("hexagono", "Hexagon")]
    fig, axes = plt.subplots(2, 4, figsize=(10.2, 6), sharex=True, sharey=True,
                             layout="constrained")
    for i, (growth, row_title) in enumerate([("externa", "Outward"), ("interna", "Inward")]):
        for j, (shape, title) in enumerate(shapes):
            ax = axes[i, j]
            for row in gaps:
                if row[:2] != [growth, shape] or (growth == "interna" and row[2] == "FK0"):
                    continue
                order, lo, hi = int(row[2][2:]), float(row[4]), float(row[5])
                color = "#a3a5a9" if "1" in row[9:11] else GAP_COLOR
                ax.add_patch(Rectangle((order - .29, lo), .58, hi - lo,
                                       facecolor=color, alpha=GAP_ALPHA,
                                       edgecolor=".25", linewidth=.35))
            for row in bands:
                if (row["growth"], row["shape"]) != (growth, shape):
                    continue
                order, omega = int(row["order"][2:]), float(row["omega0"])
                if shape == "cuadrado" and row["band"] in ("3Odd", "91"):
                    ax.plot(order + .29, omega, marker="D", color=FLAT_COLOR, ms=4.5)
                else:
                    ax.hlines(omega, order - .33, order + .33, color=FLAT_COLOR, linewidth=1.5)
            ax.axhline(2.6, color=".25", linestyle=":", linewidth=.9)
            ax.set(xlim=(-.5, 4.55), ylim=(1, 3), xticks=range(5),
                   xticklabels=[f"FK{n}" for n in range(5)])
            if i == 0:
                ax.set_title(title)
            if j == 0:
                ax.set_ylabel(rf"{row_title} growth, $\omega_r$")
    save(fig, output, "fixed_side_band_gaps")


def triangle_panels(output, growth):
    archived = {
        "inward": [
            ("inward_FK3_L1.csv", "Inward FK3", "inward", 3, 1., "0.99", "1.1176"),
            ("twin_inward_FK3_L1.csv", "Smooth twin", "smooth", 0, .6728, "0.99", "0.8955"),
        ],
        "outward": [
            ("outward_FK3_L1.csv", "Outward FK3", "outward", 3, 1., "3.39", "1.271"),
            ("twin_outward_FK3_L1.dat", "Smooth twin", "smooth", 0, 1.2439, "3.39", "1.267"),
        ],
    }
    cases = [("smooth_L1.csv", "FK0", "smooth", 0, 1., "2.19", "1.134"), *archived[growth]]
    vertices = read_csv(ROOT / "geometry/unit_side_vertices.csv")
    fig = plt.figure(figsize=(8.6, 4.65), layout="constrained")
    fig.get_layout_engine().set(w_pad=.09, wspace=.06)
    layout = fig.add_gridspec(2, 3, height_ratios=[.9, 3.1])
    for j, (data, title, geometry, order, side, filling, gap_sum) in enumerate(cases):
        outline = fig.add_subplot(layout[0, j])
        poly = np.array([(float(r["x"]), float(r["y"])) for r in vertices
                         if r["shape"] == "triangle" and r["growth"] == geometry
                         and int(r["order"]) == order])
        x, y = (poly * side + [np.pi, np.pi / 2]).T
        outline.fill(x, y, color="#364750")
        outline.hlines([0, np.pi], 0, 2*np.pi, colors=".32", linewidth=.8)
        outline.set(xlim=(0, 2*np.pi), ylim=(-.15, np.pi + .15), aspect="equal")
        outline.set_title(rf"{title}, $L={side:.3f}$" + "\n" + rf"$f={filling}\%$")
        outline.axis("off")
        if data.endswith(".csv"):
            points = np.loadtxt(BANDS / data, delimiter=",", skiprows=1)
            stream = io.StringIO()
            np.savetxt(stream, points, fmt="%.16g")
            stream.seek(0)
            k, frequencies = cargar_bandas(stream)
        else:
            k, frequencies = cargar_bandas(BANDS / data)
        ax = fig.add_subplot(layout[1, j])
        for gap in detectar_gaps_global(k, frequencies, .008):
            ax.axhspan(gap["wr_top"], gap["wr_bot"], color=GAP_COLOR, alpha=GAP_ALPHA,
                       linewidth=0, zorder=0)
        flat = {col for col, _ in detectar_modos_discretos(frequencies)}
        for col, values in enumerate(frequencies.T):
            present = np.isfinite(values)
            ax.plot(k[present], values[present], ".", markersize=3.2,
                    color=FLAT_COLOR if col in flat else ALLOWED_COLOR)
        ax.axhline(2.6, color=".16", linestyle=":", linewidth=.9)
        ax.axvline(0, color=".55", linestyle=":", linewidth=.6)
        ax.grid()
        ax.set(xlim=(-.5, .5), ylim=(1, 3), xticks=[-.5, -.25, 0, .25, .5],
               yticks=[1, 1.5, 2, 2.5, 3], xlabel=r"$K_r$", title=rf"$\Sigma={gap_sum}$")
        if j == 0:
            ax.set_ylabel(r"$\omega_r$")
        else:
            ax.tick_params(labelleft=False)
    save(fig, output, f"{growth}_triangle_band_structure")


def square_parity(output):
    rows = read_csv(BANDS / "square_parity.csv")
    fig, axes = plt.subplots(2, 2, figsize=(8.8, 6), layout="constrained")
    colors = {"Even": ALLOWED_COLOR, "Odd": "#b95824"}
    for ax, window in zip(axes.flat, ("FK2_b5", "FK3_b3", "FK3_b5", "FK4_b3")):
        ax.set_title(f"{window[:3]}, band-{window[-1]} window")
        for parity, band in [("Even", 3), ("Odd", 91)]:
            color = colors[parity]
            if window == "FK4_b3":
                file = BANDS / "fk3_fk4/traces" / f"L1n2_externa_cuadrado_FK4_b{band}.dat"
                k, w = np.loadtxt(file, usecols=(0, 1), unpack=True)
                start = 1 if parity == "Odd" else 0
                ax.plot(k[start:], w[start:], "o-", color=color, ms=3.5)
                if start:
                    ax.plot(k[:2], w[:2], "--", color=color)
                    ax.plot(k[0], w[0], "o", mfc="white", mec=color, ms=5)
            else:
                data = [(float(r["K_r"]), float(r["omega_r"])) for r in rows
                        if r["seed_window"] == window
                        and (float(r["parity"]) > 0) == (parity == "Even")]
                k, w = np.array(data).T
                ax.plot(k, w, color=color)
        if window.endswith("b5"):
            ax.axhline(2.6, color=".2", linestyle=":", linewidth=.9)
        ax.set(xlabel=r"$K_r$", ylabel=r"$\omega_r$", xlim=(-.015, .515))
        ax.yaxis.set_major_formatter(StrMethodFormatter("{x:.3f}"))
        ax.grid()
    axes[0, 1].legend(handles=[Line2D([], [], color=color, label=parity)
                              for parity, color in colors.items()])
    save(fig, output, "square_mode_parity")


def hull_bandwidths(output):
    rows = read_csv(BANDS / "inward_fk3_hull/paired_widths.csv")
    fig, (left, right) = plt.subplots(1, 2, figsize=(8.8, 4.15), layout="constrained")
    limits = (1., 1.8)
    ticks = [1., 1.2, 1.4, 1.6, 1.8]
    right.plot(limits, limits, "--", color=".4", linewidth=.9)
    for level, color, marker in [("2.19", ALLOWED_COLOR, "o"), ("5.07", "#b95824", "s"),
                                 ("8.72", "#238066", "^")]:
        group = [r for r in rows if r["filling_pct"] == level and r["resolved"] == "True"]
        width = [float(r["hull_width_m1"]) for r in group]
        standard = [float(r["ratio_m1"]) for r in group]
        refined = [float(r["ratio_m2"]) for r in group]
        style = dict(linestyle="none", marker=marker, color=color, markersize=5.4, alpha=.85)
        left.semilogx(width, standard, label=rf"$f={level}\%$", **style)
        right.plot(standard, refined, **style)
    left.set(xlabel=r"Hull bandwidth $W_{\rm hull}$",
             ylabel=r"$W_{\rm FK3}/W_{\rm hull}$", title="(a) Standard boundary",
             ylim=limits, yticks=ticks)
    right.set(xlim=limits, ylim=limits, xticks=ticks, yticks=ticks,
              xlabel="Standard ratio ($m=1$)", ylabel="Refined ratio ($m=2$)",
              title="(b) Boundary refinement")
    left.legend()
    for ax in (left, right):
        ax.grid()
    save(fig, output, "prefractal_hull_bandwidth_ratios")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    apply_style()
    fixed_side(args.output)
    triangle_panels(args.output, "inward")
    triangle_panels(args.output, "outward")
    square_parity(args.output)
    hull_bandwidths(args.output)

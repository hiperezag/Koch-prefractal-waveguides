"""Draw Fig. 1 from the archived vertices, following figura_geometria.py."""

import argparse
import csv
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Rectangle
import numpy as np
from scipy.spatial import ConvexHull

from figure_style import SHAPE_STYLE, apply_style


ROOT = Path(__file__).resolve().parents[1]
PERIOD, PLATE_SPACING = 2*np.pi, np.pi
CONDUCTOR = ".45"
INK = ".15"
HULL_COLOR = "#2a78d6"
HULL_LINE = (0, (4, 2.5))


def draw_inclusion(ax, points, show_hull=False, opacity=1):
    ax.add_patch(Polygon(points, closed=True, facecolor=CONDUCTOR, edgecolor="none",
                         alpha=opacity, zorder=2))
    if show_hull:
        hull = points[ConvexHull(points).vertices]
        ax.add_patch(Polygon(hull, closed=True, fill=False, edgecolor=HULL_COLOR,
                             linewidth=1.2, linestyle=HULL_LINE, zorder=3))


def dimension(ax, start, end, label, offset):
    ax.annotate("", end, start, arrowprops={"arrowstyle": "<->", "color": INK,
                  "linewidth": .9, "shrinkA": 0, "shrinkB": 0})
    center = (np.array(start) + np.array(end))/2 + offset
    ax.text(*center, label, ha="center", va="center", color=INK)


def guide_panel(ax, vertices):
    period, b = PERIOD, PLATE_SPACING
    x0, x1, thickness = -.62*period, 1.62*period, .22
    for y in (-thickness, b):
        ax.add_patch(Rectangle((x0, y), x1-x0, thickness, facecolor=".82",
                               edgecolor="none", hatch="////"))
    for y in (0, b):
        ax.plot([x0, x1], [y, y], color=INK, linewidth=1.6)
    polygon = vertices["triangle", "outward", 2]
    for cell in (-1, 0, 1):
        draw_inclusion(ax, polygon + [(cell+.5)*period, b/2],
                       opacity=1 if cell == 0 else .45)
    for x in (0, period):
        ax.plot([x, x], [0, b], color=INK, linewidth=1, linestyle=(0, (4, 3)))
    dimension(ax, (0, .42), (period, .42), r"$P=2\pi$", [0, .24])
    dimension(ax, (period+.45, 0), (period+.45, b), r"$b=\pi$", [.42, 0])
    ax.text(.5*period, b-.36, "unit cell", ha="center", va="center", color=INK)
    ax.text(x1, b+thickness+.14, "perfectly conducting plate", ha="right", color=INK)
    origin = (x1+.15, -thickness-.05)
    ax.annotate("", (origin[0]+.6, origin[1]), origin,
                arrowprops={"arrowstyle": "->", "color": INK, "linewidth": .9})
    ax.annotate("", (origin[0], origin[1]+.6), origin,
                arrowprops={"arrowstyle": "->", "color": INK, "linewidth": .9})
    ax.text(origin[0]+.7, origin[1], "$x$", va="center", color=INK)
    ax.text(origin[0], origin[1]+.73, "$y$", ha="center", color=INK)
    ax.set(xlim=(x0, x1+1.12), ylim=(-thickness-.2, b+thickness+.55))


def triangle_orders(ax, vertices):
    spacing, height = 1.45, 1.38
    for i, growth in enumerate(("outward", "inward")):
        y = -i*height
        ax.text(-.95, y, growth, ha="right", va="center", color=INK)
        for order in range(4):
            polygon = vertices["triangle", growth if order else "smooth", order]
            offset = (polygon[:,1].max()+polygon[:,1].min())/2 if i else 0
            x = order*spacing
            draw_inclusion(ax, polygon + [x, y-offset], show_hull=True)
            if i == 0:
                ax.text(x, .86, f"FK{order}", ha="center", va="bottom", color=INK)
    ax.plot([], [], color=HULL_COLOR, linewidth=1.2, linestyle=HULL_LINE,
            label="convex hull")
    ax.legend(loc="upper center", bbox_to_anchor=(.5, .03), handlelength=2.2)
    ax.set(xlim=(-1.9, 3*spacing+.8), ylim=(-height-.95, 1.25))


def area(points):
    x, y = points.T
    return .5*abs(np.dot(x, np.roll(y, -1))-np.dot(y, np.roll(x, -1)))


def shapes_panel(ax, vertices):
    spacing, height = 1.5, 1.45
    triangle = vertices["triangle", "smooth", 0]
    for i, growth in enumerate(("outward", "inward")):
        y = -i*height
        ax.text(-1., y, growth, ha="right", va="center", color=INK)
        for j, (shape, (label, _, _)) in enumerate(SHAPE_STYLE.items()):
            base = vertices[shape, "smooth", 0]
            side = np.sqrt(area(triangle)/area(base))
            mid_y = (base[:,1].max()+base[:,1].min())/2
            points = (vertices[shape, growth, 3]-[0, mid_y])*side + [j*spacing, y]
            draw_inclusion(ax, points, show_hull=True)
            if i == 1:
                label = label.lower().replace(". rectangle", ".\nrectangle")
                ax.text(j*spacing, y-.72, label, ha="center", va="top", color=INK)
    ax.set(xlim=(-1.9, 5*spacing+.8), ylim=(-height-1.16, .9))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    output = parser.parse_args().output
    output.mkdir(parents=True, exist_ok=True)
    points = defaultdict(list)
    with (ROOT/"geometry/unit_side_vertices.csv").open(newline="") as stream:
        for row in csv.DictReader(stream):
            key = row["shape"], row["growth"], int(row["order"])
            points[key].append((float(row["x"]), float(row["y"])))
    vertices = {key: np.array(values) for key, values in points.items()}
    apply_style()
    plt.rcParams.update({"hatch.color": ".6", "hatch.linewidth": .6})
    fig = plt.figure(figsize=(7.2, 8.9))
    grid = fig.add_gridspec(3, 1, height_ratios=[.95, 1.55, 1.25], hspace=0,
                           left=.03, right=.97, top=.985, bottom=.01)
    panels = (guide_panel, triangle_orders, shapes_panel)
    for i, (panel, label) in enumerate(zip(panels, "abc")):
        ax = fig.add_subplot(grid[i])
        panel(ax, vertices)
        ax.set_aspect("equal")
        ax.axis("off")
        ax.text(0, .98, f"({label})", transform=ax.transAxes, ha="left",
                va="top", color=INK, fontsize=13.5)
    fig.savefig(output/"waveguide_geometry.pdf")
    plt.close(fig)


if __name__ == "__main__":
    main()

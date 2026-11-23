#!/usr/bin/env python3
"""
Plot the IEM and FEM dispersion curves for the four TE validation cases.

    * allowed bands (IEM)      royalblue dots
    * narrow bands (IEM)       green dots (sampled W < 0.01 at all K_r)
    * band gaps                red strips, alpha = 0.85
    * FEM                      open black circles on top of the IEM points
    * Times serif font, STIX math, light grid, dotted line at K_r = 0

Input (validation/dispersion/):
    iem_<i>.dat, fem_<i>.dat   two columns: K_r  omega_r   (i = 1..4)
Use --data-dir to select tables exported by prepare_band_panels.py.

Band gaps are computed from the FEM data with the global-jump criterion:
all frequencies of the panel are sorted and every jump larger than 0.008
(TE threshold) is a gap, unless an allowed band covers it (sampling
artifact). The gaps are clipped at omega_r = 2.6. The script prints the
intervals it shades.

Usage:
    python3 finite_element_band_validation.py
    python3 finite_element_band_validation.py --data-dir TABLES --pdf OUTPUT.pdf
Output (next to this script):
    finite_element_band_validation.pdf, finite_element_band_validation.png

Requires: numpy, matplotlib.
"""
import argparse
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from figure_style import ALLOWED_COLOR, FLAT_COLOR, GAP_ALPHA, GAP_COLOR, apply_style

AQUI = Path(__file__).resolve().parent
DATOS = AQUI.parent / "validation" / "dispersion"

PANELES = [
    (1, "(a) Triangle FK0"),
    (2, "(b) Inward triangle FK3"),
    (3, "(c) Outward vertical rectangle FK3"),
    (4, "(d) Hexagon FK0"),
]
CORTE = 2.6            # Upper limit of the gap analysis (dotted line).
UMBRAL_GAP = 0.008     # TE threshold for the global jump
UMBRAL_PLANA = 0.01    # Width threshold for the green markers.


def cargar(archivo):
    """(K_r unique, matrix (nK, n_bands) with omega sorted per K; NaN if missing)."""
    k, w = np.loadtxt(archivo, unpack=True)
    ku = np.unique(k)
    por_k = [np.sort(w[np.isclose(k, ki)]) for ki in ku]
    m = np.full((len(ku), max(len(v) for v in por_k)), np.nan)
    for i, v in enumerate(por_k):
        m[i, :len(v)] = v
    return ku, m


def bandas_prohibidas(m, techo=CORTE):
    """Gaps (lower, upper) by global jump in the sorted spectrum, clipped at `techo`."""
    w = np.sort(m[~np.isnan(m)])
    w = w[w <= techo]
    gaps = []
    for i in np.where(np.diff(w) > UMBRAL_GAP)[0]:
        lo, hi = w[i], w[i + 1]
        cubierto = any(np.nanmin(c) <= lo + 1e-12 and np.nanmax(c) >= hi - 1e-12
                       for c in m.T)
        if not cubierto:
            gaps.append((lo, hi))
    # last gap: from the top of the last band below the cutoff up to the cutoff
    if len(w) and w[-1] < techo:
        arriba = m[m > techo]
        if arriba.size and not any(np.nanmin(c) <= w[-1] and np.nanmax(c) >= techo
                                   for c in m.T):
            gaps.append((w[-1], techo))
    return gaps


def planas(m):
    """Select complete columns whose sampled bandwidth is below the threshold."""
    return {b for b in range(m.shape[1])
            if not np.isnan(m[:, b]).any()
            and m[:, b].max() - m[:, b].min() < UMBRAL_PLANA}


def main(pdf_output=None, data_dir=DATOS):
    apply_style()
    fig, ejes = plt.subplots(2, 2, figsize=(7.5, 5.4), sharex=True, sharey=True)
    for ax, (i, titulo) in zip(ejes.flat, PANELES):
        k_i, m_i = cargar(data_dir / f"iem_{i}.dat")
        k_f, m_f = cargar(data_dir / f"fem_{i}.dat")

        gaps = bandas_prohibidas(m_f)
        print(f"{titulo}: {len(gaps)} gaps (FEM)")
        for lo, hi in gaps:
            print(f"    {lo:.8f}  {hi:.8f}")
            ax.axhspan(lo, hi, color=GAP_COLOR, alpha=GAP_ALPHA, lw=0, zorder=0)

        flat = planas(m_i)
        for b in range(m_i.shape[1]):
            ok = ~np.isnan(m_i[:, b])
            ax.plot(k_i[ok], m_i[ok, b], ".", ms=4,
                    color=FLAT_COLOR if b in flat else ALLOWED_COLOR,
                    zorder=2 if b in flat else 1)
        for b in range(m_f.shape[1]):
            ok = ~np.isnan(m_f[:, b])
            ax.plot(k_f[ok], m_f[ok, b], "o", ms=2.6, mfc="none", mec="black",
                    mew=0.4, zorder=3)

        ax.axhline(CORTE, color="black", lw=0.8, ls=":", zorder=4)
        ax.axvline(0, color="gray", lw=0.5, ls=":")
        ax.set_title(titulo)
        ax.set_xlim(-0.5, 0.5)
        ax.set_ylim(1.05, 3.0)
        ax.set_xticks([-0.5, -0.25, 0, 0.25, 0.5])
        ax.set_yticks([1.2, 1.6, 2.0, 2.4, 2.8])
        ax.grid(True, alpha=0.3)
    for ax in ejes[1]:
        ax.set_xlabel("$K_r$")
    for ax in ejes[:, 0]:
        ax.set_ylabel(r"$\omega_r$")

    leyenda = [
        Line2D([], [], ls="", marker=".", ms=7, color=ALLOWED_COLOR, label="IEM"),
        Line2D([], [], ls="", marker=".", ms=7, color=FLAT_COLOR,
               label=r"IEM, $W<0.01$"),
        Line2D([], [], ls="", marker="o", ms=4, mfc="none", mec="black",
               mew=0.6, label="FEM"),
        Patch(color=GAP_COLOR, alpha=GAP_ALPHA, lw=0, label="FEM gaps"),
        Line2D([], [], color="black", lw=0.8, ls=":", label="Analysis limit"),
    ]
    fig.legend(handles=leyenda, loc="lower center", ncol=5, frameon=False,
               fontsize=10, handletextpad=0.4, columnspacing=1.2)
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    stem = Path(__file__).stem
    outputs = [pdf_output] if pdf_output else [AQUI / f"{stem}.{ext}" for ext in ("pdf", "png")]
    for output in outputs:
        output.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output, dpi=300)
        print(f"-> {output}")
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument("--data-dir", type=Path, default=DATOS,
                        help="Directory containing the eight IEM and FEM point tables")
    parser.add_argument("--pdf", type=Path, help="Write just this PDF to the given path")
    args = parser.parse_args()
    main(pdf_output=args.pdf, data_dir=args.data_dir)

"""Reproduce the bandwidth--attenuation figure and its descriptive fits."""
import argparse
import csv
from pathlib import Path
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np


HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'figures'))
from figure_style import ALLOWED_COLOR, apply_style


GEOMETRIES = {
    'externa_triangulo_FK0': ('Smooth FK0', 'o', '#262729'),
    'externa_triangulo_FK1': ('Outward FK1–FK3', 's', '#17629d'),
    'externa_triangulo_FK2': ('Outward FK1–FK3', 's', '#17629d'),
    'externa_triangulo_FK3': ('Outward FK1–FK3', 's', '#17629d'),
    'interna_triangulo_FK1_igualn': ('Inward FK1–FK3', '^', '#a65524'),
    'interna_triangulo_FK2_igualn': ('Inward FK1–FK3', '^', '#a65524'),
    'interna_triangulo_FK3_igualn': ('Inward FK1–FK3', '^', '#a65524'),
    'interna_triangulo_FK3_igualn_envolvente': ('Hull of inward FK3', 'D', '#23785c'),
}


def read_modes():
    with (HERE / 'triangle_modes.csv').open(newline='') as stream:
        rows = list(csv.DictReader(stream))
    modes = {(r['geometry'], int(r['mode'])): (float(r['T']), float(r['W'])) for r in rows}
    assert len(modes) == 56
    return np.log([[modes[g, n] for n in range(1, 8)] for g in GEOMETRIES])


def common_fit(x, y):
    """Fit one intercept per mode and a common slope across geometries."""
    dx = x - x.mean(axis=-2, keepdims=True)
    dy = y - y.mean(axis=-2, keepdims=True)
    num = np.sum(dx * dy, axis=(-2, -1))
    den = np.sum(dx**2, axis=(-2, -1))
    with np.errstate(invalid='ignore'):
        return num / den


def draw(output, check_bootstrap=False):
    apply_style()
    x, y = read_modes().transpose(2, 0, 1)
    slope = float(common_fit(x, y))
    intercept = np.mean(y - slope * x, axis=0)
    pooled = np.polyfit(x.ravel(), y.ravel(), 1)
    corr = float(np.corrcoef(x.ravel(), y.ravel())[0, 1])
    assert abs(slope - 0.60) < 0.005 and abs(corr - 0.92) < 0.005
    if check_bootstrap:
        rng = np.random.default_rng(0)
        ix = rng.integers(0, len(x), size=(10000, len(x)))
        samples = common_fit(x[ix], y[ix])
        block_range = np.percentile(samples[np.isfinite(samples)], (2.5, 97.5))
        print(f'Pooled correlation {corr:.6f}, common slope {slope:.6f}.')
        print(f'Central 95% of geometry resamples [{block_range[0]:.6f}, {block_range[1]:.6f}].')

    fig, (left, right) = plt.subplots(1, 2, figsize=(9.55, 4.65), sharey=True, layout='constrained')
    for row, (group, marker, color) in enumerate(GEOMETRIES.values()):
        left.loglog(np.exp(x[row]), np.exp(y[row]), marker, ms=5.3, mfc='none' if marker == '^' else color,
                    mec=color, ls='none', label=group)
    grid = np.linspace(x.min(), x.max(), 80)
    left.loglog(np.exp(grid), np.exp(np.polyval(pooled, grid)), ls='--', c='.32', lw=1.1, label='Pooled fit')
    for xm, c in zip(x.T, intercept):
        v = np.linspace(xm.min(), xm.max(), 50)
        left.loglog(np.exp(v), np.exp(c + slope * v), c='.43', alpha=.72, lw=.85)
    handles, labels = left.get_legend_handles_labels()
    legend = dict(zip(labels, handles))
    legend[f'Common slope $b={slope:.2f}$'] = Line2D([0], [0], c='.43', lw=.9)
    left.legend(list(legend.values()), list(legend), loc='lower right')
    left.set_title('(a) Seven modes in eight triangular geometries')
    left.set_ylabel('Bandwidth $W$')

    smooth = np.loadtxt(HERE / 'diagnostics/acoplamiento_wt_TE.dat')
    counts = []
    for n, color, marker in zip(range(1, 5), [ALLOWED_COLOR, '#b95824', '#238066', '#7656a4'], ['o', 's', '^', 'D']):
        subset = smooth[(smooth[:, 1] == n) & (smooth[:, 3] > 0) & (smooth[:, 6] < 1)]
        counts.append(len(subset))
        u, v = np.log(subset[:, 6]), np.log(subset[:, 3])
        coeff = np.polyfit(u, v, 1)
        grid = np.linspace(u.min(), u.max(), 60)
        right.loglog(np.exp(u), np.exp(v), marker, ls='none', ms=3.3, c=color, alpha=.83)
        right.loglog(np.exp(grid), np.exp(np.polyval(coeff, grid)), c=color, lw=1.1,
                     label=f'Band {n}, $b={coeff[0]:.2f}$ ($n={len(subset)}$)')
        if n == 2:
            t0, w0 = u.min(), v[u.argmin()]
    assert counts == [50, 50, 48, 37]
    grid = np.linspace(t0, 0, 60)
    right.loglog(np.exp(grid), np.exp(w0 + 0.5 * (grid - t0)), ls='--', c='.40',
                 lw=.9, label=r'Tight-binding comparison, $b=1/2$')
    right.legend(loc='lower right')
    right.set_title(r'(b) Smooth triangle with $T<1$')
    for ax in (left, right):
        ax.axhline(.01, c='.5', ls=':', lw=.85)
        ax.set_xlabel(r'Attenuation proxy $T=e^{-2S_{\rm lower}}+e^{-2S_{\rm upper}}$')
        ax.grid(color='.91', lw=.6, which='major')
        ax.spines[['top', 'right']].set_visible(False)
    fig.savefig(output, bbox_inches='tight')
    plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    parser.add_argument('--check-bootstrap', action='store_true',
                        help='Repeat the archived resampling check using seed 0')
    args = parser.parse_args()
    draw(args.output, args.check_bootstrap)

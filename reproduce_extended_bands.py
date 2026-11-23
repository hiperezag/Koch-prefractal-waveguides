"""Check the branch classification, FK3/FK4 gap edges, and inward FK3/hull ratios."""
from bisect import bisect_right
from collections import Counter, defaultdict
import csv
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parent
BANDS = ROOT / 'bands'
sys.path.insert(0, str(ROOT / 'gaps'))
from iem_bandas import cargar_bandas


def read_csv(path):
    with path.open(newline='') as stream:
        return list(csv.DictReader(stream))


def check_inventory():
    records = read_csv(BANDS / 'inventory/classified_bands.csv')
    fixed = [r for r in records if r['ensemble'] == 'fixed_side']
    equal = [r for r in records if r['ensemble'] == 'equal_filling']
    joined = fixed + [r for r in equal if (r['shape'], r['order']) != ('triangulo', 'FK0')]
    for r in records:
        flat = float(r['width']) < .01 and float(r['slope']) < .08
        assert flat == (r['flat'] == 'True')
    bins = Counter(bisect_right((.001, .006, .01, .015), float(r['width'])) for r in joined)
    assert len(joined) == 479 and [bins[i] for i in range(5)] == [77, 87, 19, 37, 259]
    assert sum(r['flat'] == 'True' for r in fixed) == 129
    assert sum(r['flat'] == 'True' for r in equal) == 55
    print(f'{len(joined)} branches and {sum(r["flat"] == "True" for r in fixed)} fixed-side flat bands.')


def traced_gaps(root, shape, order, kind):
    suffix = 'igualn' if kind == 'prefractal' else 'igualn_envolvente'
    tag = f'interna_{shape}_FK{order}_{suffix}'
    target_k = np.linspace(0., .5, 11)
    k, sweep = cargar_bandas(root / 'sweeps' / f'bandas_TE_{tag}.dat')
    if np.min(k) < 0:
        target_k = -target_k
    samples = np.abs(k[:, None] - target_k).argmin(axis=0)
    assert np.max(np.abs(k[samples] - target_k)) <= .006
    columns = sweep[samples].T.copy()
    matched = set()
    extra = []
    files = sorted((root / 'traces').glob(f'IL{order}_{tag}_b*.dat'))
    assert len(files) == 7
    for file in files:
        trace = np.loadtxt(file)
        assert np.array_equal(trace[:, 0], np.round(np.abs(target_k), 3))
        w = trace[:, 1]
        candidates = [j for j, c in enumerate(columns) if j not in matched and np.isfinite(c[0])]
        closest = min(candidates, key=lambda j: abs(columns[j, 0] - w[0]))
        if abs(columns[closest, 0] - w[0]) > .006:
            extra.append(w)
        else:
            columns[closest] = w
            matched.add(closest)
    intervals = sorted((np.nanmin(c), np.nanmax(c)) for c in [*columns, *extra] if np.isfinite(c).any())
    upper = intervals[0][1]
    gaps = []
    for lo, hi in intervals[1:]:
        if lo - upper > .008 and upper < 2.6:
            gaps.append((upper, min(lo, 2.6)))
        upper = max(upper, hi)
    return np.asarray(gaps)


def check_fk4():
    root = BANDS / 'fk3_fk4'
    summary = read_csv(root / 'fk4_width_controls.csv')
    assert len(summary) == 113 and len({r['geometry'] for r in summary}) == 17
    for r in summary:
        trace = np.loadtxt(root / 'traces' / f"{r['ensemble']}_{r['geometry']}_b{r['band']}.dat")
        assert len(trace) == 11
        checked = trace[1:, 1] if r['endpoint'] == 'overlapping' else trace[:, 1]
        assert abs(np.ptp(checked) - float(r['trace_width'])) < 1e-10
        slope = max(abs(np.gradient(trace[:, 1], trace[:, 0])))
        assert abs(slope - float(r['trace_slope'])) < 1e-10
    physical = [r for r in summary if r['ensemble'] == 'L1n2' and r['physical'] == 'True']
    assert len(physical) == 28 and sum(r['flat'] == 'True' for r in physical) == 18
    assert sum(r['ensemble'] == 'IL4' and r['flat'] == 'True' for r in summary) == 26
    odd = np.loadtxt(root / 'traces/L1n2_externa_cuadrado_FK4_b91.dat')
    assert 0.00127 < np.ptp(odd[:, 1]) < 0.00129
    groups = defaultdict(list)
    for r in read_csv(root / 'gap_edges.csv'):
        groups[r['shape'], r['order'], r['geometry']].append((float(r['lower']), float(r['upper'])))
    assert len(groups) == 24
    for key, pairs in groups.items():
        assert np.allclose(traced_gaps(root, *key), pairs, atol=5e-9, rtol=0)
    print('FK4 widths and all 24 FK3/FK4 gap lists reproduce from archived traces.')


def check_hull_ratios():
    root = BANDS / 'inward_fk3_hull'
    base = {'2.19': 'triangulo', '5.07': 'cuadrado', '8.72': 'pentagono'}
    pairs = read_csv(root / 'paired_widths.csv')
    good = 0
    unresolved = Counter()
    for r in pairs:
        stem = f"{base[r['filling_pct']]}_{r['shape']}"
        widths = []
        for m in (1, 2):
            p = np.loadtxt(root / 'traces' / f"{stem}_prefractal_m{m}_b{r['fk3_band']}.dat")
            h = np.loadtxt(root / 'traces' / f"{stem}_envolvente_m{m}_b{r['hull_band']}.dat")
            w = np.ptp(p[:, 1])
            v = np.ptp(h[:, 1])
            assert abs(w - float(r[f'fk3_width_m{m}'])) < 1e-10
            assert abs(v - float(r[f'hull_width_m{m}'])) < 1e-10
            assert abs(w / v - float(r[f'ratio_m{m}'])) < 1e-10
            widths.extend((w, v))
        resolved = min(widths) >= 1e-5
        assert resolved == (r['resolved'] == 'True')
        if resolved:
            good += 1
            assert min(float(r['ratio_m1']), float(r['ratio_m2'])) > 1.04
        else:
            unresolved[r['filling_pct'], r['shape']] += 1
    assert len(pairs) == 71 and good == 66
    assert unresolved == {('8.72', 'rectangulo_vertical'): 5}
    print('Inward FK3 vs. hull: 66 resolved and 5 unresolved pairs.')


if __name__ == '__main__':
    check_inventory()
    check_fk4()
    check_hull_ratios()

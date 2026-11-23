"""
Shared band-gap and discrete-mode detection.

The base-case analysis, per-run plots, and aggregate sweep analysis use
the same detector (analisis.py, variaciones_L.py, and variaciones_analisis.py).

A gap candidate is a jump above the threshold in the sorted union of
frequencies over all K (0.008 for TE and 0.03 for TM). The coverage
filter discards a candidate crossed by an accepted band segment.
The matching rule links nearby K samples by frequency and accepts a
segment when |Δω| ≤ SLOPE_MAX × ΔK. SLOPE_MAX allows a margin above the
reduced vacuum slope of one.

At 100 K samples a nearly empty guide can leave sampling intervals of
about 0.0101 in a dispersive band. These intervals exceed the TE gap
threshold in some sweeps at L ≤ 0.3. The global-jump and coverage
criteria agree for the densely sampled flat TE bands at L = 1.
The slope bound also rejects long matches near the upper frequency
limit, where a branch can leave the sampled window.
"""
import numpy as np

UMBRAL_PROHIBIDA = {"TE": 0.008, "TM": 0.03}
SLOPE_MAX     = 1.6     # Margin above the reduced vacuum slope |dω_r/dK_r| = 1.


def cargar_bandas(archivo):
    """Return unique K values and sorted frequencies per K, padded with NaN."""
    datos = np.loadtxt(archivo)
    k, omega = datos[:, 0], datos[:, 1]
    k_unicos = np.unique(k)
    bandas_por_k = {}
    for ki, wi in zip(k, omega):
        bandas_por_k.setdefault(ki, []).append(wi)
    max_bandas = max(len(v) for v in bandas_por_k.values())
    bandas_mat = np.full((len(k_unicos), max_bandas), np.nan)
    for i, ki in enumerate(k_unicos):
        omegas = sorted(bandas_por_k[ki])
        bandas_mat[i, :len(omegas)] = omegas
    return k_unicos, bandas_mat


def _segmentos(k_unicos, bandas_mat):
    """Match band segments by nearest frequency at nearby K samples.

    The frequency rank of a mode can change when a sample is missing.
    Nearest-frequency matching with a slope bound reduces long false
    connections in this case. At 100 K samples the one-step bound is
    about 0.016, below the smallest physical gap of 0.0178 in the
    calibration spectra.
    """
    segs = []
    for j in range(len(k_unicos) - 1):
        for r in (1, 2):                 # Two steps bridge one missing sample.
            if j + r >= len(k_unicos):
                continue
            c1 = bandas_mat[j,     ~np.isnan(bandas_mat[j])]
            c2 = bandas_mat[j + r, ~np.isnan(bandas_mat[j + r])]
            if len(c1) == 0 or len(c2) == 0:
                continue
            dk  = k_unicos[j + r] - k_unicos[j]
            los, his = [], []
            for a, b in ((c1, c2), (c2, c1)):
                pos  = np.searchsorted(b, a)
                izq  = b[np.clip(pos - 1, 0, len(b) - 1)]
                der  = b[np.clip(pos,     0, len(b) - 1)]
                vec  = np.where(np.abs(a - izq) <= np.abs(a - der), izq, der)
                ok   = np.abs(a - vec) <= SLOPE_MAX * dk
                los.append(np.minimum(a, vec)[ok])
                his.append(np.maximum(a, vec)[ok])
            segs.append((np.concatenate(los), np.concatenate(his)))
    return segs


def detectar_gaps_global(k_unicos, bandas_mat, umbral=0.03):
    """Find global frequency jumps and exclude intervals covered by a band.

    Each record uses the keys expected by analisis.py. 'wr_top' is the
    lower gap edge at the top of the lower band. 'wr_bot' is the upper
    gap edge at the bottom of the upper band. 'centro' and 'ancho' give
    the gap center and width.
    """
    omega_all    = bandas_mat[~np.isnan(bandas_mat)].flatten()
    omega_sorted = np.sort(omega_all)
    diffs        = np.diff(omega_sorted)
    indices      = np.where(diffs > umbral)[0]
    segs = _segmentos(k_unicos, bandas_mat)
    gaps = []
    for i in indices:
        wb, wt = omega_sorted[i], omega_sorted[i + 1]
        if any(np.any((lo <= wb + 1e-12) & (hi >= wt - 1e-12))
               for lo, hi in segs):
            continue                     # A band crosses this sampling interval.
        gaps.append({'wr_top': float(wb), 'wr_bot': float(wt),
                     'centro': float(0.5 * (wb + wt)),
                     'ancho':  float(wt - wb)})
    return gaps


def detectar_modos_discretos(bandas_mat, umbral=0.01):
    """Identify complete frequency columns with a range below the threshold."""
    n_kr, n_cols = bandas_mat.shape
    modos = []
    for b in range(n_cols):
        col     = bandas_mat[:, b]
        validos = col[~np.isnan(col)]
        if len(validos) == n_kr:
            if validos.max() - validos.min() < umbral:
                modos.append((b, round(float(np.mean(validos)), 15)))
    return modos


UMBRAL_VG = 0.08   # Group-velocity threshold described in detectar_modos_vg.


def _rastrear_banda(k_unicos, bandas_mat, idx0, v0):
    """Trace a branch in both directions from (k_unicos[idx0], v0).

    Use nearest-frequency matching with the slope bound at each step.
    The result is aligned with k_unicos. NaN marks a missing sample or
    a branch that has left the frequency window.
    """
    n_k = len(k_unicos)
    traza = np.full(n_k, np.nan)
    traza[idx0] = v0
    for signo, rango in ((+1, range(idx0, n_k - 1)), (-1, range(idx0, 0, -1))):
        cur = v0
        for i in rango:
            j = i + signo
            col = bandas_mat[j]
            col = col[~np.isnan(col)]
            if len(col) == 0 or np.isnan(cur):
                cur = np.nan
                traza[j] = np.nan
                continue
            dk   = abs(k_unicos[j] - k_unicos[i])
            pos  = np.searchsorted(col, cur)
            vec  = col[max(pos - 1, 0):pos + 1]
            pick = vec[np.argmin(np.abs(vec - cur))]
            cur  = pick if abs(pick - cur) <= SLOPE_MAX * dk else np.nan
            traza[j] = cur
    return traza


def _velocidad_grupo(k_unicos, traza):
    """Estimate dω/dK along a trace with centered finite differences.

    Use one-sided differences at endpoints and next to missing samples.
    """
    vg = np.full_like(traza, np.nan)
    n = len(traza)
    for i in range(n):
        if np.isnan(traza[i]):
            continue
        izq  = i > 0 and not np.isnan(traza[i - 1])
        der  = i < n - 1 and not np.isnan(traza[i + 1])
        if izq and der:
            vg[i] = (traza[i + 1] - traza[i - 1]) / (k_unicos[i + 1] - k_unicos[i - 1])
        elif der:
            vg[i] = (traza[i + 1] - traza[i]) / (k_unicos[i + 1] - k_unicos[i])
        elif izq:
            vg[i] = (traza[i] - traza[i - 1]) / (k_unicos[i] - k_unicos[i - 1])
    return vg


def detectar_modos_vg(k_unicos, bandas_mat, umbral=UMBRAL_VG):
    """Complement the frequency-range criterion with a group-velocity test.

    Sorted frequency columns can change branch at a crossing. Trace
    each branch from the sample nearest K = 0 using nearest-frequency
    matching. Estimate v_g = dω/dK along the trace. Accept complete
    branches whose maximum absolute v_g is below the threshold.

    For the FK3 triangle, the three field-confirmed modes have maximum
    slopes of 0.002–0.065. The next dispersive branch has a maximum
    slope of 0.149. The eight FK0 hexagon candidates have maximum
    slopes at most 0.004. The default threshold is 0.08.

    Each record gives the mean frequency ('omega') and the maximum
    absolute group velocity ('max_vg') of the traced branch.
    """
    idx0 = int(np.argmin(np.abs(k_unicos)))
    col0 = bandas_mat[idx0]
    vistos, modos = set(), []
    for v0 in col0[~np.isnan(col0)]:
        key = round(float(v0), 9)
        if key in vistos:
            continue
        vistos.add(key)
        traza = _rastrear_banda(k_unicos, bandas_mat, idx0, v0)
        if np.any(np.isnan(traza)):
            continue                     # Require a trace at every K sample.
        max_vg = float(np.nanmax(np.abs(_velocidad_grupo(k_unicos, traza))))
        if max_vg < umbral:
            modos.append({'omega': round(float(np.mean(traza)), 15),
                          'max_vg': round(max_vg, 6)})
    return sorted(modos, key=lambda m: m['omega'])

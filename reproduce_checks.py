"""Replay central gap, field-energy, parity, and bandwidth comparisons using the accompanying data."""
from pathlib import Path
import csv
import math
import sys
import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT/'gaps'))
import iem_bandas

def table(name):
    with (ROOT/name).open(newline='') as f:
        return list(csv.DictReader(f))

def quadratic_crossing(filling, gaps, target):
    # Reference estimator: quadratic least squares within one filling percentage point
    # of the first linear crossing. Both the lower and upper shifts use the same rule.
    i = next(i for i in range(len(gaps)-1) if (gaps[i]-target)*(gaps[i+1]-target)<=0)
    center = filling[i] + (target-gaps[i])*(filling[i+1]-filling[i])/(gaps[i+1]-gaps[i])
    keep = abs(filling-center)<=1.
    assert sum(keep)>=5
    coeff = np.polyfit(filling[keep]-center,gaps[keep]-target,2)
    roots = [z.real for z in np.roots(coeff) if abs(z.imag)<1e-10 and abs(z.real)<=1]
    return center + min(roots,key=abs)

def widest_gap(path):
    k, bands = iem_bandas.cargar_bandas(path)
    gaps = iem_bandas.detectar_gaps_global(k,bands,iem_bandas.UMBRAL_PROHIBIDA['TE'])
    return 100*max(g['ancho']/g['centro'] for g in gaps if g['wr_bot']<=2.6)

def main():
    methods = table('gaps/saving_methods.csv')
    principal = next(r for r in methods if r['shape']=='triangle' and r['order']=='3' and float(r['filling_pct'])<3)
    finite=[r for r in methods if math.isfinite(float(r['saving']))]
    assert len(methods)==36 and len(finite)==34 and sum(float(r['saving'])<float(r['area_ratio']) for r in finite)==32
    smooth=[r for r in table('gaps/smooth_curves.csv') if r['shape']=='triangle']
    f=np.array([float(r['filling_pct']) for r in smooth])
    g=np.array([float(r['widest_gap_pct']) for r in smooth])
    f3=float(principal['filling_pct'])
    target=float(principal['target_gap_pct'])
    factors=[quadratic_crossing(f,g,target+shift)/f3 for shift in (0,-.2,.2)]
    # The dense curve is tabulated to 0.001 percentage point; raw paired files have full precision.
    assert np.allclose(factors,[float(principal[k]) for k in ('saving','shifted_low','shifted_high')],atol=.001,rtol=0)
    print('Dense curve: f_eq = %.3f%%, factor = %.3f, shifted interval [%.3f, %.3f]' % (factors[0]*f3,*factors))
    runs=table('gaps/paired_runs.csv'); curves=table('gaps/paired_curves.csv'); paired=table('gaps/paired_factors.csv')
    assert len(runs)==len(curves)==60
    for ref in paired:
        key=ref['variant']
        points=[]
        for r in sorted((r for r in runs if r['variant']==key),key=lambda r:float(r['filling_pct'])):
            width=widest_gap(ROOT/'bands'/'paired'/key/(r['geometry']+'.dat'))
            filling=float(r['filling_pct'])
            expected=next(x for x in curves if x['variant']==key and abs(float(x['filling_pct'])-filling)<1e-7)
            assert abs(width-float(expected['widest_gap_pct']))<=.00051
            if r['geometry']=='FK3': target=width
            else: points.append((filling,width))
        f,g=np.array(points).T
        factor=quadratic_crossing(f,g,target)/f3
        assert len(points)==11 and abs(factor-float(ref['quadratic_factor']))<.000051
        print(f'Paired {key}: {factor:.4f}')
    modes=table('fields/recess_by_mode.csv')
    refined=[r for r in modes if r['best_mesh']=='4']
    assert len(modes)==214 and len(refined)==21
    for r in modes:
        v,gr,mid=[float(r[k]) for k in ('volume','green','qEM_midpoint')]
        assert math.isclose((v+gr)/2,mid,rel_tol=.001)
    for r in refined:
        lo,hi=sorted(float(r[k]) for k in ('volume_standard','green_standard'))
        assert lo-1e-7<=float(r['qEM_midpoint'])<=hi+1e-7
    largest=max(modes,key=lambda r:float(r['qEM_midpoint']))
    assert largest['geometry']=='externa_triangulo_FK1'
    print('Largest central energy estimate: %.3f%%; volume %.3f%%, Green %.3f%%' % tuple(100*float(largest[k]) for k in ('qEM_midpoint','volume','green')))
    square=table('bands/square_parity.csv')
    for seed in sorted(set(r['seed_window'] for r in square)):
        widths=[]
        for branch in (0,1):
            series=[r for r in square if r['seed_window']==seed and int(r['branch'])==branch]
            k=np.array([float(r['K_r']) for r in series])
            w=np.array([float(r['omega_r']) for r in series])
            assert len(w)==51 and np.all(np.diff(k)>0)
            assert len({np.sign(float(r['parity'])) for r in series})==1
            widths.append(np.ptp(w))
            assert (np.ptp(w)<.01 and max(abs(np.gradient(w,k)))<.08)==(seed=='FK3_b3')
        print('Square %s widths: %.5f, %.5f' % (seed,*widths))
    band=table('bands/triangle_widths.csv')
    for r in band:
        assert math.isclose(float(r['FK3_W'])/float(r['hull_W']),float(r['ratio_m1']))
    ratios=[float(r[k]) for r in band for k in ('ratio_m1','ratio_m2','ratio_m4')]
    print('Triangular FK3 to hull width ratios: %.3f to %.3f' % (min(ratios),max(ratios)))
if __name__=='__main__':
    main()

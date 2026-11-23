"""Plot the triangle, the inward FK3 with the same side, and the smooth triangle of equal area."""
from pathlib import Path
import argparse
import csv
import io
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'gaps'))
import iem_bandas

def area(polygon):
    x,y=polygon.T
    return abs(np.dot(x[:-1],y[1:])-np.dot(y[:-1],x[1:]))/2

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('output',type=Path,nargs='?',default=Path('triangle_bands.pdf'))
    args=ap.parse_args()
    with (ROOT/'geometry/unit_side_vertices.csv').open() as f:
        rows=list(csv.DictReader(f))
    def vertices(growth,order):
        return np.array([[float(row['x']),float(row['y'])] for row in rows
                         if row['shape']=='triangle' and row['growth']==growth and int(row['order'])==order])
    smooth,koch=vertices('smooth',0),vertices('inward',3)
    twin=smooth*np.sqrt(area(koch)/area(smooth))
    cases=[('smooth_L1','Full triangular hull',smooth),
           ('inward_FK3_L1','Inward FK3',koch),
           ('twin_inward_FK3_L1','Same-area smooth twin',twin)]
    plt.rcParams.update({'font.family':'DejaVu Sans','mathtext.fontset':'dejavusans','font.size':9,'pdf.fonttype':42,
                         'axes.linewidth':.6,'xtick.major.size':3,'ytick.major.size':3})
    fig=plt.figure(figsize=(8.7,4.3),layout='constrained')
    gs=fig.add_gridspec(2,3,height_ratios=(.8,2.8),hspace=.01,wspace=.03)
    for i,(case,title,poly) in enumerate(cases):
        geom=fig.add_subplot(gs[0,i])
        geom.set(xlim=(0,2*np.pi),ylim=(-.15,np.pi+.15),aspect='equal')
        geom.fill(poly[:,0]+np.pi,poly[:,1]+np.pi/2,color='#384754')
        for y in (0,np.pi):geom.plot([0,2*np.pi],[y,y],color='#647381',lw=.8)
        geom.axis('off')
        geom.set_title(title+'\n'+r'$f=$'+f'{100*area(poly)/(2*np.pi**2):.2f}%',fontsize=10,pad=0)
        points=np.loadtxt(ROOT/'bands'/f'{case}.csv',delimiter=',',skiprows=1)
        # The original detector reads two whitespace-separated columns.
        stream=io.StringIO()
        np.savetxt(stream,points,fmt='%.15g')
        stream.seek(0)
        k,bands=iem_bandas.cargar_bandas(stream)
        gaps=iem_bandas.detectar_gaps_global(k,bands,iem_bandas.UMBRAL_PROHIBIDA['TE'])
        sigma=sum(min(g['wr_bot'],2.6)-g['wr_top'] for g in gaps if g['wr_top']<2.6)
        ax=fig.add_subplot(gs[1,i])
        for gap in gaps:
            ax.axhspan(gap['wr_top'],gap['wr_bot'],facecolor='#c6dfdf',linewidth=0,zorder=1)
        ax.scatter(points[:,0],points[:,1],s=1.5,color='#253964',linewidths=0,zorder=3)
        ax.axhline(2.6,color='#76501b',lw=1.1,ls='--',zorder=4)
        ax.set(xlim=(-.5,.5),ylim=(1.0,3.),xticks=[-.5,0,.5],yticks=[1,1.5,2,2.5,3],xlabel=r'$K_r$')
        if i==0:ax.set_ylabel(r'Reduced frequency $\omega_r$')
        else:ax.tick_params(labelleft=False)
        ax.set_title(r'$\Sigma=$'+f'{sigma:.3f}',fontsize=9,pad=3)
    fig.savefig(args.output,dpi=220,bbox_inches='tight',pad_inches=.02)
    plt.close(fig)
if __name__=='__main__':
    main()

"""Render the six archived intensities with a common linear scale after per-mode normalization."""
from pathlib import Path
import argparse

import matplotlib
matplotlib.use('Agg')
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
from matplotlib.path import Path as Polygon
import numpy as np

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path, nargs='?', default=Path('triangle_field_maps.pdf'))
    args = parser.parse_args()
    plt.rcParams.update({'font.size':9,'font.family':'DejaVu Sans','mathtext.fontset':'dejavusans',
                         'pdf.fonttype':42,'ps.fonttype':42,'axes.linewidth':0.6,
                         'xtick.major.size':2.5,'ytick.major.size':2.5})
    with np.load(Path(__file__).with_name('triangle_maps.npz')) as maps:
        x, y = np.meshgrid(maps['x'],maps['y'])
        xy = np.column_stack((x.ravel(),y.ravel()))
        labels = [('twin','Same-area smooth triangle'),('FK3','Inward FK3'),('hull','Convex hull')]
        fig, axs = plt.subplots(2,3,figsize=(8.2,3.88))
        fig.subplots_adjust(left=.105,right=.985,bottom=.24,top=.845,hspace=.42,wspace=.16)
        cm = plt.get_cmap('viridis').copy()
        cm.set_bad('white')
        for col,(tag,title) in enumerate(labels):
            boundary = maps[tag+'_polygon']
            closed = np.vstack((boundary,boundary[0]))
            solid = Polygon(closed).contains_points(xy).reshape(x.shape)
            for row,mode in enumerate((1,3)):
                ax=axs[row,col]
                values = np.ma.array(maps[f'{tag}_I{mode}'],mask=solid)
                im=ax.imshow(values/values.max(),origin='lower',extent=(0,2*np.pi,0,np.pi),
                             interpolation='nearest',aspect='equal',vmin=0,vmax=1,cmap=cm)
                ax.plot(closed[:,0],closed[:,1],color='#30343b',lw=.50)
                if tag=='FK3':
                    hull=maps['hull_polygon']
                    hull=np.vstack((hull,hull[0]))
                    ax.plot(hull[:,0],hull[:,1],ls='--',color='white',lw=.65,
                            path_effects=[pe.Stroke(linewidth=1.18,foreground='#373739'),pe.Normal()])
                ax.set(xticks=[0,np.pi,2*np.pi],yticks=[0,np.pi],xlim=(0,2*np.pi),ylim=(0,np.pi))
                ax.tick_params(pad=2,labelsize=8)
                ax.set_xticklabels(['0','0.5','1'] if row else [])
                ax.set_yticklabels(['0','1'] if col==0 else [])
                if col==0:
                    ax.set_ylabel(f'Mode {mode}\n'+r'$y/b$',fontsize=9)
                if row:
                    ax.set_xlabel(r'$x/P$',fontsize=9,labelpad=1)
                if not row:ax.text(.5,1.34,title,transform=ax.transAxes,ha='center',va='bottom',fontsize=10)
                ax.set_title(r'$\omega_r=$'+f'{float(maps[f"{tag}_omega{mode}"]):.4f}',fontsize=8.5,pad=3)
        cax=fig.add_axes([.34,.095,.40,.022])
        bar=fig.colorbar(im,cax=cax,orientation='horizontal',ticks=[0,.25,.5,.75,1])
        bar.ax.tick_params(labelsize=8,pad=1)
        bar.set_label(r'Vacuum intensity $|E_z|^2/\max |E_z|^2$',fontsize=9,labelpad=2)
        fig.savefig(args.output,dpi=250,bbox_inches='tight',pad_inches=.015)
        plt.close(fig)
if __name__=='__main__':
    main()

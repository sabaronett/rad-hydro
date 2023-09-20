#!/usr/bin/env python3
#==============================================================================
# vid.py
#
# Create an MPEG-4 animation of the mean radiation energy density.
#
# Author: Stanley A. Baronett
# Created: 2022-09-20
# Updated: 2023-09-20
#==============================================================================
import sys
sys.path.insert(0, '/mnt/home/sbaronett/github/PrincetonUniversity/athena/vis/python')
import athena_read
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib import animation
import numpy as np
from pathlib import Path

# Read data and configure plot
athinput = athena_read.athinput('athinput.beam')
problem_id = athinput['job']['problem_id']
outputs = sorted(list(Path('athdf').glob(problem_id + '.out1.*.athdf')))
athdf = athena_read.athdf(outputs[0])
xv, yv = athdf['x1v'], athdf['x2v']
ts, Ers = [], []
dpi = 236 # 658 x 1075 resolution
fig, ax = plt.subplots(dpi=dpi)
vmin, vmax = 0.05, 2

for output in outputs:
    athdf = athena_read.athdf(output)
    ts.append(athdf['Time'])
    Ers.append(athdf['Er'][0])

# Initialize first frame
clip = np.clip(Ers[0], vmin, vmax)
mesh = ax.pcolormesh(xv, yv, clip, norm=colors.LogNorm(vmin, vmax),
                     cmap='plasma')
cb = fig.colorbar(mesh, label='$E_\mathrm{r}/(a_\mathrm{r}T_0^4)$')
ax.minorticks_on()
ax.set(aspect='equal', title=f'$t={ts[0]:.2f}t_0$', xlabel='$x/L_0$',
       ylabel='$y/L_0$')

def animate(i):
    """
    Update frame.

    Parameters
    ----------
        i : int
            Frame number.
    """
    ax.set_title(f'$t={ts[i]:.2f}t_0$')
    clip = np.clip(Ers[i], vmin, vmax)
    mesh.set_array(clip)
    mesh.set_clim(vmin, vmax)
    print(f'  frame {i:4n}/{len(ts)}', flush=True)

# Compile and save animation
print('Processing frames...', flush=True)
title = f'{problem_id}'
anim = animation.FuncAnimation(fig, animate, frames=len(ts), repeat=False)
metadata = dict(title=(title+' mean radiation energy density'),
                artist='Stanley A. Baronett')
writer = animation.FFMpegWriter(fps=60, metadata=metadata, bitrate=-1)
anim.save(title+'_Er.mp4', writer=writer,
          savefig_kwargs={'bbox_inches': 'tight',
                          'pad_inches' : '0.01'})
print('Done.\nVideo saved.', flush=True)

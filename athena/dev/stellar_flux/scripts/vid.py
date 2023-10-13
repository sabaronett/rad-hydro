#!/usr/bin/env python3
#==============================================================================
# vid.py
#
# Create an MPEG-4 animation of the radial flux field.
#
# Author: Stanley A. Baronett
# Created: 2022-09-20
# Updated: 2023-10-13
#==============================================================================
import sys
sys.path.insert(0, '/mnt/home/sbaronett/github/PrincetonUniversity/athena/vis/python')
import athena_read
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib import animation
import numpy as np
from pathlib import Path

def animate(i):
    """
    Update frame.

    Parameters
    ----------
        i : int
            Frame number.
    """
    ax.set_title(f'$t={times[i]:.0f}t_0$')
    clip = np.clip(Fr1s[i], vmin, vmax)
    mesh.set_array(clip)
    mesh.set_clim(vmin, vmax)
    print(f'  frame {i:4n}/{len(times)}', flush=True)

# conversions
length_unit  = 5.98e14             # L_0 [cm]
au = 1.495978707e13                # [cm]
Lau = length_unit/au
rad2deg = 180/np.pi

# ID run, set constants
problem_id = 'stellar_flux'
npsi = 4
nzeta = 5
vmin, vmax = 7e-7, 0.7
fig, ax = plt.subplots(figsize=(8, 4.5), dpi=150)

# Read and plot
path = f'{problem_id}/npsi{npsi}/nzeta{nzeta}'
athinput = athena_read.athinput(f'athinput.{problem_id}')
outputs = sorted(list(Path('athdf').glob(problem_id + '.out1.*.athdf')))
athdf = athena_read.athdf(outputs[0])
xv, yv = athdf['x1v'], athdf['x2v']
times, Fr1s = [], []

for output in outputs:
    athdf = athena_read.athdf(output)
    times.append(athdf['Time'])
    Fr1s.append(athdf['Fr1'][0])

# Initialize first frame
clip = np.clip(Fr1s[0], vmin, vmax)
mesh = ax.pcolormesh(xv*Lau, yv*rad2deg, clip, norm=colors.LogNorm())
cb = plt.colorbar(mesh, label='$F_{\mathrm{r},x}/(ca_\mathrm{r}T_0^4)$')

# Format (sub)plots
ax.minorticks_on()
ax.set_title(f'$t={times[0]:.0f}t_0$')#, fontsize='medium')
ax.set(xlabel='$r$/au', ylabel=r'$\theta^\circ$')
ax.tick_params(axis='both', which='both', top=True, right=True)
plt.gca().invert_yaxis()

# Compile and save animation
print('Processing frames...', flush=True)
title = f'{problem_id}'
anim = animation.FuncAnimation(fig, animate, frames=len(times), repeat=False)
metadata = dict(title=(title+' radial flux field'),
                artist='Stanley A. Baronett')
writer = animation.FFMpegWriter(fps=60, metadata=metadata, bitrate=-1)
anim.save(title+'_Fr1.mp4', writer=writer)
print('Done.\nVideo saved.', flush=True)

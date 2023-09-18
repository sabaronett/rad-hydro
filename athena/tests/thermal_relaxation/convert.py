import sys
sys.path.append('\\wsl.localhost\Ubuntu\home\stanley\bitbucket\ccyang\athena-dust\vis\python')
import athena_read
from pathlib import Path
import pyopenvdb as vdb

path = 'C:\Users\stanl\Downloads\sup\blender\athdf'
g_outputs = sorted(list(Path(path).glob('SI.out1.*.athdf')))
p_outputs = sorted(list(Path(path).glob('SI.out2.*.athdf')))

g_data = athena_read.athdf(g_outputs[-1])
rho = g_data['rho']
g_grid = vdb.FloatGrid()
g_grid.copyFromArray(rho)
g_grid.activeVoxelCount() == rho.size
g_grid.evalActiveVoxelBoundingBox()
g_grid.name='rho'
g_grid.gridClass = vdb.GridClass.FOG_VOLUME

p_data = athena_read.athdf(p_outputs[-1])
rhop = p_data['rhop']
p_grid = vdb.FloatGrid()
p_grid.copyFromArray(rhop)
p_grid.activeVoxelCount() == rhop.size
p_grid.evalActiveVoxelBoundingBox()
p_grid.name='rhop'
p_grid.gridClass = vdb.GridClass.FOG_VOLUME

vdb.write(f'C:\Users\stanl\Downloads\sup\blender\vdb\SI.out.{3:i}', grids=[g_grid,
            p_grid])


for i, g_output in enumerate(g_outputs):
    g_data = athena_read.athdf(g_output)
    rho = g_data['rho']
    g_grid = vdb.FloatGrid()
    g_grid.copyFromArray(rho)
    g_grid.activeVoxelCount() == rho.size
    g_grid.evalActiveVoxelBoundingBox()
    g_grid.name='rho'
    g_grid.gridClass = vdb.GridClass.FOG_VOLUME

    p_data = athena_read.athdf(p_outputs[i])
    rhop = p_data['rhop']
    p_grid = vdb.FloatGrid()
    p_grid.copyFromArray(rhop)
    p_grid.activeVoxelCount() == rhop.size
    p_grid.evalActiveVoxelBoundingBox()
    p_grid.name='rhop'
    p_grid.gridClass = vdb.GridClass.FOG_VOLUME

    vdb.write(f'C:\Users\stanl\Downloads\sup\blender\vdb\SI.out.{i:03}', grids=[g_grid,
            p_grid])

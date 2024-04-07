#!/usr/bin/env python3
#===============================================================================
# problem_setup.py
#
# Creates RADMC-3D input files for the protoplanetary disk model used in our
# analytic ray-tracing solutions. This script is adapted from
# `radmc3d-2.0/examples/run_ppdisk_simple_1/problem_setup.py`.
#
# Author: Stanley A. Baronett
# Created: 2023-12-11
# Updated: 2024-04-04
#===============================================================================
import numpy as np

# BEGIN `amr_grid.inp`==========================================================
# Constants
au       = 1.495978707e13                     # Astronomical Unit       [cm]

# Grid parameters
nr       = 256
ntheta   = 512
nphi     = 1
rin      = 10*au                              # [cm]
rout     = 100*au                             # [cm]
thetamin = 0
thetamax = np.pi
phimin   = 0
phimax   = 2*np.pi

# Make the coordinates
ri       = np.linspace(rin, rout, nr+1)
thetai   = np.linspace(thetamin, thetamax, ntheta+1)
phii     = np.linspace(phimin, phimax, nphi+1)

# Write the grid file
# https://www.ita.uni-heidelberg.de/~dullemond/software/radmc-3d/manual_radmc3d/inputoutputfiles.html#input-required-amr-grid-inp
with open('amr_grid.inp', 'w+') as f:
    f.write('1\n')                            # iformat
    f.write('0\n')                            # AMR grid style (0=regular grid)
    f.write('100\n')                          # Coordinate system: spherical
    f.write('0\n')                            # gridinfo
    f.write('1 1 0\n')                        # Include r,theta coordinates
    f.write(f'{nr:d} {ntheta:d} {nphi:d}\n')  # Size of grid
    for value in ri:
        f.write(f'{value:13.6e}\n')           # X coordinates (cell walls)
    for value in thetai:
        f.write(f'{value:13.6e}\n')           # Y coordinates (cell walls)
    for value in phii:
        f.write(f'{value:13.6e}\n')           # Z coordinates (cell walls)
# END `amr_grid.inp`============================================================


# BEGIN `wavelength_micron.inp`=================================================
# Wavelength grid
lam1     = 0.1e0
lam2     = 7.0e0
lam3     = 25.e0
lam4     = 1.0e4
n12      = 20
n23      = 100
n34      = 30
lam12    = np.logspace(np.log10(lam1), np.log10(lam2), n12, endpoint=False)
lam23    = np.logspace(np.log10(lam2), np.log10(lam3), n23, endpoint=False)
lam34    = np.logspace(np.log10(lam3), np.log10(lam4), n34, endpoint=True)
lam      = np.concatenate([lam12, lam23, lam34])
nlam     = lam.size

# Write the wavelength file
# https://www.ita.uni-heidelberg.de/~dullemond/software/radmc-3d/manual_radmc3d/inputoutputfiles.html#input-required-wavelength-micron-inp
with open('wavelength_micron.inp', 'w+') as f:
    f.write(f'{nlam:d}\n')
    for value in lam:
        f.write(f'{value:13.6e}\n')
# END `wavelength_micron.inp`===================================================


# BEGIN `stars.inp`=============================================================
# Constants
ms       = 1.98892e33                         # Solar mass              [g]
ts       = 5.78388e3                          # Solar temperature       [K]
rs       = 6.9368e10                          # Solar radius            [cm]

# Star parameters
mstar    = ms                                 # unimportant as of verion 2.0
rstar    = rs                                 # for extended source treatment
tstar    = ts                                 # blackbody temperature
pstar    = np.array([0.,0.,0.])               # position

# Write the stars.inp file
# https://www.ita.uni-heidelberg.de/~dullemond/software/radmc-3d/manual_radmc3d/inputoutputfiles.html#sec-stars
with open('stars.inp', 'w+') as f:
    f.write('2\n')
    f.write(f'1 {nlam:d}\n\n')
    f.write(f'{rstar:13.6e} {mstar:13.6e} {pstar[0]:13.6e} {pstar[1]:13.6e} '\
            + f'{pstar[2]:13.6e}\n\n')
    for value in lam:
        f.write(f'{value:13.6e}\n')
    f.write(f'\n{-tstar:13.6e}\n')
# END `stars.inp`===============================================================


# BEGIN `dust_density.inp`======================================================
# Make the cell-centered grid
rc       = 0.5*(ri[0:nr] + ri[1:nr+1])
thetac   = 0.5*(thetai[0:ntheta] + thetai[1:ntheta+1])
phic     = 0.5*(phii[0:nphi] + phii[1:nphi+1])
qq       = np.meshgrid(rc, thetac, indexing='ij')
rr       = qq[0]
tt       = qq[1]
# zr       = np.pi/2.e0 - qq[1]

# Disk parameters
gm0      = 1.0                                # GM
T_unit   = 6.14e3                             # T_0 [K]
density_unit = 4.28e-14                       # \rho_0 [g/cm^3]
length_unit = 5.98e14                         # L_0 [cm]
dfloor   = 1e-12                              # minimum density [g/cm^3]
r0       = 0.425278227742474                  # disk radial normalization [L_0]
rho0     = 0.2                                # disk density normalization [\rho_0]
p0_over_r0 = 4.80e-03                         # (H/r0)^2
pslope   = -0.5                               # pressure power-law index
dslope   = -2.25                              # density power-law index

# Make the dust density model
def GetCylCoord(x1, x2, x3):
    rad = np.abs(x1*np.sin(x2))
    phi = x3
    z = x1*np.cos(x2)
    return rad, phi, z

def DenProfileCyl(rad, phi, z):
    p_over_r = PoverR(rad, phi, z)
    denmid = rho0*np.power((rad + r0)/r0, dslope)\
            /(1 + np.exp(-np.exp(np.e)*(rad - r0)/r0))
    dentem = denmid*np.exp(gm0/p_over_r*(1./np.sqrt(rad**2 + z**2) - 1./rad))
    den = dentem
    return den

def PoverR(rad, phi, z):
    poverr = p0_over_r0*np.power(rad/r0, pslope)
    return poverr

rads, phis, zs = GetCylCoord(rr/length_unit, tt, 0)
rhod     = DenProfileCyl(rads, phis, zs)
rhod[rhod < dfloor] = dfloor
rhod     = rhod*density_unit

# Write the density file
# https://www.ita.uni-heidelberg.de/~dullemond/software/radmc-3d/manual_radmc3d/inputoutputfiles.html#input-required-for-dust-transfer-dust-density-inp
with open('dust_density.inp', 'w+') as f:
    f.write('1\n')                            # Format number
    f.write(f'{nr*ntheta*nphi:d}\n')          # Nr of cells
    f.write('1\n')                            # Nr of dust species
    data = rhod.ravel(order='F')              # Fortran-style indexing
    data.tofile(f, sep='\n', format='%13.6e')
    f.write('\n')
# END `dust_density.inp`========================================================


# BEGIN `dustopac.inp`==========================================================
# Dust opacity control file
# https://www.ita.uni-heidelberg.de/~dullemond/software/radmc-3d/manual_radmc3d/inputoutputfiles.html#the-dustopac-inp-file
with open('dustopac.inp', 'w+') as f:
    f.write('2               iformat (Format number of this file)\n')
    f.write('1               nspec (Nr of dust species)\n')
    f.write('============================================================================\n')
    f.write('1               inputstyle (Way in which this dust species is read)\n')
    f.write('0               iquantum (0=Thermal grain)\n')
    # f.write('silicate        Extension of name of dustkappa_***.inp file\n')
    f.write('dsharp_abs      Extension of name of dustkappa_***.inp file\n')
    f.write('----------------------------------------------------------------------------\n')
# END `dustopac.inp`============================================================


# BEGIN `radmc3d.inp`===========================================================
# Monte Carlo parameters
nphot    = int(1e8)
countwrite = int(1e5)

# Write the radmc3d.inp control file
# https://www.ita.uni-heidelberg.de/~dullemond/software/radmc-3d/manual_radmc3d/clioptions.html#additional-arguments-general
with open('radmc3d.inp', 'w+') as f:
    f.write(f'countwrite = {countwrite}\n')   # nr photons between std. outputs
    f.write('iranfreqmode = 1\n')             # differ RNG seeds per OMP process
    f.write('istar_sphere = 0\n')             # point (0)/sphere (1) star source
    f.write('mc_weighted_photons = 0\n')      # focus photons toward grid
    f.write(f'nphot = {nphot:d}\n')           # number of photons
    f.write('scattering_mode_max = 0\n')      # no scattering (zero dust albedo)
    f.write('setthreads = 18\n')              # Intel Core i7-12700H optimized
# END `radmc3d.inp`=============================================================

#!/usr/bin/env python3
#==============================================================================
# athena_opac.py
#
# Creates multifrequency opacity tables for input into Athena++.
#
# This script takes as input (1) a RADMC-3D-formatted opacity table with the
# filename structure `dustkappa_*.inp` (see https://www.ita.uni-heidelberg.de/
# ~dullemond/software/radmc-3d/manual_radmc3d/
# inputoutputfiles.html#the-dustkappa-inp-files); (2) an Athena++-formatted
# input file with the following requred parameters:
#   <radiation>
#   n_frequency    # no. of frequency groups
#   frequency_min  # [0, \nu_l) left group boundary [k_BT_0/h]
#   frequency_max  (if n_frequency > 2) # [\nu_r, \inf) right group boundary
#
#   <problem>
#   n_temperature    # no. of temperature groups
#   temperature_min  # min mean opacity temperature [K]
#   temperature_max  # max mean opacity temperature [K]
#
# Author: Stanley A. Baronett
# Created: 2024-04-19
# Updated: 2024-04-20
#==============================================================================
import numpy as np
from pathlib import Path
from radmc3dPy import analyze
from radmc3dPy.natconst import *
from scipy import integrate
from scipy.constants import c, h, k
import sys
sys.path.insert(0,'/home/stanley/github/PrincetonUniversity/athena/vis/python')
import athena_read

def GetBnu_table(Ts, nus):
    """Computes Planck's law for a table of temperatures and frequencies

    B_\nu(\nu, T) = \frac{2h\nu^3}{c^2}\frac{1}{e^{h\nu/(k_\mathrm{B}T)} - 1}
    """
    table = np.zeros((len(Ts), len(nus)))

    for i, T in enumerate(Ts):
        for j, nu in enumerate(nus):
            exp = np.exp(h*nu/k/T)
            if exp > sys.float_info.max:
                exp = sys.float_info.max
            table[i][j] = 2*h*nu**3/c**2/(exp - 1)
    
    table = np.where(table < 5e-324, 5e-324, table)

    return table

def GetdBnu_dT_table(Ts, nus):
    """Partial derivative of Planck's law with respect to temperature

    \frac{\partial B_\nu}{\partial T} =
        \frac{2h^2\nu^4e^{h\nu/(k_\mathrm{B}T)}}{c^2k_\mathrm{B}T^2}
        \frac{1}{(e^{h\nu/(k_\mathrm{B}T)} - 1)^2}
    """
    table = np.zeros((len(Ts), len(nus)))

    for i, T in enumerate(Ts):
        for j, nu in enumerate(nus):
            exp = np.exp(h*nu/k/T)
            if exp > sys.float_info.max:
                exp = sys.float_info.max
                denom = sys.float_info.max
            else:
                denom = c**2*k*T**2*(exp - 1)**2
                if denom > sys.float_info.max:
                    denom = sys.float_info.max
            numer = 2*h**2*nu**4*exp
            table[i][j] = numer/denom

    table = np.where(table < 5e-324, 5e-324, table)

    return table

def BinarySearchIncreasing(arr, low, high, target):
    """Iterative binary search on a strictly increasing array.

    Iteratively use binary search on a strictly increasing array to find the
    index that right-brackets the target, arr[mid-1] < target < arr[mid]
    """    
    while (low <= high):
        mid = int(low + (high - low)//2)
        if ((arr[mid-1] < target) and (target < arr[mid])):
            return mid
        elif (arr[mid] < target):
            low = mid
        else:
            high = mid

    raise Exception("Array may not be strictly increasing")

def RosselandMeanOpacities(kappa_nu, dBnu_dT, nu):
    numer = integrate.simpson(dBnu_dT/kappa_nu, x=nu)
    denom = integrate.simpson(dBnu_dT, x=nu)
    kappa = denom/numer
    return kappa

def PlanckMeanOpacities(kappa_nu, Bnu, nu, temp_table):
    numer = integrate.simpson(kappa_nu*Bnu, x=nu)
    denom = integrate.simpson(Bnu, x=nu)
    kappa = numer/denom
    return kappa

# Read absorption coefficient as a function of frequency
fname = list(Path('./').glob(f'dustkappa_*.inp'))[0].parts[0]
ext = fname[10:-4]
opac = analyze.readOpac(ext=['dsharp'])
opac_freq = np.flip(1e6*c/opac.wav[0])
opac_kabs = np.flip(opac.kabs[0])
# ksca = np.flip(opac.ksca[0])

# Make tables to compute and save mean opacities
fname = list(Path('./').glob(f'athinput.*'))[0].parts[0]
athinput = athena_read.athinput(fname)
n_frequency = athinput['radiation']['n_frequency']
T_unit = athinput['radiation']['T_unit']                          # [K]
frequency_min = athinput['radiation']['frequency_min']*k*T_unit/h # [Hz]
density_unit = athinput['radiation']['density_unit']              # [g/cm^3]
length_unit = athinput['radiation']['length_unit']                # [cm]
n_temperature = athinput['problem']['n_temperature']
temperature_min = athinput['problem']['temperature_min']          # [K]
temperature_max = athinput['problem']['temperature_max']          # [K]
ff = np.asarray(frequency_min)       # frequency group f interfaces [Hz]
temp_table = np.logspace(np.log10(temperature_min), np.log10(temperature_max),
                         n_temperature)
Bnu_table = GetBnu_table(temp_table, opac_freq)
dBnu_dT_table = GetdBnu_dT_table(temp_table, opac_freq)
kappa_af_table = np.zeros((n_temperature, n_frequency))
kappa_pf_table = np.zeros((n_temperature, n_frequency))

if n_frequency > 2:
    frequency_max = athinput['radiation']['frequency_max']*k*T_unit/h # [Hz]
    ff = np.logspace(np.log10(frequency_min), np.log10(frequency_max),
                             n_frequency-1)

ff = np.insert(ff, 0, 0)
ff = np.append(ff, float('inf'))
i_nu0 = 0
i_nu1 = BinarySearchIncreasing(opac_freq, 0, len(opac_freq)-1, ff[1])

for i in range(n_frequency):
    kappa_af_table[:, i] = RosselandMeanOpacities(opac_kabs[i_nu0:i_nu1],
                                                  dBnu_dT_table[:,i_nu0:i_nu1],
                                                  opac_freq[i_nu0:i_nu1])
    kappa_pf_table[:, i] = PlanckMeanOpacities(opac_kabs[i_nu0:i_nu1],
                                               Bnu_table[:,i_nu0:i_nu1],
                                               opac_freq[i_nu0:i_nu1],
                                               temp_table)
    i_nu0 = i_nu1
    if i < (n_frequency - 2):  # intermediate frequency group
        i_nu1 = BinarySearchIncreasing(opac_freq, 0, len(opac_freq)-1, ff[i+2])
    else:                      # (next-to-) last frequency group
        i_nu1 = len(opac_freq)-1

# Convert units from cgs to code
temp_table /= T_unit                                             # [T_0]
kappa_af_table *= density_unit*length_unit                       # [\rho_0*L_0]
kappa_pf_table *= density_unit*length_unit                       # [\rho_0*L_0]

# Save tables to text files for Athena++ input
np.savetxt('temp_table.txt', temp_table)
# np.savetxt('kappa_sf_table.txt', kappa_sf_table)
np.savetxt('kappa_rf_table.txt', kappa_af_table)
np.savetxt('kappa_pf_table.txt', kappa_pf_table)

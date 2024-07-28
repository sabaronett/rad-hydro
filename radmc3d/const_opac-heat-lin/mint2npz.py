#!/usr/bin/env python3
#===============================================================================
# mint2npz.py
#
# Uses the radmc3dPy module to read `./mean_intensity.out`, integrate it over
# all frequencies (in decreasing order), and save it to tot_mint.npz.
#
# Author: Stanley A. Baronett
# Created: 2024-07-25
# Updated: 2024-07-27
#===============================================================================
import numpy as np
from radmc3dPy import *
from scipy import integrate

data = analyze.readData(mint=True)
print('Integrating over all frequencies', flush=True)
tot_mint = np.abs(integrate.trapezoid(data.meanint, x=data.freq, axis=-1))
print('Saving to tot_mint.npz', flush=True)
np.savez_compressed('tot_mint', tot_mint=tot_mint)

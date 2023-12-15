#!/usr/bin/env python3
#===============================================================================
# midplane_pressure.py
#
# Finds and saves the time-varying pressure along the disk midplane.
#
# Author: Stanley A. Baronett
# Created: 2023-12-15
# Updated: 2023-12-15
#===============================================================================
import sys
sys.path.insert(0, '/mnt/home/sbaronett/github/sabaronett/athena/vis/python')
import athena_read
import numpy as np
from pathlib import Path

# Collect Athena++ inputs, outputs, and sim constants
athinput = athena_read.athinput('athinput.si')
j_mid = athinput['mesh']['nx2']//2
outputs = sorted(list(Path('athdf').glob(athinput['job']['problem_id']+
                                         '.out1.*.athdf')))
press, t = [], []

print(f'Compiling data...', flush=True)

for i, output in enumerate(outputs):
    athdf = athena_read.athdf(output)
    press.append(athdf['press'][0, j_mid])
    t.append(athdf['Time'])
    print('\t{:.2%}'.format(i/len(outputs)), flush=True)

print(f'\tDone.\nSaving results...', flush=True)
np.savez_compressed('midplane_pressure', press=press, t=t)

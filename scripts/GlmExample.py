import numpy as np
import matplotlib.pyplot as plt

from map_tools import *
from map_cycle import *
from map_run import *


### edit plot parameters here:
plt.rcParams["font.size"]=10
plt.rcParams["lines.markersize"]=3
plt.rcParams["figure.figsize"]=[12, 6]
plt.rcParams["figure.dpi"]=200

### create a Run object with a run number, a range of cycles, the current data format, ask it to calibrate, and ask it to cut off data taken when not moving:
B0_map = Run(1926, np.arange(0, 150), dtformat='v4_PhiZswap', calib=True, trim=True)

### /!\ note on maps directory: Run looks for maps in "../maps_bin/"
### /!\ note on calibration: meant to be used when the run contains a calibration ring (90° or 180° rotation).
### If calib==True and no calibration ring is found, the z probe will calibrate to an arbitrary value

### Cycle object can be extracted from a run by its [rho, z] index:
some_ring = B0_map.cycles['rings'][40, -41]

### plot example for some ring:
fig, ax = some_ring.simplePlot('phi', 'Brho')
plt.show()

### extract Glms coefs and their uncertainties up to order l=3 from the z probe:
lmax=3
G, G_err = B0_map.getGlm(lmax, source='z', norm=True)
### G[l][lmax + 1 + m] is G_lm
print("G_10 = {:.2e} +/- {:.2e}".format(G[1][lmax+1+0], G_err[1][lmax+1+0]))

### plot Glm bars of 0.1 width for l=1 and all m
B0_map.plotGm(l=1, w=0.1)
plt.show()

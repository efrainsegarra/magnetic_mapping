import numpy as np
import matplotlib.pyplot as plt

from map_tools import *
from map_cycle import *


### edit plot parameters here:
plt.rcParams["font.size"]=10
plt.rcParams["lines.markersize"]=3
plt.rcParams["figure.figsize"]=[12, 6]
plt.rcParams["figure.dpi"]=100

### create a Cycle object with a run number, a cycle number, and the current data format.
some_cycle = Cycle(run=1926, cyc=20, dtformat='v4_PhiZswap')
### /!\ note on data format: use 'v4' for the original fluxgate orientation (up to runs 18**)
### and 'v4_PhiZswap' for the 90° flip of the original that is currently being used (runs 19**+).
### /!\ note on maps directory: Cycle looks for maps in "../maps_bin/"

### fluxgate data is stored in numy arrays some_cycle.x with x = t, rho, phi, z, Brho, Bphi, Bz.
### unit is (s) for t, (cm) for rho.., (pT) for Brho....
print("This is a {} scan at rho={:.0f}cm, z={:.0f}cm".format(some_cycle.scantype, some_cycle.rho[0], some_cycle.z[0]))

### use the simplePlot method for a straightforward plot:
fig, ax = some_cycle.simplePlot(xlabel='phi', ylabel='Bz', fontsize=15)
### possible xlabels and ylabels: 't', 'rho', 'phi', 'z', 'Brho', 'Bphi', 'Bz'.
plt.show()

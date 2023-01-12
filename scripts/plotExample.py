import numpy as np
import matplotlib.pyplot as plt

from map_tools import *
from map_cycle import *

# 1963 - standard
#   z scan cyc 31
# 1964
#   cyc 31 = rho 50, z -30
#   z scan cyc 28

### edit plot parameters here:
plt.rcParams["font.size"]=10
plt.rcParams["lines.markersize"]=3
plt.rcParams["figure.figsize"]=[12, 6]
plt.rcParams["figure.dpi"]=100

### create a Cycle object with a run number, a cycle number, and the current data format.
#x_cycle = Cycle(run=1963, cyc=22, dtformat='v4_PhiZswap')
cycle = Cycle(run=1977, cyc=11, dtformat='v4_PhiZswap')
#y_cycle = Cycle(run=1963, cyc=23, dtformat='v4_PhiZswap')
#z_cycle = Cycle(run=1963, cyc=24, dtformat='v4_PhiZswap')

# new degaussing response is cyc 22
# old degaussing response is cyc 22-24 (that's why after degaussing the mapping is old = new + 3

# new degaussing cyc 31 is z=-30, r=50 ring (after degauss)
# old degaussing should be 34

### /!\ note on data format: use 'v4' for the original fluxgate orientation (up to runs 18**)
### and 'v4_PhiZswap' for the 90° flip of the original that is currently being used (runs 19**+).
### /!\ note on maps directory: Cycle looks for maps in "../maps_bin/"

### fluxgate data is stored in numy arrays before_cycle.x with x = t, rho, phi, z, Brho, Bphi, Bz.
### unit is (s) for t, (cm) for rho.., (pT) for Brho....
print("This is a {} scan at rho={:.0f}cm, z={:.0f}cm".format(cycle.scantype, cycle.rho[0], cycle.z[0]))

### use the simplePlot method for a straightforward plot:
#fig, ax, t_x, Bphi_x  = x_cycle.simplePlot(xlabel='t', ylabel='Bz', fontsize=15)
#fig, ax, t_y, Bphi_y  = y_cycle.simplePlot(xlabel='t', ylabel='Bz', fontsize=15)
#fig, ax, t_z, Bphi_z  = z_cycle.simplePlot(xlabel='t', ylabel='Bz', fontsize=15)
#t = np.concatenate((t_x,t_y,t_z),axis=None)
#Bphi = np.concatenate((Bphi_x,Bphi_y,Bphi_z),axis=None)
#plt.figure(0)
#plt.plot(t_x,Bphi_x,alpha=0.5,color='blue',label='L6-X')
#plt.plot(t_y,Bphi_y,alpha=0.5,color='red',label='L6-Y')
#plt.plot(t_z,Bphi_z,alpha=0.5,color='green',label='L6-Z')
#plt.legend(numpoints=1,loc='best',fontsize=13)
#ax = plt.gca()
#ax.set_xlabel(r"{}".format((r'$t$ (s)')), size=15)
#ax.set_ylabel(r"{}".format((r'$B_z$ (pT)')), size=15)
#plt.savefig('1963_degaussresponse_Bz.png',bbox_inches='tight')
#plt.xticks(fontsize=14)
#plt.yticks(fontsize=14)
#plt.savefig('1964_degaussresponse_Bphi.png',bbox_inches='tight')
#fig, ax, t, Brho  = cycle.simplePlot(xlabel='phi', ylabel='Bz', fontsize=15)
#plt.savefig('1964_degaussresponse_Brho.png',bbox_inches='tight')
fig, ax, t, Bz    = cycle.simplePlot(xlabel='z', ylabel='Bz', fontsize=15)
#plt.savefig('1964_degaussresponse_Bz.png',bbox_inches='tight')
#sBphi = np.abs(np.max(Bphi) - np.min(Bphi))
#sBrho = np.abs(np.max(Brho) - np.min(Brho))
#sBz   = np.abs(np.max(Bz  ) - np.min(Bz  ))
#print(np.round(cycle.rho[0]), np.round(cycle.z[0]), sBphi, sBrho, sBz )


### possible xlabels and ylabels: 't', 'rho', 'phi', 'z', 'Brho', 'Bphi', 'Bz'.

plt.show()

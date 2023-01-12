
import numpy as np
import matplotlib.pyplot as plt



from map_cycle import *


def grabCycleData(Run,Cyc):
    cycle = Cycle(run=Run, cyc=Cyc, dtformat='v4_PhiZswap')
    return  cycle.scantype,\
            cycle.rho[0],\
            cycle.phi[0],\
            cycle.z[0],\
            getattr(cycle,'t'),\
            getattr(cycle,'rho'),\
            getattr(cycle,'phi'),\
            getattr(cycle,'z'),\
            getattr(cycle,'Brho'),\
            getattr(cycle,'Bphi'),\
            getattr(cycle,'Bz')

def transXYZ(rho,phi,z):
    x = rho*np.cos(phi)
    y = rho*np.sin(phi)
    return x,y,z

ax = plt.axes(projection='3d')

scantype, initrho, initphi, initz, t,rho,phi,z,Brho,Bphi,Bz = grabCycleData(2037,11)
x1,y1,z1 = transXYZ(rho,phi,z)
scantype, initrho, initphi, initz, t,rho,phi,z,Brho,Bphi,Bz = grabCycleData(2037,13)
x2,y2,z2 = transXYZ(rho,phi,z)
scantype, initrho, initphi, initz, t,rho,phi,z,Brho,Bphi,Bz = grabCycleData(2037,15)
x3,y3,z3 = transXYZ(rho,phi,z)
scantype, initrho, initphi, initz, t,rho,phi,z,Brho,Bphi,Bz = grabCycleData(2037,17)
x4,y4,z4 = transXYZ(rho,phi,z)
scantype, initrho, initphi, initz, t,rho,phi,z,Brho,Bphi,Bz = grabCycleData(2037,19)
x5,y5,z5 = transXYZ(rho,phi,z)
scantype, initrho, initphi, initz, t,rho,phi,z,Brho,Bphi,Bz = grabCycleData(2037,4)
x6,y6,z6 = transXYZ(rho,phi,z)

ax.scatter(x1,y1,z1)
ax.scatter(x2,y2,z2)
ax.scatter(x3,y3,z3)
ax.scatter(x4,y4,z4)
ax.scatter(x5,y5,z5)
ax.scatter(x6,y6,z6)

ax.view_init(30, 50)
ax.xaxis.set_tick_params(labelsize=13)
ax.yaxis.set_tick_params(labelsize=13)
ax.zaxis.set_tick_params(labelsize=13)
ax.set_xlabel('x (cm)',fontsize=15)
ax.set_ylabel('y (cm)',fontsize=15)
ax.set_zlabel('z (cm)',fontsize=15)

plt.tight_layout()
plt.savefig("mapper_movement.pdf",bbox_inches='tight',pad_inches=0.5)
plt.show()

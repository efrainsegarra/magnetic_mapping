import numpy as np
import matplotlib.pyplot as plt

#from map_tools import *
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

def getSpread(Brho,Bphi,Bz):
    sBrho = np.abs(np.max(Brho) - np.min(Brho))
    sBphi = np.abs(np.max(Bphi) - np.min(Bphi))
    sBz   = np.abs(np.max(Bz  ) - np.min(Bz  ))
    return sBrho, sBphi, sBz

def makePlot(ylab,Data,ifig):
    if ylab == r'$B_\rho$':
        ind = 1
    elif ylab == r'$B_\phi$':
        ind = 2
    elif ylab == r'$B_z$':
        ind = 3
    else:
        print("Incorrect index for plotting data. Please check.\n\texiting...")
        exit()

    plt.figure(ifig)

    for key in Data:
        plt.plot(Data[key][0],Data[key][ind],linestyle='-',color=Data[key][4], alpha=0.7  ,label='Amp = {}'.format( Data[key][5]), zorder=-key)

    plt.grid(True)
    plt.xticks(fontsize=13)
    plt.yticks(fontsize=13)
    plt.xlabel(r'$t$ (s)',fontsize=15)
    plt.ylabel(ylab+r' (nT) ',fontsize=15)
    plt.legend(numpoints=1,loc='best',fontsize=13)
    plt.xlim([0,600])
    plt.tight_layout()




CycleList = [26]
RunList = [1978,1979,1980,1981,1982]
COLS = ['blue','red','green','grey','purple']
Amps = ['0.3','0.36','0.42','0.48','0.54']
Data = {}
for r in range(len(RunList)):
    Run = RunList[r]
    for Cyc in CycleList:
        scantype, initrho, initphi, initz, t,rho,phi,z,Brho,Bphi,Bz = grabCycleData(Run,Cyc)
        if scantype != 'static':
            print("Cycle list is wrong, please check.\n\texiting...")
            exit()
        
        if Run not in Data.keys():
            Data[Run] = [ [], [], [], [], [], [] ]
        Data[Run][0] = t    
        Data[Run][1] = Brho/1e3
        Data[Run][2] = Bphi/1e3
        Data[Run][3] = Bz/1e3
        Data[Run][4] = COLS[ r ]
        Data[Run][5] = Amps[ r ]


makePlot(r'$B_\rho$',Data,1)
makePlot(r'$B_\phi$',Data,2)
makePlot(r'$B_z$'   ,Data,3)

plt.figure(1)
plt.savefig('response_degauss_Brho.pdf',bbox_inches='tight')
plt.figure(2)
plt.savefig('response_degauss_Bphi.pdf',bbox_inches='tight')
plt.figure(3)
plt.savefig('response_degauss_Bz.pdf',bbox_inches='tight')

plt.show()

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

def avgData(Data):
    newData = {}
    scan = 4
    avg_ind = 20
    for Run in Data.keys():
        newData[Run] = [ [] for i in range(len(Data[Run])) ]
        for dat_ind in range(len(Data[Run])):
            i = 0
            while i < len(Data[Run][dat_ind][scan])-avg_ind:
                avg = np.average(Data[Run][dat_ind][scan][i:i+avg_ind])
                newData[Run][dat_ind].append(avg)
                i+=avg_ind
        newData[Run][4] = Data[Run][4]
        newData[Run][5] = Data[Run][5]
    return newData

def avgArr(arr):
    newarr = []
    avg_ind = 20
    i = 0
    while i < len(arr)-avg_ind:
        avg = np.average(arr[i:i+avg_ind])
        newarr.append(avg)
        i+=avg_ind
    return np.asarray(newarr)


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

    scan = 4
    for Run in Data.keys():
        style = '--'
        if Run == 1977:
            style = ':'
        plt.plot(Data[Run][0],Data[Run][ind],marker='o',linestyle=style,color=Data[Run][4] ,label=Data[Run][5] )

    #plt.plot(Data[1977][0],Data[1977][ind],marker='o',linestyle=':',color='black'  ,label='Standard')
    #plt.plot(Data[1978][0],Data[1978][ind],marker='o',linestyle='--',color='red'  ,label='XsYcZs, A=0.30, DT=250')
    #plt.plot(Data[1979][0],Data[1979][ind],marker='o',linestyle='--',color='green'  ,label='XsYcZs, A=0.36, DT=300')
    #plt.plot(Data[1980][0],Data[1980][ind],marker='o',linestyle='--',color='orange'  ,label='XsYcZs, A=0.42, DT=350')
    #plt.plot(Data[1981][0],Data[1981][ind],marker='o',linestyle='--',color='purple'  ,label='XsYcZs, A=0.48, DT=400')
    #plt.plot(Data[1982][0],Data[1982][ind],marker='o',linestyle='--',color='magenta'  ,label='XsYcZs, A=0.54, DT=450')
    #plt.plot(Data[1964][0],Data[1964][ind],marker='o',linestyle='--',color='green' ,label='XsYsZs')
    #plt.plot(Data[1965][0],Data[1965][ind],marker='o',linestyle='--',color='red'   ,label='XsYsZc')
    #plt.plot(Data[1966][0],Data[1966][ind],marker='o',linestyle='--',color='blue',label='XsYcZs, A=0.6, DT=500')
    #plt.plot(Data[1967][0],Data[1967][ind],marker='o',linestyle='--',color='purple',label='XsYcZc')

    plt.grid(True)
    plt.xticks(fontsize=13)
    plt.yticks(fontsize=13)
    plt.xlabel(r'$z$ (cm)',fontsize=15)
    plt.ylabel(ylab+' (pT) with $r=0,\phi=0$',fontsize=15)
    plt.legend(numpoints=1,loc='best',fontsize=12)
    plt.tight_layout()



# NEW DEGAUSSING (1964-1967) 
#   cyc 24 = zscan 0
#   cyc 25 = zscan 1
#   cyc 26 = zscan 2
#   cyc 27 = zscan 3
#   cyc 28 = zscan 4
# STANDARD DEGAUSSING REDO (1972-1977)
#   cyc 4 = zscan 0
#   cyc 5 = zscan 1
#   cyc 6 = zscan 2
#   cyc 7 = zscan 3
#   cyc 8 = zscan 4
# DEGAUSSING TESTS (1978-1982) with varying amp/dt
#   cyc 30 = zscan 0
#   cyc 31 = zscan 1
#   cyc 32 = zscan 2
#   cyc 33 = zscan 3
#   cyc 34 = zscan 4

def getOffset(Run):
    CalibCycle = -1
    if Run < 1977:
        CalibCycle = 40
    if Run == 1977:
        CalibCycle = 2
    if Run >= 1978:
        CalibCycle = 28
    scantype, initrho, initphi, initz, t,rho,phi,z,Brho,Bphi,Bz = grabCycleData(Run,CalibCycle)
    if scantype != 'phi':
        print("Cycle list is wrong, please check.\n\texiting...")
        exit()

    sortBy = np.argsort(phi)
    phi = phi[sortBy]
    Brho = Brho[sortBy]
    Bphi = Bphi[sortBy]

    phi = avgArr(phi)
    Brho = avgArr(Brho)
    Bphi = avgArr(Bphi)

    ind0 = np.absolute(phi-0).argmin()
    ind180 = np.absolute(phi-180).argmin()
    rhoOff = (Brho[ind0] + Brho[ind180])/2.
    phiOff = (Bphi[ind0] + Bphi[ind180])/2.

    return rhoOff,phiOff


CycleList = [24,25,26,27,28]
RunList = [1977,1966]
COLS = ['black','blue']
LABS = ['Standard','XsYcZs']

Data = {}
for iRun in range(len(RunList)):
    Run = RunList[iRun]
    for Cyc in CycleList:
        if Run > 1967 and Run < 1978:
            Cyc = Cyc - 20
        elif Run >= 1978:
            Cyc = Cyc + 6
        rhoOff, phiOff = getOffset(Run)
        scantype, initrho, initphi, initz, t,rho,phi,z,Brho,Bphi,Bz = grabCycleData(Run,Cyc)
        sortBy = np.argsort(z)
        z = z[sortBy]

        Brho = Brho[sortBy] - rhoOff
        Bphi = Bphi[sortBy] - phiOff
        Bz   = Bz[sortBy] - np.average(Bz)
        if scantype != 'z':
            print("Cycle list is wrong, please check.\n\texiting...")
            exit()
        
        if Run not in Data.keys():
            Data[Run] = [ [], [], [], [], '', '' ]
        Data[Run][0].append( z    )
        Data[Run][1].append( Brho )
        Data[Run][2].append( Bphi )
        Data[Run][3].append( Bz   )
        Data[Run][4] =  COLS[iRun] 
        Data[Run][5] =  LABS[iRun] 

AvgData = avgData(Data)


makePlot(r'$B_\rho$',AvgData,1)
makePlot(r'$B_\phi$',AvgData,2)
makePlot(r'$B_z$'   ,AvgData,3)

plt.figure(1)
plt.savefig('zscan_Brho.png',bbox_inches='tight')
plt.figure(2)
plt.savefig('zscan_Bphi.png',bbox_inches='tight')
plt.figure(3)
plt.savefig('zscan_Bz.png',bbox_inches='tight')

plt.show()

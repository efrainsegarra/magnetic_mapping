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

def avgArr(arr):
    newarr = []
    avg_ind = 20
    i = 0
    while i < len(arr)-avg_ind:
        avg = np.average(arr[i:i+avg_ind])
        newarr.append(avg)
        i+=avg_ind
    return np.asarray(newarr)

def getOffset(Run):
    CalibCycle = []
    AvgRho = 0
    AvgPhi = 0

    if Run > 2000 and Run < 2011 or Run > 2016 and Run < 2020:
        CalibCycle.append(135)
        CalibCycle.append(154)
    if Run == 1977:
        CalibCycle.append(2)
        CalibCycle.append(21)
    if Run > 2011 and Run < 2016:
        CalibCycle.append(3)
        CalibCycle.append(22)
    if Run == 2020 or Run == 2021:
        CalibCycle.append(51)
        CalibCycle.append(70)
    if Run == 2024:
        CalibCycle.append(7)
        CalibCycle.append(26)
    if Run >= 2028 and Run <= 2037:
        CalibCycle.append(2)
        CalibCycle.append(21)
    if Run >= 1964 and Run <= 1967:
        CalibCycle.append(23)
        CalibCycle.append(40)


    for Cyc in CalibCycle:
        scantype, initrho, initphi, initz, t,rho,phi,z,Brho,Bphi,Bz = grabCycleData(Run,Cyc)
        if scantype != 'phi':
            print("Calibration list is wrong, please check.\n\texiting...")
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

        plt.figure(0)
        plt.plot(phi,Brho-rhoOff,label=r'$B_\rho$ with Cyc '+str(Cyc),color='blue')
        plt.plot(phi,Bphi-phiOff,label=r'$B_\phi$ with Cyc '+str(Cyc),color='red')

        AvgRho += rhoOff
        AvgPhi += phiOff
    plt.figure(0)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.xlim([0,360])
    plt.xlabel(r'$\phi$ (deg)',fontsize=15)
    plt.ylabel(r'$B_i$ (pT)',fontsize=15)
    plt.title("Run "+str(Run)+r" $\phi,\rho$ Calibration at $z=0$",fontsize=17)
    plt.legend(numpoints=1,loc='best')
    plt.axhline(y=0,linestyle='--',color='black',linewidth=2)
    plt.tight_layout()
    plt.savefig("Run_"+str(Run)+"_PhiRhoCalib.pdf",bbox_inches='tight')

    plt.gcf()
    plt.close()
    return AvgRho/len(CalibCycle),AvgPhi/len(CalibCycle)

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
    for Run in Data.keys():
        style = '--'
        if Run == 1977:
            style = ':'
        plt.plot(Data[Run][0],Data[Run][ind],marker='o',linestyle=style,color=Data[Run][4] ,label="Run "+str(Run)+" "+Data[Run][5] )

    plt.grid(True)
    plt.xticks(fontsize=13)
    plt.yticks(fontsize=13)
    plt.xlabel(r'$z$ (cm)',fontsize=15)
    plt.ylabel(r'Spread of '+ylab+' (pT)',fontsize=15)
    plt.legend(numpoints=1,loc='best',fontsize=12)
    plt.xlim([-32,32])
    plt.ylim([0,300])
    plt.tight_layout()



RunList = [1977,2034]
COLS = ["blue","green","orange","green","purple","blue"]
LABS = ["Serial","Optimized","3.0hr","4.5hr"]
Data = {}



for iRun in range(len(RunList)):
    Run = RunList[iRun]
    phiCyc = -1; stopCyc = -1
    if Run > 2000 and Run < 2011 or Run > 2016 and Run < 2020:
        phiCyc = 144
    if Run == 1977:
        phiCyc = 11
    if Run > 2011 and Run < 2016:
        phiCyc = 12
    if Run == 2020 or Run == 2021:
        phiCyc = 60
    if Run == 2024:
        phiCyc = 16
    if Run >= 2028 and Run <= 2037:
        phiCyc = 11
    if Run >= 1964 and Run <= 1967:
        phiCyc = 31

    stopCyc = phiCyc + 8

    # Get the phi offset for this Run and save the figure
    rhoOff,phiOff = getOffset(Run)

    while phiCyc <= stopCyc:
        scantype, initrho, initphi, initz, t,rho,phi,z,Brho,Bphi,Bz = grabCycleData(Run,phiCyc)
        print(Run,phiCyc,scantype)
        if scantype != 'phi':
            print("Cycle list is wrong, please check.\n\texiting...")
            exit()
        
        # Subtrack off offsets
        Brho -= rhoOff
        Bphi -= phiOff
        Bz -= np.average(Bz)

        # Plot the individual rings to investigate later
        plt.plot(phi,Brho,label=r'$B_\rho$')
        plt.plot(phi,Bphi,label=r'$B_\phi$')
        plt.plot(phi,Bz  ,label=r'$B_z$')
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.xlabel(r'$\phi$ (deg)',fontsize=15)
        plt.ylabel(r'$B_i$ (pT)',fontsize=15)
        plt.title("Run "+str(Run)+r" $\phi$-ring at $z=$"+str(int(np.round(initz)))+"$,r=$"+str(int(np.round(initrho))),fontsize=17)
        plt.axhline(y=0,linestyle='--',color='black',linewidth=2)
        plt.legend(numpoints=1,loc='best',fontsize=12)
        plt.ylim([-200,200])
        plt.xlim([0,360])
        plt.tight_layout()
        plt.savefig("Run_"+str(Run)+r"_PhiRing_z_"+str(int(np.round(initz)))+"_r_"+str(int(np.round(initrho)))+".pdf",bbox_inches='tight')
        plt.gcf()
        plt.close()
    
        # Calculate the spread in B_i
        sBrho,sBphi,sBz = getSpread(Brho,Bphi,Bz)
        
        # Save this data to now plot global picture
        if Run not in Data.keys():
            Data[Run] = [ [], [], [], [] , "", "" ]
        Data[Run][0].append( int(np.round(initz)) )
        Data[Run][1].append( sBrho )
        Data[Run][2].append( sBphi )
        Data[Run][3].append( sBz   )
        Data[Run][4] =  COLS[iRun] 
        Data[Run][5] =  LABS[iRun] 

        # Continue to the next phi cycle
        phiCyc += 2

makePlot(r'$B_\rho$',Data,1)
makePlot(r'$B_\phi$',Data,2)
makePlot(r'$B_z$'   ,Data,3)

plt.figure(1)
plt.savefig('spread_Brho.pdf',bbox_inches='tight')
plt.figure(2)
plt.savefig('spread_Bphi.pdf',bbox_inches='tight')
plt.figure(3)
plt.savefig('spread_Bz.pdf',bbox_inches='tight')

plt.show()

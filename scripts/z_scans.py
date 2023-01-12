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

def avgData(Data):

    # First we can average the individual z-scans over z:
    newData = {}
    avg_ind = 20
    for Run in Data.keys():
        newData[Run] = [ [] for dat_ind in range(len(Data[Run])) ]
        for dat_ind in range(len(Data[Run])):
            newData[Run][dat_ind] = [ [] for scan in range(len(Data[Run][dat_ind])) ]
            for scan in range(len(Data[Run][dat_ind])):
                i = 0
                while i < len(Data[Run][dat_ind][scan])-avg_ind:
                    avg = np.average(Data[Run][dat_ind][scan][i:i+avg_ind])
                    newData[Run][dat_ind][scan].append(avg)
                    i+=avg_ind
        newData[Run][4] = Data[Run][4]
        newData[Run][5] = Data[Run][5]

    # Now we want to loop over the scans, and grab the avg + std-dev to make ONE plot
    finalData = {}
    CURRLEN = -1
    for Run in newData.keys():
        finalData[Run] = [ [] for i in range( (len(newData[Run])-2)*2 + 2 ) ]
        PREVLEN = len(newData[Run][0][0])
        for dat_ind in range(len(newData[Run])-2):
            scan_avg = []
            scan_err = []
            for val_ind in range(len(newData[Run][0][0])):

                spread = []
                for scan in range(len(newData[Run][dat_ind])):
                    CURRLEN = len(newData[Run][dat_ind][scan])
                    if CURRLEN != PREVLEN:
                        print("Cannot do the averaging over z-scans!!!")
                        print("Exiting...")
                        exit()
                    PREVLEN = CURRLEN
                    
                    spread.append( newData[Run][dat_ind][scan][val_ind] ) 
                    
                scan_avg.append( np.average(spread) )
                scan_err.append( np.std(spread) )
                
            finalData[Run][dat_ind] = scan_avg
            finalData[Run][dat_ind + 4] = scan_err

        finalData[Run][8] = Data[Run][4]
        finalData[Run][9] = Data[Run][5]
    
    return finalData



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
    save = ''
    if ylab == r'$B_\rho$':
        ind = 1
        save = 'Brho'
    elif ylab == r'$B_\phi$':
        ind = 2
        save = 'Bphi'
    elif ylab == r'$B_z$':
        ind = 3
        save = 'Bz'
    else:
        print("Incorrect index for plotting data. Please check.\n\texiting...")
        exit()
    plt.figure(ifig)
    for Run in Data.keys():
        style = '--'
        if Run == 1977:
            style = ':'
        plt.errorbar(Data[Run][0],Data[Run][ind],yerr=Data[Run][ind+4],marker='o',linestyle=style,color=Data[Run][8] ,label="Run "+str(Run)+" "+Data[Run][9] )

    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.xlabel(r'$z$ (cm)',fontsize=15)
    plt.ylabel(ylab+r' (pT)',fontsize=15)
    #plt.title(r"Average of $z$-scan at $\phi=0,r=0$",fontsize=17)
    plt.axhline(y=0,linestyle='--',color='black',linewidth=2)
    plt.legend(numpoints=1,loc='best',fontsize=12)
    plt.ylim([-100,100])
    plt.xlim([-30,30])
    plt.tight_layout()
    plt.savefig("average_zscan_"+save+".pdf",bbox_inches='tight')



#RunList = [2002,2005,2009,2019,1977]
#COLS = ["red","green","orange","magenta","blue"]
#LABS = ["X+Y-Z- & X+Y+Z+","X+Y-Z- & X+Y+Z-","X+Y-Z- & X+Y-Z+"]
#LABS = ["X+Y+Z+ & X+Y-Z+","X+Y+Z- & X+Y-Z+","X+Y-Z+ & X+Y-Z-","X+Y-Z- & X+Y-Z+","Standard"]
#RunList = [2012,2013,2014,2015,1977]
#COLS = ["red","green","orange","magenta","blue"]
#LABS = ["X+Y-Z-","X+Y+Z-","X+Y+Z+","X+Y-Z+","Standard"]
#RunList = [2032,2033,2034,2019,1977]
#COLS = ["red","green","orange","magenta","blue"]
#LABS = ["2032","2033","2034","X+Y-Z- & X+Y-Z+ PREV","Standard"]
RunList = [1977,2034]
COLS = ["blue","green","orange","green","purple","blue"]
LABS = ["Serial","Optimized","3.0hr","4.5hr"]
Data = {}


for iRun in range(len(RunList)):
    Run = RunList[iRun]
    zCyc = -1; stopCyc = -1
    if Run > 2000 and Run < 2011 or Run > 2016 and Run < 2020:
        zCyc = 137
    if Run == 1977:
        zCyc = 4
    if Run > 2011 and Run < 2016:
        zCyc = 5
    if Run >= 2028 and Run <= 2037:
        zCyc = 4
    if Run >= 1964 and Run <= 1967:
        zCyc = 24
    stopCyc = zCyc + 5

    # Get the phi offset for this Run and save the figure
    rhoOff,phiOff = getOffset(Run)
    ctr = 0
    while zCyc < stopCyc:
        scantype, initrho, initphi, initz, t,rho,phi,z,Brho,Bphi,Bz = grabCycleData(Run,zCyc)
        if scantype != 'z':
            print("Cycle list is wrong, please check.\n\texiting...")
            exit()
        
        # Subtrack off offsets
        Brho -= rhoOff
        Bphi -= phiOff
        Bz -= np.average(Bz)

        # Plot the individual scans to investigate later
        plt.plot(z,Brho,label=r'$B_\rho$')
        plt.plot(z,Bphi,label=r'$B_\phi$')
        plt.plot(z,Bz  ,label=r'$B_z$')
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.xlabel(r'$z$ (cm)',fontsize=15)
        plt.ylabel(r'$B_i$ (pT)',fontsize=15)
        plt.title("Run "+str(Run)+r" $z$-scan at $\phi=0,r=0$",fontsize=17)
        plt.axhline(y=0,linestyle='--',color='black',linewidth=2)
        plt.legend(numpoints=1,loc='best',fontsize=12)
        plt.ylim([-100,100])
        plt.xlim([-30,30])
        plt.tight_layout()
        plt.savefig("Run_"+str(Run)+r"_zScan_"+str(ctr)+"_r_0_phi_0.pdf",bbox_inches='tight')
        plt.gcf()
        plt.close()
    
        sortBy = np.argsort(z)
        z = z[sortBy]
        Brho = Brho[sortBy]
        Bphi = Bphi[sortBy]
        Bz = Bz[sortBy]

        # Save this data to now plot global picture
        if Run not in Data.keys():
            Data[Run] = [ [], [], [], [] , "", "" ]
        Data[Run][0].append( z )
        Data[Run][1].append( Brho )
        Data[Run][2].append( Bphi )
        Data[Run][3].append( Bz   )
        Data[Run][4] =  COLS[iRun] 
        Data[Run][5] =  LABS[iRun] 

        # Continue to the next phi cycle
        zCyc += 1
        ctr += 1

AvgData = avgData(Data)

makePlot(r'$B_\rho$',AvgData,1)
makePlot(r'$B_\phi$',AvgData,2)
makePlot(r'$B_z$'   ,AvgData,3)

#plt.figure(1)
#plt.savefig('spread_Brho.pdf',bbox_inches='tight')
#plt.figure(2)
#plt.savefig('spread_Bphi.pdf',bbox_inches='tight')
#plt.figure(3)
#plt.savefig('spread_Bz.pdf',bbox_inches='tight')

plt.show()

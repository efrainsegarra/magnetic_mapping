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

def avgArr(arr):
    newarr = []
    avg_ind = 20
    i = 0
    while i < len(arr)-avg_ind:
        avg = np.average(arr[i:i+avg_ind])
        newarr.append(avg)
        i+=avg_ind
    return np.asarray(newarr)

def avgZData(Data):

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
                
            finalData[Run][dat_ind] = np.asarray(scan_avg)
            finalData[Run][dat_ind + 4] = np.asarray(scan_err)

        finalData[Run][8] = Data[Run][4]
        finalData[Run][9] = Data[Run][5]
    
    return finalData

def makeAverages(Data,Run):
    PHI = []
    BRHO = []
    BPHI = []
    BZ = []
    Z = []
    for z in range(len(Data[Run][0])):
        phi  = Data[Run][0][z]
        Brho = Data[Run][1][z]
        Bphi = Data[Run][2][z]
        Bz   = Data[Run][3][z]

        sortBy = np.argsort(phi)
        phi = phi[sortBy]
        Brho = Brho[sortBy]
        Bphi = Bphi[sortBy]
        Bz   = Bz[sortBy]

        phi  = avgArr(phi)
        Brho = avgArr(Brho)
        Bphi = avgArr(Bphi)
        Bz   = avgArr(Bz)

        PHI.append(phi)
        BRHO.append(Brho)
        BPHI.append(Bphi)
        BZ.append(Bz)
        Z.append( Data[Run][4][z] )
    return PHI, BRHO, BPHI, BZ, Z

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
    if Run >= 2038 and Run <= 2041:
        CalibCycle.append(5)
        CalibCycle.append(24)


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

def getPhiRings(RunList):
    DataPhi = {}
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
        if Run >= 2038 and Run <= 2041:
            phiCyc = 14
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
            
            # Save this data to now plot global picture
            if Run not in DataPhi.keys():
                DataPhi[Run] = [ [], [], [], [], [] ]

            # Save the individual B_i(rho,phi,z):
            DataPhi[Run][0].append( phi  )
            DataPhi[Run][1].append( Brho )
            DataPhi[Run][2].append( Bphi )
            DataPhi[Run][3].append( Bz   )
            DataPhi[Run][4].append( int(np.round(initz,0) ) )

            # Continue to the next phi cycle
            phiCyc += 2

    return DataPhi

def getZScans(RunList):
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
        if Run >= 2038 and Run <= 2041:
            zCyc = 7
        stopCyc = zCyc + 5

        # Get the phi offset for this Run and save the figure
        rhoOff,phiOff = getOffset(Run)
        ctr = 0
        while zCyc < stopCyc:
            scantype, initrho, initphi, initz, t,rho,phi,z,Brho,Bphi,Bz = grabCycleData(Run,zCyc)
            print(Run,zCyc,scantype)
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
            plt.ylim([-2500,2500])
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

            # Continue to the next phi cycle
            zCyc += 1
            ctr += 1
    return Data

def gatherData(AvgDataZ,Run):
    return AvgDataZ[Run][0],\
            AvgDataZ[Run][1],\
            AvgDataZ[Run][2],\
            AvgDataZ[Run][3]

def getSpread(Brho,Bphi,Bz):
    sBrho = np.abs(np.max(Brho) - np.min(Brho))
    sBphi = np.abs(np.max(Bphi) - np.min(Bphi))
    sBz   = np.abs(np.max(Bz  ) - np.min(Bz  ))
    return sBrho, sBphi, sBz


# We want to take the map results from (B0 up) + (B0 down) 
# and look at residual


RunList = [2038,2039,2040,2041]
DataPhi = getPhiRings(RunList)
DataZ = getZScans(RunList)


AvgDataZ = avgZData(DataZ)
U1_z,U1_Brho,U1_Bphi,U1_Bz = gatherData(AvgDataZ,2038)
D1_z,D1_Brho,D1_Bphi,D1_Bz = gatherData(AvgDataZ,2039)
D2_z,D2_Brho,D2_Bphi,D2_Bz = gatherData(AvgDataZ,2040)
U2_z,U2_Brho,U2_Bphi,U2_Bz = gatherData(AvgDataZ,2041)

plt.figure(0)
plt.plot(U1_z, (U1_Brho + D1_Brho)/2. , label='2038+2039')
plt.plot(U1_z, (U1_Brho + D2_Brho)/2. , label='2038+2040')
plt.plot(U2_z, (U2_Brho + D1_Brho)/2. , label='2041+2039')
plt.plot(U2_z, (U2_Brho + D2_Brho)/2. , label='2041+2040')
plt.axhline(y=0,linestyle='--',color='black',linewidth=2)
plt.xticks(fontsize=13)
plt.yticks(fontsize=13)
plt.xlabel(r'$z$ (cm)',fontsize=15)
plt.ylabel(r'$B_\rho$ (pT)',fontsize=15)
plt.legend(numpoints=1,loc='best',fontsize=12)
plt.xlim([-30,30])
plt.ylim([-100,100])
plt.tight_layout()
plt.savefig("sumB0_Brho_Zscans.pdf",bbox_inches='tight')
plt.close()

plt.figure(1)
plt.plot(U1_z, (U1_Bphi + D1_Bphi)/2. , label='2038+2039')
plt.plot(U1_z, (U1_Bphi + D2_Bphi)/2. , label='2038+2040')
plt.plot(U2_z, (U2_Bphi + D1_Bphi)/2. , label='2041+2039')
plt.plot(U2_z, (U2_Bphi + D2_Bphi)/2. , label='2041+2040')
plt.axhline(y=0,linestyle='--',color='black',linewidth=2)
plt.xticks(fontsize=13)
plt.yticks(fontsize=13)
plt.xlabel(r'$z$ (cm)',fontsize=15)
plt.ylabel(r'$B_\phi$ (pT)',fontsize=15)
plt.legend(numpoints=1,loc='best',fontsize=12)
plt.xlim([-30,30])
plt.ylim([-100,100])
plt.tight_layout()
plt.savefig("sumB0_Bphi_Zscans.pdf",bbox_inches='tight')
plt.close()

plt.figure(0)
plt.plot(U1_z, (U1_Bz + D1_Bz)/2. , label='2038+2039')
plt.plot(U1_z, (U1_Bz + D2_Bz)/2. , label='2038+2040')
plt.plot(U2_z, (U2_Bz + D1_Bz)/2. , label='2041+2039')
plt.plot(U2_z, (U2_Bz + D2_Bz)/2. , label='2041+2040')
plt.axhline(y=0,linestyle='--',color='black',linewidth=2)
plt.xticks(fontsize=13)
plt.yticks(fontsize=13)
plt.xlabel(r'$z$ (cm)',fontsize=15)
plt.ylabel(r'$B_z$ (pT)',fontsize=15)
plt.legend(numpoints=1,loc='best',fontsize=12)
plt.xlim([-30,30])
plt.ylim([-100,100])
plt.tight_layout()
plt.savefig("sumB0_Bz_Zscans.pdf",bbox_inches='tight')
plt.close()




U1_phi,U1_Brho,U1_Bphi,U1_Bz,U1_Z = makeAverages(DataPhi,2038)
D1_phi,D1_Brho,D1_Bphi,D1_Bz,D1_Z = makeAverages(DataPhi,2039)
D2_phi,D2_Brho,D2_Bphi,D2_Bz,D2_Z = makeAverages(DataPhi,2040)
U2_phi,U2_Brho,U2_Bphi,U2_Bz,U2_Z = makeAverages(DataPhi,2041)

Zs = []
U1_Spread_Brho = []; U1_Spread_Bphi = []; U1_Spread_Bz = []
U2_Spread_Brho = []; U2_Spread_Bphi = []; U2_Spread_Bz = []
D1_Spread_Brho = []; D1_Spread_Bphi = []; D1_Spread_Bz = []
D2_Spread_Brho = []; D2_Spread_Bphi = []; D2_Spread_Bz = []

for i in range(len(U1_Z)):
    plt.figure(i)
    plt.title(r'Phi-ring at $\rho=50,z=$'+str(U1_Z[i]),fontsize=17)
    plt.plot(U1_phi[i], (U1_Brho[i] + D1_Brho[i])/2. , label='2038+2039')
    plt.plot(U1_phi[i], (U1_Brho[i] + D2_Brho[i])/2. , label='2038+2040')
    plt.plot(U2_phi[i], (U2_Brho[i] + D1_Brho[i])/2. , label='2041+2039')
    plt.plot(U2_phi[i], (U2_Brho[i] + D2_Brho[i])/2. , label='2041+2040')
    plt.axhline(y=0,linestyle='--',color='black',linewidth=2)
    plt.xticks(fontsize=13)
    plt.yticks(fontsize=13)
    plt.xlabel(r'$\phi$ (deg.)',fontsize=15)
    plt.ylabel(r'$B_\rho$ (pT)',fontsize=15)
    plt.legend(numpoints=1,loc='best',fontsize=12)
    plt.xlim([0,360])
    plt.ylim([-200,200])
    plt.tight_layout()
    plt.savefig("sumB0_Brho_PhiRings_%i.pdf" % i,bbox_inches='tight')
    plt.close()


    plt.figure(i+5)
    plt.title(r'Phi-ring at $\rho=50,z=$'+str(U1_Z[i]),fontsize=17)
    plt.plot(U1_phi[i], (U1_Bphi[i] + D1_Bphi[i])/2. , label='2038+2039')
    plt.plot(U1_phi[i], (U1_Bphi[i] + D2_Bphi[i])/2. , label='2038+2040')
    plt.plot(U2_phi[i], (U2_Bphi[i] + D1_Bphi[i])/2. , label='2041+2039')
    plt.plot(U2_phi[i], (U2_Bphi[i] + D2_Bphi[i])/2. , label='2041+2040')
    plt.axhline(y=0,linestyle='--',color='black',linewidth=2)
    plt.xticks(fontsize=13)
    plt.yticks(fontsize=13)
    plt.xlabel(r'$\phi$ (deg.)',fontsize=15)
    plt.ylabel(r'$B_\phi$ (pT)',fontsize=15)
    plt.legend(numpoints=1,loc='best',fontsize=12)
    plt.xlim([0,360])
    plt.ylim([-200,200])
    plt.tight_layout()
    plt.savefig("sumB0_Bphi_PhiRings_%i.pdf" % i,bbox_inches='tight')
    plt.close()

    plt.figure(i+10)
    plt.title(r'Phi-ring at $\rho=50,z=$'+str(U1_Z[i]),fontsize=17)
    plt.plot(U1_phi[i], (U1_Bz[i] + D1_Bz[i])/2. , label='2038+2039')
    plt.plot(U1_phi[i], (U1_Bz[i] + D2_Bz[i])/2. , label='2038+2040')
    plt.plot(U2_phi[i], (U2_Bz[i] + D1_Bz[i])/2. , label='2041+2039')
    plt.plot(U2_phi[i], (U2_Bz[i] + D2_Bz[i])/2. , label='2041+2040')
    plt.axhline(y=0,linestyle='--',color='black',linewidth=2)
    plt.xticks(fontsize=13)
    plt.yticks(fontsize=13)
    plt.xlabel(r'$\phi$ (deg.)',fontsize=15)
    plt.ylabel(r'$B_z$ (pT)',fontsize=15)
    plt.legend(numpoints=1,loc='best',fontsize=12)
    plt.xlim([0,360])
    plt.ylim([-200,200])
    plt.tight_layout()
    plt.savefig("sumB0_Bz_PhiRings_%i.pdf" % i,bbox_inches='tight')
    plt.close()



    U1_sBrho,U1_sBphi,U1_sBz = getSpread( (U1_Brho[i] + D1_Brho[i])/2. ,(U1_Bphi[i] + D1_Bphi[i])/2. ,(U1_Bz[i] + D1_Bz[i])/2. )
    U2_sBrho,U2_sBphi,U2_sBz = getSpread( (U1_Brho[i] + D2_Brho[i])/2. ,(U1_Bphi[i] + D2_Bphi[i])/2. ,(U1_Bz[i] + D2_Bz[i])/2. )
    D1_sBrho,D1_sBphi,D1_sBz = getSpread( (U2_Brho[i] + D1_Brho[i])/2. ,(U2_Bphi[i] + D1_Bphi[i])/2. ,(U2_Bz[i] + D1_Bz[i])/2. )
    D2_sBrho,D2_sBphi,D2_sBz = getSpread( (U2_Brho[i] + D2_Brho[i])/2. ,(U2_Bphi[i] + D2_Bphi[i])/2. ,(U2_Bz[i] + D2_Bz[i])/2. )

    Zs.append( U1_Z[i] )
    U1_Spread_Brho.append( U1_sBrho ); U2_Spread_Brho.append( U2_sBrho )
    D1_Spread_Brho.append( D1_sBrho ); D2_Spread_Brho.append( D2_sBrho )

    U1_Spread_Bphi.append( U1_sBphi ); U2_Spread_Bphi.append( U2_sBphi )
    D1_Spread_Bphi.append( D1_sBphi ); D2_Spread_Bphi.append( D2_sBphi )

    U1_Spread_Bz.append( U1_sBz ); U2_Spread_Bz.append( U2_sBz )
    D1_Spread_Bz.append( D1_sBz ); D2_Spread_Bz.append( D2_sBz )
    
    
plt.figure(1)
plt.plot(Zs,U1_Spread_Brho,marker='o',linestyle='--',label='2038+2039',color='red')
plt.plot(Zs,U2_Spread_Brho,marker='o',linestyle='--',label='2038+2040',color='orange')
plt.plot(Zs,D1_Spread_Brho,marker='o',linestyle='--',label='2041+2039',color='green')
plt.plot(Zs,D2_Spread_Brho,marker='o',linestyle='--',label='2041+2040',color='magenta')
plt.grid(True)
plt.xticks(fontsize=13)
plt.yticks(fontsize=13)
plt.xlabel(r'$z$ (cm)',fontsize=15)
plt.ylabel(r'Spread of $B_\rho$ (pT)',fontsize=15)
plt.legend(numpoints=1,loc='best',fontsize=12)
plt.xlim([-32,32])
plt.ylim([0,300])
plt.tight_layout()
plt.savefig('spread_Brho.pdf',bbox_inches='tight')

plt.figure(2)
plt.plot(Zs,U1_Spread_Bphi,marker='o',linestyle='--',label='2038+2039',color='red')
plt.plot(Zs,U2_Spread_Bphi,marker='o',linestyle='--',label='2038+2040',color='orange')
plt.plot(Zs,D1_Spread_Bphi,marker='o',linestyle='--',label='2041+2039',color='green')
plt.plot(Zs,D2_Spread_Bphi,marker='o',linestyle='--',label='2041+2040',color='magenta')
plt.grid(True)
plt.xticks(fontsize=13)
plt.yticks(fontsize=13)
plt.xlabel(r'$z$ (cm)',fontsize=15)
plt.ylabel(r'Spread of $B_\phi$ (pT)',fontsize=15)
plt.legend(numpoints=1,loc='best',fontsize=12)
plt.xlim([-32,32])
plt.ylim([0,300])
plt.tight_layout()
plt.savefig('spread_Bphi.pdf',bbox_inches='tight')

plt.figure(3)
plt.plot(Zs,U1_Spread_Bz,marker='o',linestyle='--',label='2038+2039',color='red')
plt.plot(Zs,U2_Spread_Bz,marker='o',linestyle='--',label='2038+2040',color='orange')
plt.plot(Zs,D1_Spread_Bz,marker='o',linestyle='--',label='2041+2039',color='green')
plt.plot(Zs,D2_Spread_Bz,marker='o',linestyle='--',label='2041+2040',color='magenta')
plt.grid(True)
plt.xticks(fontsize=13)
plt.yticks(fontsize=13)
plt.xlabel(r'$z$ (cm)',fontsize=15)
plt.ylabel(r'Spread of $B_z$ (pT)',fontsize=15)
plt.legend(numpoints=1,loc='best',fontsize=12)
plt.xlim([-32,32])
plt.ylim([0,300])
plt.tight_layout()
plt.savefig('spread_Bz.pdf',bbox_inches='tight')

plt.show()

import numpy as np
from iminuit import Minuit
from iminuit.cost import LeastSquares
import matplotlib.pyplot as plt
import os.path

from map_tools import *
from map_cycle import *

# plot parameters
plt.rcParams["font.size"]=10
plt.rcParams["lines.markersize"]=3
plt.rcParams["figure.figsize"]=[12,8]
plt.rcParams["figure.dpi"]=200


def dumpRunVars(run, rn, filepath=None):
    if filepath==None:
        filepath = "data/vars_run{}".format(rn)
    C_data = run.C_data
    C_data_err = run.C_data_err
    G = run.G
    G_err = run.G_err
    Gs = run.Gs
    Gs_err = run.Gs_err
    Gsource = run.Gsource
    lmax = run.lmax
    runvars = np.array([rn, C_data, C_data_err, G, G_err, Gs, Gs_err, Gsource, lmax])
    if os.path.exists(filepath):
        print("filepath already exists, overwrite? (1/0)")
        yes = int(input())
        if yes==True:
            np.save(filepath, runvars, allow_pickle=True)
            print("file overwritten")
    else:
        np.save(filepath, runvars, allow_pickle=True)
        print("file written")
        
def loadRunVars(run, rn, filepath=None):
    if filepath==None:
        filepath = "data/vars_run{}".format(rn)
    runvars = np.load(filepath, allow_pickle=True)
    rn, C_data, C_data_err, G, G_err, Gs, Gs_err, Gsource, lmax = runvars
    run.C_data = C_data
    run.C_data_err = C_data_err
    run.G = G
    run.G_err = G_err
    run.Gs = Gs
    run.Gs_err = Gs_err
    run.Gsource = Gsource
    run.lmax = lmax
    print("vars restored")
    
def getPhantom(G, lmax=7):
    L = np.array([Geometry.L1, Geometry.L3, Geometry.L5, Geometry.L7])
    C = np.array([Geometry.C1, Geometry.C3, Geometry.C5, Geometry.C7])
    coefs = L/C
    phanG = np.array([ coefs[l] * G[2*l+1][lmax+1] for l in range((lmax-1)//2 + 1) ])
    return phanG


class Run:
    
    def __init__(self, run, cycrange=None, runtype=None, trim=False, calib=False, dtformat='v3', advcalib_nfit=0):
        """
        runtype = 'rings', 'zscans' or 'rhoscans'
        calcycles = <None> for no calibration
        """      
        self.dtformat = dtformat
        self.runnumber = run
        self.runtype = runtype
        self.logpath = '../maps_bin/{}_run.log'.format(str(run).zfill(6))
        if os.path.exists(self.logpath):
            with open(self.logpath) as file:
                lines = file.readlines()
                self.date_time = lines[0][1:17]
        cycles = {}
        cycNo = []
        cycTypes = {}
        self.Gs = {} # for keeping future Glms
        self.Gs_err = {}
        if runtype=='zscans':
            for cyc in cycrange:
                cycle = Cycle(run, cyc, dtformat=dtformat, trim=trim)
                if cycle.file_exists:
                    cycTypes[cyc] = cycle.scantype
                    if cycle.scantype=='z':
                        cycles[cyc] = cycle
                        cycNo.append(cyc)
            self.cycles = cycles
        elif runtype=='rhoscans':
            for cyc in cycrange:
                cycle = Cycle(run, cyc, dtformat=dtformat, trim=trim)
                if cycle.file_exists:
                    cycTypes[cyc] = cycle.scantype
                    if cycle.scantype=='rho':
                        cycles[cyc] = cycle
                        cycNo.append(cyc)
            self.cycles = cycles
        elif runtype=='rings':
            rings = []
            calcycles = {}
            for cyc in cycrange:
                cycle = Cycle(run, cyc, dtformat=dtformat, trim=trim)
                if cycle.file_exists:
                    cycTypes[cyc] = cycle.scantype
                    if cycle.scantype=='phi':
                        ring = cycle.ring
                        if ring not in rings:
                            cycles[ring] = cycle
                            rings.append(ring)
                        if ring==(0, 0):
                            calcycles[cyc] = cycle
                        cycNo.append(cyc)
            self.calcycles = calcycles
            self.rings = np.array(rings)
            self.cycles = cycles
            self.getRingsCoord()
            self.getRings_float()
            self.nPhi = len(self.cycles[self.ringsRho[0], self.ringsZ[0]].t)
        else:
            cycle_time = {}
            rings = []
            lines = []
            calcycles = {}
            calcycles2 = {}
            calcycles3 = {}
            cycles['rings'] = {}
            cycles['lines'] = {}
            for cyc in cycrange:
                cycle = Cycle(run, cyc, dtformat=dtformat, trim=trim)
                if cycle.file_exists:
                    cycTypes[cyc] = cycle.scantype
                    cycNo.append(cyc)
                    if cycle.scantype=='phi':
                        ring = cycle.ring
                        cycle_time[ring] = cycle.t_abs[0]
                        if ring not in rings:
                            cycles['rings'][ring] = cycle
                            rings.append(ring)
                        if ring==(0, 0):
                            calcycles[cyc] = cycle
                        if ring[0]==0:
                            calcycles2[ring] = cycle
                            calcycles3[cyc] = cycle
                    elif cycle.scantype=='z':
                        line = (round(cycle.rho[0]), round(cycle.phi[0]))
                        cycle_time[line] = cycle.t_abs[0]
                        if line not in lines:
                            cycles['lines'][line] = cycle
                            lines.append(line)
            # calibrate cycle timestamps
            t0 = cycle_time[min(cycle_time, key=cycle_time.get, default=(0, 0))]
            for k in cycle_time:
                cycle_time[k] -= t0
            self.cycle_time = cycle_time
            self.calcycles = calcycles
            self.calcycles2 = calcycles2
            self.calcycles3 = calcycles3
            self.rings = np.array(rings)
            self.lines = np.array(lines)
            self.cycles = cycles
            self.getRings_float()
        self.cycTypes = cycTypes
        self.cycNo = np.array(cycNo)
        if runtype!=None:
            self.nc = len(self.cycles)
            print("found {} cycles of '{}' scantype".format(self.nc, runtype))
        else:
            self.nrings = len(self.cycles['rings'])
            self.nlines = len(self.cycles['lines'])
            self.nc = self.nlines + self.nrings
            print("found {} cycles".format(self.nc))
        self.calibrated = False
        if calib==True and advcalib_nfit==0:
            self.calibrate()
        elif calib==True and advcalib_nfit!=0:
            self.calibrate_adv(advcalib_nfit)
        
    def saveVars(self):
        
        V = vars(self)
        oldVars = {}
        for var in V:
            oldVars[var] = getattr(self, var)
        self.oldVars = oldVars
        
        if self.runtype==None:
            for cyc in self.cycles['rings']:
                self.cycles['rings'][cyc].saveVars()
            for cyc in self.cycles['lines']:
                self.cycles['lines'][cyc].saveVars()
        else:
            cycles = self.cycles
            for cyc in cycles:
                cycles[cyc].saveVars()
            
    def restoreVars(self):
        
        try:
            oldVars = self.oldVars
        except AttributeError:
            print("no saved variables to restore")
            return 0
        
        for var in oldVars:
            setattr(self, var, oldVars[var])
            
        if self.runtype==None:
            for cyc in self.cycles['rings']:
                self.cycles['rings'][cyc].restoreVars()
            for cyc in self.cycles['lines']:
                self.cycles['lines'][cyc].restoreVars()
        else:
            cycles = self.cycles
            for cyc in cycles:
                cycles[cyc].restoreVars()
    
    def getRingsCoord(self):
        Rho = []
        Z = []
        rings = self.rings
        for (rho, z) in rings:
            if rho not in Rho:
                Rho.append(rho)
            if z not in Z:
                Z.append(z)
        self.ringsRho = np.array(Rho)
        self.ringsZ = np.array(Z)
        print("got rings")
        return self.ringsRho, self.ringsZ
    
    def getRings_float(self):
        if self.runtype==None:
            cycles = self.cycles['rings']
        else:
            cycles = self.cycles 
        rings = []
        for cyc in cycles:
            c = cycles[cyc]
            rings.append([c.R[0][0], c.R[2][0]])
        self.rings_float = np.array(rings)
        return self.rings_float
    
    def getMeanField(self):
        Bmean = {}
        Bstd = {}
        if self.runtype==None:
            cycles = self.cycles['rings']
        else:
            cycles = self.cycles 
        for cyc in cycles:
            Bmean[cyc] = np.array([ np.mean(cycles[cyc].B[i]) for i in range(3) ])
            Bstd[cyc] = np.array([ np.std(cycles[cyc].B[i]) for i in range(3) ])
        self.Bmean = Bmean
        self.Bstd = Bstd           
        return Bmean, Bstd
    
    def substractMeanField(self):
        if self.runtype==None:
            cycles = self.cycles['rings']
        else:
            cycles = self.cycles 
        self.getMeanField()
        self.saveVars()
        for cyc in cycles:
            B = cycles[cyc].B
            Bmean = self.Bmean[cyc]
            cycles[cyc].B = np.array([B[i]-Bmean[i] for i in range(3)])
            cycles[cyc].Brho, cycles[cyc].Bphi, cycles[cyc].Bz = cycles[cyc].B
        
    def arrayForm(self, dic):
        """
        converts a dictionary with ring keys to array of arrays
        """
        rho, z = self.ringsRho, self.ringsZ
        nr, nz = len(rho), len(z)
        dim = np.ndim(dic[(rho[0], z[0])])
        if dim==0:
            arr = np.array([ [dic[rho[i], z[j]] for j in range(nz)] for i in range(nr) ])
        elif dim==1:
            n = len(dic[(rho[0], z[0])])
            arr = np.array([ [ [dic[rho[i], z[j]][k] for j in range(nz)] for i in range(nr) ] for k in range(n)])
        return arr
    
    def calibrate(self):
        
        if self.calibrated==False:
            
            if self.runtype==None:
                ringCycles = self.cycles['rings']
                lineCycles = self.cycles['lines']
            else:
                cycles = self.cycles
                
            if self.dtformat=='v4':
                calcycles = self.calcycles3
                cal0, cal90, cal180  = 0, 0, 0
                for c in calcycles:
                    cyc = calcycles[c]
                    pinrot = int(cyc.pinmech[3][0])
                    if pinrot==0 and cal0==0:
                        cyc.cycFit(nfit=2)
                        Brho_off = cyc.par[0][0]
                        Bphi_off = cyc.par[1][0]
                        Bz_plus = cyc.par[2][0]
                        if cal90==0:
                            Bz_off = cyc.par[2][0]
                        cal0 = 1
                    if pinrot==90 and cal90==0:
                        cyc.cycFit(nfit=2)
                        Bz_off = cyc.par[2][0]
                        cal90 = 1
                    if pinrot==180:
                        cyc.cycFit(nfit=2)
                        Bz_minus = cyc.par[2][0]
                        cal180 = 1
                if cal180==1:
                    Bz_off = 0.5*(Bz_plus - Bz_minus)
                if cal90==0 and cal180==0:
                    print("/!\ B_z not properly calibrated /!\ ")
            else:
                calcycles = self.calcycles
                cs = []
                for c in calcycles:
                    cyc = calcycles[c]
                    cyc.cycFit(nfit=2)
                    cs.append(c)
                Brho_off = calcycles[cs[0]].par[0][0]
                Bphi_off = calcycles[cs[0]].par[1][0]
                if len(cs)>1:
                    Bz_off = calcycles[cs[1]].par[2][0]
                else:
                    Bz_off = calcycles[cs[0]].par[2][0]
                    print("/!\ B_z not properly calibrated /!\ ")
                    
            if self.runtype==None:
                for c in ringCycles:
                    B = ringCycles[c].B
                    ringCycles[c].B_nonCal = B
                    ringCycles[c].B = np.array([ B[0]-Brho_off, B[1]-Bphi_off, B[2]-Bz_off ])
                    ringCycles[c].Brho, ringCycles[c].Bphi, ringCycles[c].Bz = ringCycles[c].B
                for c in lineCycles:
                    B = lineCycles[c].B
                    lineCycles[c].B_nonCal = B
                    lineCycles[c].B = np.array([ B[0]-Brho_off, B[1]-Bphi_off, B[2]-Bz_off ])
                    lineCycles[c].Brho, lineCycles[c].Bphi, lineCycles[c].Bz = lineCycles[c].B
            else:
                for c in cycles:
                    B = cycles[c].B
                    cycles[c].B_nonCal = B
                    cycles[c].B = np.array([ B[0]-Brho_off, B[1]-Bphi_off, B[2]-Bz_off ])
                    cycles[c].Brho, cycles[c].Bphi, cycles[c].Bz = cycles[c].B
            print("calibrated")
            self.calibrated=True
            self.B_off = Brho_off, Bphi_off, Bz_off
            
        else:
            print("already calibrated")
    
    def calibrate_adv(self, nfit=0):
        
        calcycles = self.calcycles2
        if self.runtype==None:
            cycles = self.cycles['rings']
            print("calibration method unfinished for lines")
        else:
            print("/!\ not calibrated")
            print("calibration method unfinished for ringtype runs")
            return None
        
        # fluxgate offsets from rho=0 rings
        Brho_off, Bphi_off = {}, {}
        for c in calcycles:
            t = round(self.cycle_time[c])
            cyc = cycles[c]
            cyc.cycFit(nfit=nfit)
            Brho_off[t] = cyc.par[0][0]
            Bphi_off[t] = cyc.par[1][0]
        Brho_off = dict(sorted(Brho_off.items()))
        Bphi_off = dict(sorted(Bphi_off.items()))

#         print(Brho_off)
#         print(Bphi_off)
        
        # fluxgate slopes
        time = [*Brho_off]
        nt = len(time)
        Brho_slope, Bphi_slope = {}, {}
        for it in range(nt-1):
            t0 = time[it]
            t1 = time[it+1]
            Brho_slope[t0] = (Brho_off[t1]-Brho_off[t0]) / (t1-t0)
            Bphi_slope[t0] = (Bphi_off[t1]-Bphi_off[t0]) / (t1-t0)

        print(time)
#         print(Brho_slope)
#         print(Bphi_slope)
            
        # /!\ temporary Bz fix
        cycles[0, 0].cycFit(nfit=2)
        zoff = cycles[0, 0].par[2][0]
        # calibration with slope from nearest points
        for c in cycles:
            ct = self.cycle_time[c]
            it = 1
#             print(c, ct)
            while ct>=time[it]:
                it += 1
                if it==(nt-1): #avoid checking while condition if index reaches the end of "time" (can't use or in while)
                    break
#             print(it)
            st = time[it-1]
#             print(st)
            #st = time[np.argmin(abs(time[:-1]-ct))]
            B = cycles[c].B
            roff = Brho_off[st] + Brho_slope[st]*(ct - st)
            poff = Bphi_off[st] + Bphi_slope[st]*(ct - st)
            cycles[c].B_nonCal = B
            cycles[c].B = np.array([ B[0]-roff, B[1]-poff, B[2]-zoff ])
            cycles[c].Brho, cycles[c].Bphi, cycles[c].Bz = cycles[c].B
#             print(roff, poff)
    
    # deprecated
    def calibrate_old(self, calcycles):
        """
        calcycles = ('calcycle_rho', 'calcycle_phi', 'calcycle_z')
        """
        cycles = self.cycles
        try:
            Bmean = self.Bmean
            Bstd = self.Bstd
        except AttributeError:
            Bmean, Bstd = self.getMeanField()
        
        if self.calibrated==False:
            B0_rho = Bmean[calcycles[0]][0]
            B0_phi = Bmean[calcycles[1]][1]
            if len(calcycles)>2:
                B0_z = Bmean[calcycles[2]][2]
            else: 
                B0_z = Bmean[min(abs(self.ringsRho)), min(abs(self.ringsZ))][2]
            for cyc in cycles:
                B = cycles[cyc].B
                cycles[cyc].B_nonCal = B
                cycles[cyc].B = np.array([ B[0]-B0_rho, B[1]-B0_phi, B[2]-B0_z ])
            self.calibrated = True
            self.getMeanField()
            print("calibrated")
            return 1
        else:
            print("already calibrated")
            return 0
    
    def getPars(self, nfit):
        if self.runtype==None:
            cycles = self.cycles['rings']
        else:
            cycles = self.cycles 
        Par = []
        Err = []
        for cyc in cycles:
            cycle = cycles[cyc]
            cycle.cycFit(bounds=None, nfit=nfit, norm=True)
            Par.append(cycle.par)
            Err.append(cycle.err)
        self.pars = np.array(Par)
        self.errs = np.array(Err)
        self.getRingCoefs()
        return self.pars, self.errs
    
    def getRingCoefs(self):
        try:
            pars, errs = self.pars, self.errs
        except AttributeError:
            print("get ring fit parameters first")
            return None
        if self.runtype==None:
            cycles = self.cycles['rings']
        else:
            cycles = self.cycles 
        C = {}
        C_err = {}
        for cyc in cycles:
            cycle = cycles[cyc]
            C[cyc] = cycle.par
            C_err[cyc] = cycle.err
        self.C_data = C
        self.C_data_err = C_err
        return C, C_err
        
    def getGl0_lines(self, nfit=9):
        return None
    
    def getGl0(self, nfit=9):
        try:
            pars, errs = self.pars, self.errs
        except AttributeError:
            pars, errs = self.getPars(nfit)
        if self.runtype==None:
            nc = self.nrings
        else:
            nc = self.nc
        self.G0 = np.array([ [pars[i][2][l] for i in range(nc)] for l in range(nfit+1) ])
        self.G0err = np.array([ [errs[i][2][l] for i in range(nc)] for l in range(nfit+1) ])
            
    
    def getParvst(self, ipar, iB, nfit=6): # old, replace?
        try:
            pars = self.pars
        except AttributeError:
            self.getPars(nfit)
            pars = self.pars
        if self.runtype==None:
            nc = self.nrings
        else:
            nc = self.nc
        par = np.array([ pars[i][iB][ipar] for i in range(nc) ])
        t = [0]
        if self.runtype==None:
            cycles = self.cycles['rings']
        else:
            cycles = self.cycles 
        for cyc in cycles:
            cycle = cycles[cyc]
            dt = cycle.t[-1]-cycle.t[0]
            t.append(t[-1] + dt)
        t = np.array(t[:-1])
        return par, t
    
    def fourierCoef(self, iB, m, nfit=4):
        try:
            pars = self.pars
            errs = self.errs
        except AttributeError:
            self.getPars(nfit=nfit)
            pars = self.pars
            errs = self.errs
        if self.runtype==None:
            nc = self.nrings
        else:
            nc = self.nc
        a = np.array([ pars[i][iB][m] for i in range(nc) ])
        a_err = np.array([ errs[i][iB][m] for i in range(nc) ])
        return a, a_err     
            
    def getAm(self, iB, m, lmax=3):
        m0 = Tools.m_to_m(m, iB, lmax, mode='positive')
        a, a_err = self.fourierCoef(iB, m0, nfit=lmax+1)
        return a, a_err
    
    def getGm(self, iB, m, lmin=0, lmax=3):
        try:
            rings = self.rings_float
        except AttributeError:
            self.getRings()
            rings = self.rings_float
        a, a_err = self.getAm(iB=iB, m=m, lmax=lmax)
        #print(a, a_err, rings, m, lmin, lmax)
        g, g_err = Fits.GlmPoly(iB, a, a_err, rings, m, lmin, lmax)
        return g, g_err

    def getGlm(self, lmax, source='rho', norm=True):
        """
        source = "rho", "phi" or "z"
        """
        try:
            G = self.G
            G_err = self.G_err
            source0 = self.Gsource
            lmax0 = self.lmax
            if lmax==lmax0 and source==source0:
                return G, G_err
            else:
                self.getPars(nfit=lmax+1)
        except AttributeError:
            self.getPars(nfit=lmax+1)
                
        G = np.zeros((lmax+1, 2*lmax+3))
        G_err = np.zeros((lmax+1, 2*lmax+3))

        for m in range(-lmax-1, lmax+2):
            #print("m = {}".format(m))
            if m==0:
                if source=='rho':
                    lmin=1
                    iB=0
                else:
                    lmin=0
                    iB=2
            else:            
                if source=='rho':
                    lmin = abs(m)-1
                    iB=0
                elif source=='phi':
                    lmin = abs(m)-1
                    iB=1
                else:
                    lmin = abs(m)
                    iB=2

            g, g_err = self.getGm(iB=iB, m=m, lmin=lmin, lmax=lmax)
            #print(g)
            for l in range(lmin, lmax+1):
                G[l][lmax+1+m] = g[l-lmin]
                G_err[l][lmax+1+m] = g_err[l-lmin]
                
        self.G, self.G_err = G, G_err
        self.Gsource = source
        self.lmax = lmax
        print("got G, G_err")
        if norm==True:
            self.normGlm()
            
        self.Gs[source] = G
        self.Gs_err[source] = G_err
        self.getFittedCoefs()
        self.getCycleFit()
            
        return self.G, self.G_err
    
    def combineGlm(self, mode='z'):
        """
        mode = 'rho' or 'z'
        combines Glm from mainly rho or z source to give a full set of Glms
        saves the combination to self.Gs["mainly_{}".format(mode)]
        """
        try:
            Grho, Grho_err = self.Gs["rho"], self.Gs_err["rho"]
            Gz, Gz_err = self.Gs["z"], self.Gs_err["z"]
            lmax = self.lmax
        except KeyError:
            print("get Glm with both rho and z sources first")
            print(help(self.combineGlm))
            return 0
            
        newG = np.zeros((lmax+1, 2*lmax+3))
        newG_err = np.zeros((lmax+1, 2*lmax+3))
        for l in range(0, lmax+1):
            for m in range(0, 2*lmax+3):
                ## relative m
                m0 = m-lmax-1
                if m0==0:
                    if l==0:
                        newG[l][m] = Gz[l][m]
                        newG_err[l][m] = Gz_err[l][m]
                    else:
                        if mode=='z':
                            newG[l][m] = Gz[l][m]
                            newG_err[l][m] = Gz_err[l][m]
                        elif mode=='rho':
                            newG[l][m] = Grho[l][m]
                            newG_err[l][m] = Grho_err[l][m]
                elif abs(m0)>l:
                    newG[l][m] = Grho[l][m]
                    newG_err[l][m] = Grho_err[l][m]
                else:
                    if mode=='z':
                        newG[l][m] = Gz[l][m]
                        newG_err[l][m] = Gz_err[l][m]
                    elif mode=='rho':
                        newG[l][m] = Grho[l][m]
                        newG_err[l][m] = Grho_err[l][m]
        
        source = "mainly_{}".format(mode)
        self.G, self.G_err = newG, newG_err
        self.Gsource = source
        print("Glm combined according to {} mode".format(source))
        
        self.Gs[source] = newG
        self.Gs_err[source] = newG_err
        self.getFittedCoefs()
        self.getCycleFit()
        
        return self.G, self.G_err
            
    
    def normGlm(self):
        
        try:
            G = self.G
            G_err = self.G_err
            source = self.Gsource
        except AttributeError:
            print("get Glm first")
            return 0
        if self.runtype==None:
            nc = self.nrings
        else:
            nc = self.nc
        
        lmax = np.shape(self.G)[0] - 1
        A_data = []
        A_err = []
        A_model = []
        
        for m in range(-lmax-1, lmax+2):
            #print("m = {}".format(m))
            if m==0:
                lmin = 0
                iB=2
            else:
                lmin = abs(m)-1
                if source=='rho':
                    iB=0
                elif source=='phi':
                    iB=1
                elif source=='z':
                    if m==-lmax-1 or m==lmax+1: # default for boundary m values is Brho
                        iB=0
                    else:
                        iB=2
                        
            a_model = self.fourierCoefFit(iB=iB, m=m, lmax=lmax)
            a_data, a_err = self.getAm(iB=iB, m=m, lmax=lmax)
            redchi2 = Tools.chi2stat(a_data, a_model, a_err) / nc
            
            for l in range(lmin, lmax):
                G_err[l][lmax+1+m] *= np.sqrt(redchi2)
                
        self.G_err = G_err
        print("G_err normalized")
        
        return G, G_err
    
    def fourierCoefFit(self, iB, m, lmax):
        if self.runtype==None:
            nc = self.nrings
        else:
            nc = self.nc
        a_model = np.zeros(nc)
        if m==0:
            lmin=0
        else:
            lmin = abs(m)-1
        for i in range(nc):
            ring = self.rings_float[i]
            for l in range(lmin, lmax+1):
                a_model[i] += self.G[l][lmax+1+m] * Harmonic.PiI(iB, ring[0], ring[1], l=l, m=m)
        return a_model       
    
    def getFittedCoefs(self):
        try:
            lmax = self.lmax
        except AttributeError:
            self.getGlm(lmax=3)
            lmax = self.lmax
            
        if self.runtype==None:
            cycles = self.cycles['rings']
            nc = self.nrings
        else:
            cycles = self.cycles
            nc = self.nc
        
        c_grid = np.zeros((2*lmax+3, 3, nc))
        for m in range(2*lmax+3):
            for iB in range(3):
                # array of each cycle's B[iB] m-fourier coefficient from Glm fit
                m0 = Tools.m_to_m(m, iB, lmax, mode='relative')
                #print("m = {}".format(m))
                #print("m0 = {}".format(m0))
                c_grid[m][iB] = self.fourierCoefFit(iB, m0, lmax)
                    
        C_model = {}
        for ic in range(nc):
            ring = self.rings[ic]
            C_model[ring[0], ring[1]] = np.array([ [c_grid[m][iB][ic] for m in range(2*lmax+3)] for iB in range(3) ])
            
        self.c_grid = c_grid
        self.C_model = C_model
        
        return C_model
    
    def fourierCoefPlot(self, iB, m, var='rho', plane=0):
        """
        var = 'rho' or 'z'
        plane = rho or z plane coordinate (int)
        """
        if self.runtype==None:
            cycles = self.cycles['rings']
        else:
            cycles = self.cycles
        A_fit = []
        A_dat = []
        x = []
        if var=='rho':
            for cyc in cycles:
                if cyc[1]==plane:
                    x.append(cyc[0])
                    A_fit.append(self.C_model[cyc][iB][m])
                    A_dat.append(self.C_data[cyc][iB][m])                    
        else:
            print("only implemented var='rho'")
        x = np.array(x)
        A_dat = np.array(A_dat)
        A_fit = np.array(A_fit)
        return x, A_dat, A_fit
        
    def fourierCoefFit_g(self, g, iB, m, lmin, lmax): #old
        a_model = np.zeros(self.nc)
        for i in range(self.nc):
            ring = self.rings_float[i]
            for l in range(lmin, lmax):
                a_model[i] += g[l-lmin] * Harmonic.PiI(iB, ring[0], ring[1], l=l, m=m)
        return a_model

    def getCycleFit(self):
        try:
            C_model = self.C_model
            lmax = self.lmax
            if len(C_model[[*C_model][0]][0])!=(2*lmax+3):
                self.getFittedCoefs()
                C_model = self.C_model
                lmax = self.lmax
        except AttributeError:
            self.getFittedCoefs()
            C_model = self.C_model
            lmax = self.lmax
        
        if self.runtype==None:
            cycles = self.cycles['rings']
        else:
            cycles = self.cycles
            
        for cyc in cycles:
            cycle = cycles[cyc]
            cycle.c_model = C_model[cyc]
            a_model = np.array([ [C_model[cyc][iB][m] for m in range(lmax+2)] for iB in range(3) ])
            b_model = np.array([ [C_model[cyc][iB][m] for m in range(lmax+2, 2*lmax+3)] for iB in range(3) ])
            cycle.B_model = np.array([ Functions.fourier(cycle.phi, a_model[iB], b_model[iB], T=360)  for iB in range(3) ])
            
        print("fitted B in self.cycles[cycle].B_model")           
    
    def getFitRMSs(self, nfit=4):
        RMS = {}
        if self.runtype==None:
            cycles = self.cycles['rings']
        else:
            cycles = self.cycles 
        for cyc in cycles:
            cycle = cycles[cyc]
            RMS[cyc] = cycle.getFitRMS(nfit=nfit)
        self.RMS = RMS
        return RMS
    
    def getPolyFitRMS(self):
        self.getCycleFit()
        cycles = self.cycles['rings']
        Bres = {}
        RMS = {}
        for c in cycles:
            bres = np.array([cycles[c].B_model[i] - cycles[c].B[i] for i in range(3)])
            RMS[c] = np.array([np.mean(bres[i]**2)**0.5 for i in range(3)])
            Bres[c] = bres
        self.RMS = RMS
        self.Bres = Bres
        return Bres, RMS
    
    def ringPlot(self, data='RMS', iB=0, nfit=4, sf=30, cmap='magma'):
        """
        data = 'RMS', 'meanB'
        """
        if data=='RMS':
            try:
                RMS = self.RMS
            except AttributeError:
                RMS = self.getFitRMSs(nfit=nfit)
            Z = self.arrayForm(RMS)[iB]
            label = 'RMS'
        elif data=='meanB':
            try:
                meanB = self.Bmean
                stdB = self.Bstd
            except AttributeError:
                meanB, stdB = self.getMeanField(nfit=nfit)
            Z = self.arrayForm(meanB)[iB]
            label = 'mean field'
        else:
            print("wrong data. data = 'RMS', 'meanB'")
    
        color = np.transpose(Z)
        size = (sf * np.ones_like(Z))**2
        x, y = np.meshgrid(self.ringsRho, self.ringsZ)
        
        ax = plt.scatter(x, y, s=size, c=color, alpha=0.5, cmap=cmap)
        cb = plt.colorbar(ax)
    
        plt.xlabel(r"$\rho$ (cm)", size=15)
        plt.ylabel(r"$z$ (cm)", size=15)
        if iB==0:
            cb.set_label(r"$B_\rho$ {} (pT)".format(label), size=15)
        elif iB==1:
            cb.set_label(r"$B_\varphi$ {} (pT)".format(label), size=15)
        elif iB==2:
            cb.set_label(r"$B_z$ {} (pT)".format(label), size=15)
        plt.grid()
        plt.show()
    
    def runPlot(self, iB=2, cmap='magma', zratio=6, elev=15, azim=30):
        
        if self.runtype==None:
            cycles = self.cycles['rings']
        else:
            cycles = self.cycles 
        nphi = self.nPhi

        Rho, Phi, Z = np.array([ [ [cycles[cyc].R[iR][iphi] for iphi in range(nphi)] for cyc in cycles] for iR in range(3) ])
        Brho, Bphi, Bz = np.array([ [ [cycles[cyc].B[iB][iphi] for iphi in range(nphi)] for cyc in cycles] for iB in range(3) ])
        
        X = Rho*np.cos(Phi*180/np.pi)
        Y = Rho*np.sin(Phi*180/np.pi)
        
        fig = plt.figure(figsize=plt.figaspect(1))
        ax = fig.add_subplot(projection='3d')
        
        if iB==0:
            p = ax.scatter(X, Y, Z, c=Brho, cmap=cmap)
            cb = fig.colorbar(p, ax=ax, shrink=0.5)
            cb.set_label(r"$B_\rho$ (pT)", size=15)
        elif iB==1:
            p = ax.scatter(X, Y, Z, c=Bphi, cmap=cmap)
            cb = fig.colorbar(p, ax=ax, shrink=0.5)
            cb.set_label(r"$B_\varphi$ (pT)", size=15)
        elif iB==2:
            p = ax.scatter(X, Y, Z, c=Bz, cmap=cmap)
            cb = fig.colorbar(p, ax=ax, shrink=0.5)
            cb.set_label(r"$B_z$ (pT)", size=15)

        ax.set_xlabel(r"$x$ (cm)", size=12)
        ax.set_ylabel(r"$y$ (cm)", size=12)
        ax.set_zlabel(r"$z$ (cm)", size=12)

        ax.view_init(elev, azim)
        ax.set_box_aspect((4, 4, zratio))
        #plt.show()
        
        return ax, fig
    
    def magPos(self, iR, iB):
        """
        z scans only
        """
        zc = []
        if self.runtype==None:
            cycles = self.cycles['rings']
        else:
            cycles = self.cycles 
        cn = self.cycNo
        for i in range(cn[0], cn[-1], 2):
            #x_do, B_do, Bmodel_do = cycles[i].getBfit(nfit=4)
            #x_up, B_up, Bmodel_up = cycles[i+1].getBfit(nfit=4)
            #L = min(len(Bmodel_do[iB]), len(Bmodel_up[iB]))
            #deltaB = Bmodel_do[iB][:L-1] - Bmodel_up[iB][:L-1]
            #i_zc = Tools.zero_crossing(deltaB, up=True)
            #zc.append(0.5 * (cycles[i].R[iR][i_zc-1] + cycles[i].R[iR][i_zc]))
            L = min(len(cycles[i].B[2]), len(cycles[i+1].B[2]))
            z = cycles[i].R[2][:L-1]
            Bdiff = 0.5 * ( cycles[i].B[2][:L-1] - cycles[i+1].B[2][:L-1] )
            par, err = Fits.polynomial(Bdiff, z, N=1)
            zc.append(-par[0]/par[1])
        zc = np.array(zc)
        return np.mean(zc), np.std(zc)
    
    def rhoZTF(self, F, argsF, G, argsG):
        """
        Consider new coordinate system rho2 = (rho, F1(rho)), z2 = (G(z), z)
        with F and G given functions
        Applies transformation of coordinates (rho, z) --> (rho + G(z), z + F(rho))
        to all cycle motor positions
        """
        
        self.saveVars()
        if self.runtype==None:
            cycles = self.cycles['rings']
        else:
            cycles = self.cycles 
        
        for cyc in cycles:
            
            cycle = cycles[cyc]
            cycle.saveVars()
            
            rho = cycle.R[0]
            z = cycle.R[2]
            
            rho2 = rho + G(z, argsG)
            z2 = z + F(rho, argsF)
            
            cycle.R[0] = rho2
            cycle.R[2] = z2
            
        self.getRings_float()
        
        return self.rings_float
    
    def fieldTilt(self, F_alpha, F_beta=lambda x : 0, F_gamma=None, scan='rho'):
        """
        F_alpha, F_beta, F_gamma are tilt functions that take scan coordinates DURING 'rho' or 'z' scan and return angles (in degrees) of the magnetic field axes, resp around rho, phi, z axes. 
        Takes angle-arrays where each entry maps from a mapper position.
        Transforms field projections accordingly.
        """
        try:
            getattr(self, transformed)
        except NameError:
            self.transformed=False
        if self.transformed==False:
            self.saveVars()
            
        if self.runtype==None:
            cycles = self.cycles['rings']
        else:
            cycles = self.cycles 
        
        for cyc in cycles:
            c = cycles[cyc]
            
            Brho, Bphi, Bz = c.B
            if scan=='rho':
                cx = c.rho
            elif scan=='z':
                cx = c.z
            alpha = F_alpha(cx) * np.pi / 180 # in radians
            beta = F_beta(cx) * np.pi / 180 # in radians
            #gamma = F_gamma(cx) * np.pi / 180 # in radians
            
            # at first order (?)
            Brho2 = Brho*np.cos(beta) + Bz*np.sin(beta)
            Bphi2 = Bphi*np.cos(alpha) + Bz*np.sin(alpha)
            Bz2 = Bz*np.cos(alpha)*np.cos(beta) - Bphi*np.sin(alpha) - Brho*np.sin(beta)
            
            c.Brho, c.Bphi, c.Bz = Brho2, Bphi2, Bz2
            c.B = c.Brho, c.Bphi, c.Bz
        
        self.transformed=True
        print("field tilted")
        
    
    def coordTF(self, F1, F2, G0, G2, H0, H1):
        """
        Consider new coordinate system x2 = (x, F1(x), F2(x)), y2 = (G0(y), y, G2(y)), z2 = (H0(z), H1(z), z)
        with F, G, H functions of the form F = lambda x : Functions.polynomial(x, args)
        Applies transformation of coordinates (x, y, z) --> (x + G0(y) + H0(z), y + F1(x) + H1(z), z + F2(x) + G2(y))
        to all motor positions
        """
        try:
            getattr(self, transformed)
        except NameError:
            self.transformed = False
        if self.transformed==False:
            self.saveVars()
            
        if self.runtype==None:
            cycles = self.cycles['rings']
        else:
            cycles = self.cycles 
        
        for cyc in cycles:
            c = cycles[cyc]
            
            # cartesian coordinates in (u_rho, u_phi, u_z) frame
            x, y, z = Tools.cartesian(c.rho, 0, c.z)
            
            x2 = x + G0(y) + H0(z)
            y2 = y + F1(x) + H1(z)
            z2 = z + F2(x) + G2(y)
            
            # back to cylindrical coordinates (angle resetted)
            rho, phi, z = Tools.cylindrical(x2, y2, z2)
            c.rho, c.phi, c.z = rho, c.phi + phi, z
            c.R = c.rho, c.phi, c.z
            
        self.getRings_float()
        self.transformed=True
        
        return self.rings_float
    
    def plotGm(self, l, lmax=3, w=0.1):
    
        try:
            G = self.G
            G_err = self.G_err
            lmax = self.lmax
        except AttributeError:
            G, G_err = self.getGlm(lmax=lmax)
        
        fig, ax = plt.subplots(1, 1)
        M = np.arange(-l-1, l+2)

        for m in M:
            ax.bar(m, G[l][lmax+1+m], yerr=G_err[l][lmax+1+m], align='center', alpha=0.5, ecolor='black', capsize=5, color='tab:blue', width=w)

        ax.grid()
        ax.set_xticks(M)
        ax.set_xlabel(r"$m$", size=15)
        if l==0:
            ax.set_ylabel(r"$G_{{ {}m }}$ (pT)".format(l), size=15)
        elif l==1:
            ax.set_ylabel(r"$G_{{ {}m }}$ (pT/cm)".format(l), size=15)
        elif l>1:
            ax.set_ylabel(r"$G_{{ {}m }}$ (pT/cm$^{{{}}}$)".format(l, l), size=15)

        return fig, ax
    
        
    def vecFieldPlot(self, cycles=None, di0=10, zratio=4, elev=20, azim=30, arrow_length=1, arrow_width=1.5, headbody_ratio=0.4, normalize=True, cmap='magma'):

        if cycles==None:
            if self.runtype==None:
                cycles = self.cycles['rings']
            else:
                cycles = self.cycles 
            
        fig = plt.figure(figsize=plt.figaspect(1))
        ax = fig.add_subplot(projection='3d')

        for cyc in cycles:

            cycle = self.cycles[cyc]
            
            bo = Tools.slopeChange(cycle.phi)
            ring_rho, ring_z = cyc 
            if ring_rho<50:
                di = 2*di0
            else:
                di = di0
            imax = bo[0] - di

            Rc = Tools.cartesian(cycle.rho, cycle.phi, cycle.z)
            Bc = Tools.cartesianField(cycle.phi, cycle.Brho, cycle.Bphi, cycle.Bz)

            # average over surrounding points
            x, y, z = np.array([[Rc[k][i] for i in range(0, imax, di)] for k in range(3)]) 
            Bx, By, Bz = np.array([[np.mean(Bc[k][i:i+di]) for i in range(0, imax, di)] for k in range(3)])

            norm = np.sqrt(Bx**2 + By**2 + Bz**2)
            head = np.array([ norm[i//2] for i in range(2*len(norm)) ])
            colors = np.concatenate([norm, head])

            q = ax.quiver(x, y, z, Bx, By, Bz, cmap=cmap, length=arrow_length, arrow_length_ratio=headbody_ratio, alpha=0.8, linewidths=arrow_width, normalize=normalize)
            q.set_array(colors)

        cb = fig.colorbar(q, ax=ax, shrink=0.5)
        cb.set_label(r"$|B|$ (pT)", size=15)

        ax.set_xlabel(r"$x$ (cm)", size=12)
        ax.set_ylabel(r"$y$ (cm)", size=12)
        ax.set_zlabel(r"$z$ (cm)", size=12)

        ax.set_zlim(-50, 50)

        ax.view_init(elev, azim)
        ax.set_box_aspect((4, 4, zratio))
        
        return fig, ax

        #plt.savefig('plots/vecField_0.png', dpi=200)
        
    def getFourierCoefs(self, ab, n, T=360):
        C = {}
        if self.runtype==None:
            cycles = self.cycles['rings']
        else:
            cycles = self.cycles 
        for cyc in cycles:
            coef = cycles[cyc].getFourierCoef(ab, n, T)
            C[cyc] = coef
        return C
    
    def gridding(self, rho0, iB=0, scans='all'):
        # scans = 'rings', 'lines', 'all'
        rings = self.rings
        if self.runtype==None:
            lines = self.lines
            ring_cycles = self.cycles['rings']
            line_cycles = self.cycles['lines']
        else:
            scans='rings'
            ring_cycles = self.cycles
        z_range = []
        for ring in rings:
            if ring[0]==rho0:
                z_range.append(ring[1])
        phi_range = []
        if self.runtype==None:
            for line in lines:
                if line[0]==rho0 and line[1]>=0:
                    phi_range.append(line[1])
            phi_len = len(phi_range)
            phi_step = (max(phi_range) - min(phi_range)) / (phi_len-1)
        phi_range = np.array(phi_range)
        z_range = np.array(z_range)
        z_len = len(z_range) 
        z_step = (max(z_range) - min(z_range)) / (z_len-1)

        Phi, Z = np.meshgrid(phi_range, z_range)

        ### initiate (phi, z) grid
        vsurf = {}
        for phi in phi_range:
            for z in z_range:
                vsurf[phi, z] = []

        ### fill (phi, z) grid with list of field values inside each cell
        if scans=='rings' or scans=='all':
            for ring in ring_cycles:
                if ring[0]==rho0:
                    z0 = ring[1]
                    phi = ring_cycles[ring].phi
                    B = ring_cycles[ring].B[iB]
                    for i in range(len(phi)):
                        for j in range(phi_len):
                            if phi[i]>(phi_range[j]-0.5*phi_step) and phi[i]<(phi_range[j]+0.5*phi_step):
                                vsurf[phi_range[j], z0].append(B[i])
        if scans=='lines' or scans=='all':
            for line in line_cycles:
                if line[0]==rho0 and line[1]>=0:
                    phi0 = line[1]
                    z = line_cycles[line].z
                    B = line_cycles[line].B[iB]
                    for i in range(len(z)):
                        for j in range(z_len):
                            if z[i]>(z_range[j]-0.5*z_step) and z[i]<(z_range[j]+0.5*z_step):
                                vsurf[phi0, z_range[j]].append(B[i])
        
        ### take average of each cell
        vsurf_avg = {}
        for cell in vsurf:
            vsurf_avg[cell] = np.mean(vsurf[cell])

        ### construct field meshgrid with vsurf
        B = np.zeros((z_len, phi_len))
        for i in range(z_len):
            for j in range(phi_len):
                B[i][j] = vsurf_avg[phi_range[j], z_range[i]]

        return Phi, Z, B
    
    
    def getSpectrum(self, metric='false edm'):
        """
        takes: metric ('false edm' or 'non-uniformity')
        returns: spectrum on metric as array
        """
        try:
            G = self.G
            G_err = self.G_err
            lmax = self.lmax
        except AttributeError:
            print("extract Glm first")
            return None
        
        if metric=='false edm':
            
            L = np.array([Geometry.L1, Geometry.L3, Geometry.L5, Geometry.L7])
            C = np.array([Geometry.C1, Geometry.C3, Geometry.C5, Geometry.C7])
            coefs = L/C
            
            G1 = np.array([ coefs[l] * G[2*l+1][lmax+1] for l in range((lmax-1)//2 + 1) ])
            G1_err = np.array([ coefs[l] * G_err[2*l+1][lmax+1] for l in range((lmax-1)//2 + 1) ])
            
            self.phantomG = G1
            self.phantomG_err = G1_err
            
            return G1, G1_err
            
        elif metric=='non-uniformity':
            
            lmax0=2
            print("non-uniformity metric currently up to l=2")
            
            Gz = np.zeros((lmax+1, 2*lmax+3))
            Gz_err = np.zeros((lmax+1, 2*lmax+3))
            
            for l in range(lmax+1):
                for m in range(-lmax-1, lmax+2):
                    # for now stops at lmax0
                    if l<=lmax0 and abs(m)<=(lmax0+1):
                        bz = Geometry.nonUniformity(l, m)**0.5
                        Gz[l][lmax+1+m] = G[l][lmax+1+m] * bz
                        Gz_err[l][lmax+1+m] = G_err[l][lmax+1+m] * bz
                    
            return Gz, Gz_err
        
        else:
            print("usage error:")
            print(help(self.getSpectrum))
            return None
        
    def plotSpectrum(self, metric='false edm', w=0.1, plot_l=1):
        
        try:
            G = self.G
            G_err = self.G_err
            lmax = self.lmax
        except AttributeError:
            print("extract Glm first")
            return None
        
        if metric=='false edm':
            
            G1, G1_err =  self.getSpectrum(metric)
            
            fig, ax = plt.subplots(1, 1)
            L = np.arange(1, lmax+1, 2)

            for l in L:
                ax.bar(l, G1[l//2], yerr=G1_err[l//2], align='center', alpha=0.5, ecolor='black', capsize=5, color='tab:blue', width=w)

            ax.grid()
            ax.set_xticks(L)
            ax.set_xlabel(r"$l$", size=15)
            ax.set_ylabel(r"Effective gradient $G_{l}$ (pT/cm)", size=15)
            
            return fig, ax
        
        elif metric=='non-uniformity':
            
            Gz, Gz_err =  self.getSpectrum(metric)
            l=plot_l

            fig, ax = plt.subplots(1, 1)
            M = np.arange(-l-1, l+2)

            for m in M:
                ax.bar(m, Gz[l][lmax+1+m], yerr=Gz_err[l][lmax+1+m], align='center', alpha=0.5, ecolor='black', capsize=5, color='tab:blue', width=w)

            ax.grid()
            ax.set_xticks(M)
            ax.set_xlabel(r"$m$", size=15)
            ax.set_ylabel(r"Field non-uniformity $\left< \left(B_z - \left<B_z\right>\right)^2 \right>$ (pT)", size=15)
            
import numpy as np
import os.path


def factorial(n):
    f=1
    for k in range(1, n+1):
        f *= k
    return f

def binomialCoef(n, k):
    return factorial(n) / (factorial(k) * factorial(n-k))

def numericalIntegral(x, y):
    n = min(len(x), len(y))
    dx = x[1:n]-x[:n-1]
    S = sum(y[:n-1]*dx)
    return S

def fourierCoef(ab, n, x, y, T):
    """
    takes: ab (string, 'a' or 'b'), y data (array), n order of coef (int), T period (float)
    returns: Fourier coef (float)
    """        
    if ab=='a':
        C = (2/T)*numericalIntegral(x, y * np.cos(2*np.pi*n*x/T))
    elif ab=='b':
        C = (2/T)*numericalIntegral(x, y * np.sin(2*np.pi*n*x/T))
    if n==0:
        C *= 0.5
    return C
        
def Gres(run, rlim=41, zlim=16):
    run.getRingCoefs()
    cycles = run.cycles['rings']
    nc = len(cycles)
    n = 0
    S = 0
    V = 0
    rings_in_volume = []
    for cyc in cycles:
        cycle = cycles[cyc]
        rho, z = cycle.rho[0], cycle.z[0]
        if rho<rlim and z>-zlim and z<zlim:
            print(round(rho), round(z))
            a0_rho_fit = run.C_model[cyc][0][0]
            a0_rho = run.C_data[cyc][0][0]
            a0_rho_res = a0_rho - a0_rho_fit
            print("{:.2f} {:.2f} {:.2f}".format(a0_rho_fit, a0_rho, a0_rho_res))
            S += rho**2 * a0_rho_res
            V += rho
            n += 1
            rings_in_volume.append(cyc)
    return S, S/V, S/V * 4/Geometry.R**2

def generateCS(run):
    cs_positions = []
    cs_fields = []
    for rho in [16, 32]:
        for z in [-32, -16, 16, 32]:
            cycle = run.cycles['rings'][rho, z]
            for i in range(0, 350, 25):
                r_cart = Tools.cartesian(cycle.rho[i], cycle.phi[i], cycle.z[i], degrees=True)
                #b_cart = Tools.cartesianField(cycle.phi[i], cycle.Brho[i], cycle.Bphi[i], cycle.Bz[i], degrees=True)
                Bm_rho, Bm_phi, Bm_z = cycle.B_model           
                b_cart = Tools.cartesianField(cycle.phi[i], Bm_rho[i], Bm_phi[i], Bm_z[i], degrees=True)
                cs_positions.append(r_cart)
                cs_fields.append(b_cart)
    cs_positions = np.array(cs_positions)
    cs_fields = np.array(cs_fields)
    cs_fields_err = 1*np.ones_like(cs_fields)  ## cs measurement error equal to RMSE of 1pT (error on poly fit)
    return cs_positions, cs_fields, cs_fields_err

def generateCS_dipole(run, rdip, mdip):
    cs_positions = []
    cs_fields = []
    for rho in [16, 32]:
        for z in [-32, -16, 16, 32]:
            cycle = run.cycles['rings'][rho, z]
            for i in range(0, 350, 25):
                r_cart = np.array(Tools.cartesian(cycle.rho[i], cycle.phi[i], cycle.z[i], degrees=True))
                #b_cart = Tools.cartesianField(cycle.phi[i], cycle.Brho[i], cycle.Bphi[i], cycle.Bz[i], degrees=True)
#                 Bm_rho, Bm_phi, Bm_z = cycle.B_model
#                 print(r_cart, rdip, mdip)
                b_cart = Dipole.field(r_cart, rdip, mdip)
#                 b_cart = Tools.cartesianField(cycle.phi[i], Bm_rho[i], Bm_phi[i], Bm_z[i], degrees=True)
                cs_positions.append(r_cart)
                cs_fields.append(b_cart)
    cs_positions = np.array(cs_positions)
    cs_fields = np.array(cs_fields)
    cs_fields_err = 1*np.ones_like(cs_fields)  ## cs measurement error equal to RMSE of 1pT (error on poly fit)
    return cs_positions, cs_fields, cs_fields_err

def numFedm(run, rlim=41, zlim=16):
    """
    takes: run and chamber limits
    returns: effective gradient (4/R**2) * <rho Brho>
    """
    cycles = run.cycles['rings']
    nc = len(cycles)
    n = 0
    S = 0
    V = 0
    rings_in_volume = []
    for cyc in cycles:
        cycle = cycles[cyc]
        rho0, z0 = cycle.rho[0], cycle.z[0]
        if rho0<rlim and z0>-zlim and z0<zlim:
            print(round(rho0), round(z0))
            S += np.sum(cycle.rho**2 * cycle.Brho)
            V += np.sum(cycle.rho)
            n += 1
            rings_in_volume.append(cyc)
    print(S, V)
    effG = (S/V )* (4/Geometry.R**2)
    return effG

def numFedm_dipole(run, rlim=41, zlim=16, rdip=[], mdip=[]):
    """
    takes: run and chamber limits
    returns: effective gradient (4/R**2) * <rho Brho>
    """
    if len(rdip)<1:
        rdip = np.array([35, 85, -41])
    if len(mdip)<1:
        mdip = 4000*1e4*np.array([1, 1, 1])/np.sqrt(3)
    cycles = run.cycles['rings']
    nc = len(cycles)
#     n = 0
    S = 0
    S_cart = 0
    V = 0
    V_cart = 0
    rings_in_volume = []
    for cyc in cycles:
        cycle = cycles[cyc]
        rho0, z0 = cycle.rho[0], cycle.z[0]
        if rho0<rlim and z0>-zlim and z0<zlim:
            Rcart = Tools.cartesian(cycle.rho, cycle.phi, cycle.z, degrees=True)
#             print(round(rho0), round(z0))
            Bdip = Dipole.fields(Rcart, rdip, mdip)
            Bdip_cyl = Tools.cylindricalField(cycle.phi, Bdip[0], Bdip[1], Bdip[2], degrees=True)
            S += np.sum(cycle.rho**2 * Bdip_cyl[0])
            V += np.sum(cycle.rho)
            S_cart += np.sum(Rcart[0]*Bdip[0] + Rcart[1]*Bdip[1])
            V_cart += len(cycle.phi)
            rings_in_volume.append(cyc)
    print(S/V)
    print(S_cart/V_cart)
#     print(rings_in_volume)
    effG = (S/V )* (-4/Geometry.R**2)
    return effG

def fedm_dipole(zd, m, R=40, H=12):
    """
    valid for a vertical dipole m e_z of strength m
    calculates <rho Brho> for a dipole field produced by this dipole
    and returns effective gradient <rho Brho> / (-R**2/4)
    """
    C = 1.256637062*10**(-1)/(4*np.pi) # pT.cm/nA
    p = C*m
    v = 2*p/(R**2*H) * ( 2*H + (R**2+2*zd**2)/(R**2+zd**2)**0.5 - (R**2+2*(zd+H)**2)/(R**2+(zd+H)**2)**0.5 )
    print(v)
    effG = v * (-4/Geometry.R**2)
    return effG

def runPlot(run, cycs=[], dphi=10, iB=2, cmap='magma', zratio=6, elev=15, azim=30, clim=[]):

    if run.runtype==None:
        cycles = run.cycles['rings']
    else:
        cycles = run.cycles
    phi = cycles[0, 0].phi
    if len(cycs)<1:
        cycs = cycles
    
    i=0
    while phi[i]<1:
        i+=1
    imin=i
    while phi[i]<359:
        i+=1
    imax=i

    Rho, Phi, Z = np.array([ [ [cycles[cyc].R[iR][iphi] for iphi in range(imin+dphi, imax-dphi, 2*dphi)] for cyc in cycs] for iR in range(3) ])
    Brho, Bphi, Bz = np.array([ [ [np.mean(cycles[cyc].B[iB][iphi-dphi:iphi+dphi]) for iphi in range(imin+dphi, imax-dphi, 2*dphi)] for cyc in cycs] for iB in range(3) ])

    X = Rho*np.cos(Phi*np.pi/180)
    Y = Rho*np.sin(Phi*np.pi/180)

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
    if len(clim)>1:
        p.set_clim(clim[0], clim[1])
    #plt.show()
    return X, Y, Z, Rho, Phi, Z

### number of G_lm up to order l (general case)
def card(l):
    return (l+1)*(l+3)

### number of G_lm up to order l for the Bz probe
def cardz(l):
    return (l+1)**2

class Geometry:
    
    R = 40
    H = 12
    Hp = 18
    
    ### coefs a_n = <rho B_rho> / (-R**2/4) in double chamber (odd) or single chamber (even)
    a1 = 1
    a2 = Hp # used to be Hp/2
    a3 = -R**2/2 + (H**2 + 3*Hp**2)/4
    a4 = (Hp**3 + Hp*H**2)/2 - R**2*Hp # used to be Hp/2 * ( -R**2 + (H**2 + Hp**2)/2 )
    a5 = (5/8) * ( 8*R**4 - 2*R**2*(H**2 + 3*Hp**2)/3 + Hp**2*H**2 + (5*Hp**4 + H**4)/10 )
    a6 = 0
    a7 = 0
    
    ### L^n coeffs from (<B_z>_TOP - <B_z>_BOT)/Hp for odd modes
    L1 = 1
    L3 = 32.9**2
    L5 = 32.7**4
    L7 = 32.5**6
    expL7 = -(1/16)*( (1/4)*(Hp**6 + 7*Hp**4*H**2 + 7*Hp**2*H**4 + H**6)
                     - (21/12)*R**2*(3*Hp**4 + 10*Hp**2*H**2 + 3*H**4)
                     + (70/3)*R**4*(Hp**2 + H**2)
                     - (35/4)*R**6)
    
    ### coefs C_n = 1 / (1 + a_n/L^n) for odd modes 
    C1 = 1
    C3 = 4*L3/(R**2 + 2*Hp**2)
    C5 = 48*L5/(15*R**4 + 10*R**2*(3*Hp**2-H**2) - 4*Hp**2*(3*Hp**2+5*H**2))
    R7 = (7/128)*R**2*( -(1/14)*(H**6 + 21*H**4 + 35*H**2*Hp**4 + 7*Hp**6) 
                       + R**2*(H**4 + 10*H**2*Hp**2 + 5*Hp**4) 
                       - (5*R**4/2)*(H**2 + 3*Hp**2) 
                       + R**6 )
    C7 = 1/(1 - (4*R7/(L7*R**2)))
    
    ### even-l (L_2n^2n-1 / c_2n)
    LC2 = Hp
    LC4 = (Hp**3 + Hp*H**2)/2 - R**2*Hp
    
    ### norm ceof
    Leff = np.array([0, L1/C1, L1/C1])
    
    def volAvg(i, j, k, R=40, H=12, Hp=18, chamber='double'):
        """
        Takes integers i, j, k, chamber geometry (R, H, H'), and chamber type ('single' or 'double')
        Returns < x^i y^j z^k > in given chamber geometry 
        """
        ### treat a = <xy> and b = <z> seperately

        if chamber=='single':
            if k%2==1:
                b = 0
            else:
                b = H**k / (2**k * (k+1))
        elif chamber=='double':
            b = 0
            for n in range(k//2 + 1):
                #print(n)
                b += binomialCoef(k, 2*n) * Hp**(k-2*n) * H**(2*n) / (2*n+1)
            b *= 1/(2**k)

        l, m = min(i, j), max(i, j)

        if l%2==1 or m%2==1:
            a = 0
        elif l==0:
            if m==0:
                a = 1
            elif m==2:
                a = R**2/4
            elif m==4:
                a = R**4/8
            elif m==6:
                a = 5*R**6/64
            elif m==8:
                a = 7*R**8/128
            else:
                print("i+j too high")
        elif l==2:
            if m==2:
                a = R**4/24
            elif m==4:
                a = R**6/64
            elif m==6:
                a = R*8/128
            else:
                print("i+j too high")
        elif l==4:
            if m==4:
                a = 3*R**8/640
            else:
                print("i+j too high")
        print(a, b)
        return a*b


    def bz(l, m):
        ### <B_z> for l, m mode
        if l==0:
            if m==0:
                b = 1
            else:
                b = 0
        elif l==1:
            if m==0:
                b = Geometry.volAvg(0, 0, 1)
            elif abs(m)==1:
                b = Geometry.volAvg(1, 0, 0)
            else:
                b = 0
        elif l==2:
            if m==-2:
                b = 2*Geometry.volAvg(1, 1, 0)
            elif m==-1:
                b = 2*Geometry.volAvg(0, 1, 1)
            elif m==0:
                b = Geometry.volAvg(0, 0, 2) - 0.5*(Geometry.volAvg(1, 0, 0)+Geometry.volAvg(0, 1, 0))
            elif m==1:
                b = 2*Geometry.volAvg(1, 0, 1)
            elif m==2:
                b = Geometry.volAvg(2, 0, 0) - Geometry.volAvg(0, 2, 0)
            else:
                b = 0
        return b

    def bz2(l, m):
        ### <B_z^2> for l, m mode
        if l==0:
            if m==0:
                b2 = 1
            else:
                b2 = 0
        elif l==1:
            if m==0:
                b2 = Geometry.volAvg(0, 0, 2)
            elif abs(m)==1:
                b2 = Geometry.volAvg(2, 0, 0)
            else:
                b2 = 0
        elif l==2:
            if m==0:
                b2 = Geometry.volAvg(0, 0, 4) - 2*Geometry.volAvg(2, 0, 2) + 0.5*(Geometry.volAvg(4, 0, 0) + Geometry.volAvg(2, 2, 0))
                
            elif abs(m)==1:
                b2 = 4*Geometry.volAvg(0, 2, 2)
            elif abs(m)==2:
                b2 = 2*Geometry.volAvg(4, 0, 0) - 2*Geometry.volAvg(2, 2, 0)
            else:
                b2 = 0
        return b2
    
    def nonUniformity(l, m):
        return Geometry.bz2(l, m) - Geometry.bz(l, m)**2

    def dratio1(g1, g3):
        return g1 / (a3*g3)

    def dratio2(g2, g4):
        return g2 / (a4*g4)

    def dratio2(g3, g5):
        return g3 / (a5*g5)
    

class Tools:
    
    def cylindrical(x, y, z):
        rho, phi, z = (x**2 + y**2)**0.5, np.arctan2(y, x), z
        return rho, phi, z
    
    def cylindricalField(phi, Bx, By, Bz, degrees=True):
        if degrees==True:
            phi = phi*np.pi/180
        Brho = Bx*np.cos(phi) + By*np.sin(phi)
        Bphi = -Bx*np.sin(phi) + By*np.cos(phi)
        return Brho, Bphi, Bz

    def cartesian(rho, phi, z, degrees=True):
        if degrees==True:
            phi = np.pi*phi/180
        x, y, z = rho*np.cos(phi), rho*np.sin(phi), z
        return x, y, z
    
    def cartesianField(phi, Brho, Bphi, Bz, degrees=True):
        if degrees==True:
            phi = phi*np.pi/180
        Bx = Brho*np.cos(phi) - Bphi*np.sin(phi)
        By = Brho*np.sin(phi) + Bphi*np.cos(phi)
        return Bx, By, Bz
    
    def gradScan(x):
        imax=[]
        imin=[]
        if (x[1]-x[0])>0:
            imin.append(0)
        elif (x[1]-x[0])<0:
            imax.append(0)
        for i in range(len(x)-2):
            if np.sign(x[i+1]-x[i])!=np.sign(x[i+2]-x[i+1]):
                if (x[i+1]-x[i])>0:
                    imax.append(i+1)
                elif (x[i+1]-x[i])<0:
                    imin.append(i+1)
        return np.array(imin), np.array(imax)
    
    def slopeChange(x):
        I = []
        for i in range(len(x)-2):
            if np.sign(x[i+1]-x[i])!=np.sign(x[i+2]-x[i+1]):
                I.append(i+1)
        return np.array(I)
    
    def chi2stat(ydata, ymodel, yerr):
        c = np.sum(((ydata-ymodel)/yerr)**2)
        return c
    
    def chi2test(ydata, ymodel):
        return np.sum((ydata-ymodel)**2/ymodel)
    
    def corr_coef(x, y):
        return np.mean( (x-np.mean(x))*(y-np.mean(y)) ) / ( np.std(x)*np.std(y) )
    
    def zero_crossing(Y, up=True):
        i=0
        if up==True:
            while Y[i]>0:
                i+=1
        else:
            while Y[i]<0:
                i+=1
        return i
    
    def smooth_extremum(x, mode='max', Dx=None, mrange=None):
        """
        mode = 'min' or 'max'
        Dx = rejection criteria: rejects x[i] if x[i]-x[i-1] < Dx
        mrange = tuple [min, max]
        """
        n = len(x)
        if Dx==None:
            Dx = np.std(x)
        if mrange==None:
            mrange = [min(x), max(x)]
        # initialize extremum
        i = 0
        while abs(x[i+1]-x[i])>=Dx:
            i += 1
        m = x[i]
        argm = i
        # search for extremum
        for i in range(1, n):
            if mode=='max':
                if x[i]>m and abs(x[i]-x[i-1])<Dx and x[i]<=mrange[1]:
                    m = x[i]
                    argm = i
            elif mode=='min':
                if x[i]<m and abs(x[i]-x[i-1])<Dx and x[i]>=mrange[0]:
                    m = x[i]
                    argm = i
        return argm, m
        
    def m_to_m(m, iB, lmax, mode='relative'):
        """
        mode = 'relative' or 'positive'
        'relative' transforms positive m indices (0, 2*lmax+3)
        to relative m indices (-lmax-1, 0, -1) U (0, lmax+2)
        """
        if iB not in [0, 1, 2]:
            print("error: iB must be 0, 1, or 2 (rho, phi, or z)")
            return None
        if mode=='positive':
            if iB==0 or iB==2:
                if m>=0:
                    m0 = m
                else:
                    m0 = lmax+1-m
            elif iB==1:
                if m>=0:
                    m0 = lmax+1+m
                else:
                    m0 = -m
        elif mode=='relative':
            if iB==0 or iB==2:
                if m<=(lmax+1):
                    m0 = m
                else:
                    m0 = lmax+1-m
            elif iB==1:
                if m<=(lmax+1):
                    m0 = -m
                else:
                    m0 = m-lmax-1
        else:
            print("mode = 'positive' or 'relative'")
            return None
        return m0         
    
        
class Functions:
    
    def polynomial(t, a):
        """
        a = array of polynomial coefficients
        """
        f=0
        for n in range(len(a)):
            f += a[n]*t**n
        return f
    
    def fourier(t, a, b, T=2*np.pi):
        """
        takes: t scalar or array, fourier coefs a n-dim array /!\, b (n-1)-dim array /!\
        returns: fourier series
        """
        f = 0
        for n in range(len(a)):
            f += a[n]*np.cos(n*t*2*np.pi/T) 
            if n>0:
                f += b[n-1]*np.sin(n*t*2*np.pi/T)
        return f
    
    
class Fits:
    
    def polynomial(y0, t, N, sigma=None):
        """
        takes: y0 array (values to be fitted), t array (fit points), N order of the fit, variance of the measure points
        returns: N size array of polynomial coefficients (fit parameters) [a_0, a_1, ..., a_N], with their associated errors
        """

        if np.any(sigma==None):
            sigma = np.ones_like(y0)
        
        u = np.zeros(N+1)
        U = np.zeros((N+1, N+1))
        I = len(t)

        for i in range(I):
            for n in range(N+1):
                u[n] += y0[i]*t[i]**n / sigma[i]**2
                for m in range(N+1):
                    U[n][m] += t[i]**(n+m) / sigma[i]**2

        Uinv = np.linalg.inv(U)

        par = np.dot(Uinv, u)
        err = np.array([Uinv[i][i]**0.5 for i in range(N+1)])

        return par, err
    
    def fourier(y0, t, N, T=2*np.pi, sigma=None):
        """
        takes: y0 array (values to be fitted), t array (fit points), N order of the fit, variance of the measure points
        returns: 2*N+1 size array of fit parameters [a_0, a_1, ..., a_N, b_1, ..., b_N], along with errors on the parameters
        """

        if np.any(sigma==None):
            sigma = np.ones_like(y0)
            
        u = np.zeros(2*N+1)
        U = np.zeros((2*N+1, 2*N+1))
        I = len(t)

        for i in range(I):
            for n in range(2*N+1):

                if n<=N:
                    u[n] += y0[i]*np.cos(n*t[i]*2*np.pi/T) / sigma[i]**2
                    for m in range(2*N+1):
                        if m<=N:
                            U[n][m] += np.cos(m*t[i]*2*np.pi/T)*np.cos(n*t[i]*2*np.pi/T) / sigma[i]**2
                        else:
                            U[n][m] += np.sin((m-N)*t[i]*2*np.pi/T)*np.cos(n*t[i]*2*np.pi/T) / sigma[i]**2

                else:
                    u[n] += y0[i]*np.sin((n-N)*t[i]*2*np.pi/T) / sigma[i]**2
                    for m in range(2*N+1):
                        if m<N:
                            U[n][m] += np.cos(m*t[i]*2*np.pi/T)*np.sin((n-N)*t[i]*2*np.pi/T) / sigma[i]**2
                        else:
                            U[n][m] += np.sin((m-N)*t[i]*2*np.pi/T)*np.sin((n-N)*t[i]*2*np.pi/T) / sigma[i]**2

        Uinv = np.linalg.inv(U)

        par = np.dot(Uinv, u)
        err = np.array([Uinv[i][i]**0.5 for i in range(2*N+1)])

        return par, err
    
    def GlmPoly(iB, a, a_err, ring, m, lmin, lmax):
        """
        takes: iB field proj, a array of fourier coefs of order m to be fitted, ring array of tuples (rho_ring, z_ring), m order of fourier coef and Glms, L order of the fit, sigma std of the measured values
        returns: l-polynomial coefficient of Pi_lm coefficients (fixed m)
        """

        u = np.zeros(lmax+1-lmin)
        U = np.zeros((lmax+1-lmin, lmax+1-lmin))
        I = len(ring)

        for i in range(I):
            rho = ring[i][0]
            z =  ring[i][1]
            for l1 in range(lmin, lmax+1):
                u[l1-lmin] += Harmonic.PiI(iB, rho, z, l1, m) * a[i] / a_err[i]**2
                for l2 in range(lmin, lmax+1):
                    U[l1-lmin][l2-lmin] += Harmonic.PiI(iB, rho, z, l1, m) * Harmonic.PiI(iB, rho, z, l2, m) / a_err[i]**2

        Uinv = np.linalg.inv(U)

        g = np.dot(Uinv, u)
        g_err = np.array([Uinv[i][i]**0.5 for i in range(lmax+1-lmin)])

        return g, g_err
    
    def GlmPoly_cartesian(iB, B, B_err, R, lmax, unpack=True):
        """
        takes: iB cartesian field proj, cartesian field and its error, cartesian coordinates of where field is measured, max l-order of fit
        returns: G_lm coefficients
        """

        def unpack_g(g, lmax, iB=2):
            ### only works for iB=2 (z probe)
            G = np.zeros((lmax+1, 2*lmax+1))
            for l in range(lmax+1):
                for m in range(-l, l+1):
                    n = cardz(l-1) + m + l
                    m0 = m + lmax
                    G[l][m0] = g[n]
            return G
        
        ### /!\ only working case /!\
        if iB!=2:
            print("only works for Bz probe")
            return None

        L = cardz(lmax)
        u = np.zeros(L)
        U = np.zeros((L, L))
        I = len(B)

        for i in range(I):
            x, y, z = R[i]
            b = B[i][iB]
            b_err = B_err[i][iB]
            for l1 in range(lmax+1):
                for m1 in range(-l1, l1+1):
    #             for m1 in range(-l1-1, l1+2):
                    n1 = cardz(l1-1) + m1 + l1 # +1
                    u[n1] += Harmonic.PiI_cart(iB, x, y, z, l1, m1) * b / b_err**2
                    for l2 in range(lmax+1):
                        for m2 in range(-l2, l2+1):
    #                     for m2 in range(-l2-1, l2+2):
                            n2 = cardz(l2-1) + m2 + l2 # +1
                            U[n1][n2] += Harmonic.PiI_cart(iB, x, y, z, l1, m1) * Harmonic.PiI_cart(iB, x, y, z, l2, m2) / b_err**2 

        Uinv = np.linalg.inv(U)

        g = np.dot(Uinv, u)
        g_err = np.array([Uinv[i][i]**0.5 for i in range(L)])

        if unpack==False:
            return g, g_err
        else:
            g = unpack_g(g, lmax, iB)
            g_err = unpack_g(g_err, lmax, iB)
            return g, g_err

        
class Dipole:
    
    def field(r, rd, m):
        "r, rd, m shape(3) in cartesian coordinates"
        C = 1.256637062*10**(-1)/(4*np.pi) # pT.cm/nA
        N = np.dot(r-rd,r-rd)**(1/2)
        u = (r-rd)/N 
        Bd = C*(3*np.dot(m,u)*u - m)/N**3
        return Bd
    
    def fields(R, rd, md):
        "R array shape (3, n) in cartesian coordinates"
        RT = np.transpose(R)
        BT = []
        for r in RT:
            BT.append(Dipole.field(r, rd, md))
        B = np.transpose(BT)
        return B
    
    def multiFields(R, Rd, Md):
        "R, Rd, M arrays shape (3, n) in cartesian coordinates"
        RdT, MdT = np.transpose(Rd), np.transpose(Md)
        ndip = len(RdT)
        for i in range(ndip):
            if i==0:
                B = Dipole.fields(R, RdT[i], MdT[i])
            else:
                B += Dipole.fields(R, RdT[i], MdT[i])
        return B 
    
    def matrix(r, rd, inv=False):
        "r and rd shape (3) in cartesian coordinates"
        C = 1.256637062*10**(-1)/(4*np.pi) # pT.cm/nA
        x, y, z = r-rd
        N = np.dot(r-rd,r-rd)**(1/2)
        A = -N**(-3)*np.identity(3) + 3*N**(-5)*np.outer([x, y, z], [x, y, z])
        Ainv = -N**3*np.identity(3) + 1.5*N*np.outer([x, y, z], [x, y, z])
        if inv==False:
            return C*A
        else:
            return (1/C)*Ainv
        
    def multiMatrix(r, Rd):
        "r of shape (3), Rd of shape (3, n) in cartesian coordinates"
        RdT = np.transpose(Rd)
        F = np.zeros((3, 3))
        for i in range(len(RdT)):
            F += Dipole.matrix(r, RdT[i])
        return F
    
    def fit(R, B, rd, sigma=1):
        "R, B of shape (3, n), rd of shape(3) in cartesian coordinates"
        RT = np.transpose(R)
        BT = np.transpose(B)
        F2sum = np.zeros((3,3))
        bFsum = np.zeros(3)
        for i in range(len(RT)):
            F = Dipole.matrix(RT[i], rd)
            F2sum += np.dot(F, F)
            bFsum += np.dot(F, BT[i])
        F2suminv = np.linalg.inv(F2sum)
        mpar = np.dot(F2suminv, bFsum)
        merr = sigma*np.array([F2suminv[i][i]**0.5 for i in range(3)])
        return mpar, merr
    
    def multiFit(R, B, Rd, sigma=1):
        "R, B, Rd of shape (3, n) in cartesian coordinates"
        RT = np.transpose(R)
        BT = np.transpose(B)
        F2sum = np.zeros((3,3))
        bFsum = np.zeros(3)
        for i in range(len(RT)):
            F = Dipole.multiMatrix(RT[i], Rd)
            F2sum += np.dot(F, F)
            bFsum += np.dot(F, BT[i])
        F2suminv = np.linalg.inv(F2sum)
        mpar = np.dot(F2suminv, bFsum)
        merr = sigma*np.array([F2suminv[i][i]**0.5 for i in range(3)])
        return mpar, merr
    
    
class Harmonic:  
    
    def Bharmonic_cart(r, G):
        lmax = len(G)-1
        x, y, z = r
        Bx, By, Bz = 0
        for l in range(lmax+1):
            for m in range(-l-1, l+2):
                Bx += G[l][lmax+1+m] * Harmonic.PiI_cart(0, x, y, z, l, m)
                By += G[l][lmax+1+m] * Harmonic.PiI_cart(1, x, y, z, l, m)
                Bz += G[l][lmax+1+m] * Harmonic.PiI_cart(2, x, y, z, l, m)
        B = np.array([Bx, By, Bz])
        return B
    
    def Bharmonic_cyl(r, G):
        lmax = len(G)-1
        rho, phi, z = r
        phi = phi * np.pi / 180 # in radians for sine functions
        Brho, Bphi, Bz = 0, 0, 0
        for l in range(lmax+1):
            for m in range(-l-1, l+2):
    #         for m in range(-l, l+1):
    #             print(l, m)
                if m>=0:
                    Brho += G[l][lmax+1+m] * np.cos(m*phi) * Harmonic.PiI(0, rho, z, l, m)
                    Bphi += G[l][lmax+1+m] * np.sin(m*phi) * Harmonic.PiI(1, rho, z, l, m)
                    Bz += G[l][lmax+1+m] * np.cos(m*phi) * Harmonic.PiI(2, rho, z, l, m)
                elif m<0:
                    Brho += G[l][lmax+1+m] * np.sin(-m*phi) * Harmonic.PiI(0, rho, z, l, m)
                    Bphi += G[l][lmax+1+m] * np.cos(m*phi) * Harmonic.PiI(1, rho, z, l, m)
                    Bz += G[l][lmax+1+m] * np.sin(-m*phi) * Harmonic.PiI(2, rho, z, l, m)
        B = np.array([Brho, Bphi, Bz])
        return B
    
    def PiI(iB, rho, z, l, m):
        if iB==0:
            return Harmonic.PiRho(rho, z, l, m)
        elif iB==1:
            return Harmonic.PiPhi(rho, z, l, m)
        elif iB==2:
            return Harmonic.PiZ(rho, z, l, m)
    
    def PiRho(rho, z, l, m):
        
        if l==0:
            
            if m==-1:
                P = 1
            elif m==0:
                P = 0
            elif m==1:
                P = 1
            else:
                print("error: m must match condition -l-1 < m < l+1")
                return None
                
        elif l==1:
            
            if m==-2:
                P = rho
            elif m==-1:
                P = z
            elif m==0:
                P = -0.5*rho
            elif m==1:
                P = z
            elif m==2:
                P = rho
            else:
                print("error: m must match condition -l-1 < m < l+1")
                return None
                
        elif l==2:
            
            if m==-3:
                P = rho**2
            elif m==-2:
                P = 2*rho*z
            elif m==-1:
                P = 0.25*(4*z**2 - 3*rho**2)
            elif m==0:
                P = -rho*z
            elif m==1:
                P = 0.25*(4*z**2 - 3*rho**2)
            elif m==2:
                P = 2*rho*z
            elif m==3:
                P = rho**2
            else:
                print("error: m must match condition -l-1 < m < l+1")
                return None
            
        elif l==3:
            
            if m==-4:
                P = rho**3
            elif m==-3:
                P = 3*rho**2*z
            elif m==-2:
                P = rho*(3*z**2 - rho**2)
            elif m==-1:
                P = 0.25*z*(4*z**2 - 9*rho**2)
            elif m==0:
                P = (3/8)*rho*(rho**2 - 4*z**2)
            elif m==1:
                P = 0.25*z*(4*z**2 - 9*rho**2)
            elif m==2:
                P = rho*(3*z**2 - rho**2)
            elif m==3:
                P = 3*rho**2*z
            elif m==4:
                P = rho**3
            else:
                print("error: m must match condition -l-1 < m < l+1")
                return None
            
        elif l==4:
            
            if m==-5:
                P = rho**4
            elif m==-4:
                P = 4*rho**3*z
            elif m==-3:
                P = 0.25*(24*rho**2*z**2 - 5*rho**4)
            elif m==-2:
                P = 4*(rho*z**3 - rho**3*z)
            elif m==-1:
                P = 0.125*(8*z**4 - 36*rho**2*z**2 + 5*rho**4)
            elif m==0:
                P = 0.5*(3*rho**3*z - 4*rho*z**3)
            elif m==1:
                P = 0.125*(8*z**4 - 36*rho**2*z**2 + 5*rho**4)
            elif m==2:
                P = 4*(rho*z**3 - rho**3*z)
            elif m==3:
                P = 0.25*(24*rho**2*z**2 - 5*rho**4)
            elif m==4:
                P = 4*rho**3*z
            elif m==5:
                P = rho**4
            else:
                print("error: m must match condition -l-1 < m < l+1")
                return None
        
        elif l==5:
            
            if m==-6:
                P = rho**5
            elif m==-5:
                P = 5*rho**4*z
            elif m==-4:
                P = 0.5*(20*rho**3*z**2 - 3*rho**5)
            elif m==-3:
                P = 1.25*(8*rho**2*z**3 - 5*rho**4*z)
            elif m==-2:
                P = (5/16)*(16*rho*z**4 - 32*rho**3*z**2 + 3*rho**5)
            elif m==-1:
                P = 0.125*(8*z**5 - 60*rho**2*z**3 + 25*rho**4*z)
            elif m==0:
                P = (5/16)*(-8*rho*z**4 + 12*rho**3*z**2 - rho**5)
            elif m==1:
                P = 0.125*(8*z**5 - 60*rho**2*z**3 + 25*rho**4*z)
            elif m==2:
                P = (5/16)*(16*rho*z**4 - 32*rho**3*z**2 + 3*rho**5)
            elif m==3:
                P = 1.25*(8*rho**2*z**3 - 5*rho**4*z)
            elif m==4:
                P = 0.5*(20*rho**3*z**2 - 3*rho**5)
            elif m==5:
                P = 5*rho**4*z
            elif m==6:
                P = rho**5
            else:
                print("error: m must match condition -l-1 < m < l+1")
                return None
        
        elif l==6:
            
            if m==-7:
                P = rho**6
            elif m==-6:
                P = 6*rho**5*z
            elif m==-5:
                P = 0.25*rho**4*(60*z**2 - 7*rho**2)
            elif m==-4:
                P = rho**3*z*(20*z**2 - 9*rho**2)
            elif m==-3:
                P = (3/16)*rho**2*(80*z**4 - 100*rho**2*z**2 + 7*rho**4)
            elif m==-2:
                P = 0.125*rho*z*(48*z**4 - 160*rho**2*z**2 + 45*rho**4)
            elif m==-1:
                P = (1/64)*(64*z**6 - 720*rho**2*z**2 + 600*rho**4*z**2 - 35*rho**6)
            elif m==0:
                P = (3/8)*rho*(-8*z**5 + 20*rho**2*z**3 - 5*rho**4*z)
            elif m==1:
                P = (1/64)*(64*z**6 - 720*rho**2*z**2 + 600*rho**4*z**2 - 35*rho**6)
            elif m==2:
                P = 0.125*rho*z*(48*z**4 - 160*rho**2*z**2 + 45*rho**4)
            elif m==3:
                P = (3/16)*rho**2*(80*z**4 - 100*rho**2*z**2 + 7*rho**4)
            elif m==4:
                P = rho**3*z*(20*z**2 - 9*rho**2)
            elif m==5:
                P = 0.25*rho**4*(60*z**2 - 7*rho**2)
            elif m==6:
                P = 6*rho**5*z
            elif m==7:
                P = rho**6
            else:
                print("error: m must match condition -l-1 < m < l+1")
                return None
            
        elif l==7:
            
            if m==-8: 
                P = pow(rho, 7)
            elif m==-7: 
                P = 7*pow(rho, 6)*z
            elif m==-6: 
                P = 21*pow(rho, 5)*pow(z, 2) - 2*pow(rho, 7)
            elif m==-5: 
                P = (7./4) * (20*pow(rho, 4)*pow(z, 3) - 7*pow(rho, 6)*z)
            elif m==-4: 
                P = (7./4) * (20*pow(rho, 3)*pow(z, 4) - 18*pow(rho, 5)*pow(z, 2) + pow(rho, 7))
            elif m==-3: 
                P = (7./16) * (48*pow(rho, 2)*pow(z, 5) - 100*pow(rho, 4)*pow(z, 3) + 21*pow(rho, 6)*z)
            elif m==-2: 
                P = (7./16) * (16*rho*pow(z, 6) - 80*pow(rho, 3)*pow(z, 4) + 45*pow(rho, 5)*pow(z, 2) - 2*pow(rho, 7))
            elif m==-1: 
                P = (1./64) * (64*pow(z, 7) - 1008*pow(rho, 2)*pow(z, 5) + 1400*pow(rho, 4)*pow(z, 3) - 245*pow(rho, 6)*z)
            elif m==0: 
                P = (7./128) * (-64*rho*pow(z, 6) + 240*pow(rho, 3)*pow(z, 4) - 120*pow(rho, 5)*pow(z, 2) + 5*pow(rho, 7))
            elif m==1: 
                P = (1./64) * (64*pow(z, 7) - 1008*pow(rho, 2)*pow(z, 5) + 1400*pow(rho, 4)*pow(z, 3) - 245*pow(rho, 6)*z)
            elif m==2: 
                P = (7./16) * (16*rho*pow(z, 6) - 80*pow(rho, 3)*pow(z, 4) + 45*pow(rho, 5)*pow(z, 2) - 2*pow(rho, 7))
            elif m==3: 
                P = (7./16) * (48*pow(rho, 2)*pow(z, 5) - 100*pow(rho, 4)*pow(z, 3) + 21*pow(rho, 6)*z)
            elif m==4: 
                P = (7./4) * (20*pow(rho, 3)*pow(z, 4) - 18*pow(rho, 5)*pow(z, 2) + pow(rho, 7))
            elif m==5: 
                P = (7./4) * (20*pow(rho, 4)*pow(z, 3) - 7*pow(rho, 6)*z)
            elif m==6: 
                P = 21*pow(rho, 5)*pow(z, 2) - 2*pow(rho, 7)
            elif m==7: 
                P = 7*pow(rho, 6)*z
            elif m==8: 
                P = pow(rho, 7)
            else:
                print("error: m must match condition -l-1 < m < l+1")
                return None

        else:
            print("l too high")
            return None
        
        return P
    
    def PiPhi(rho, z, l, m):
        
        if l==0:
            
            if m==-1:
                P = 1
            elif m==0:
                P = 0
            elif m==1:
                P = -1
        elif l==1:
            
            if m==-2:
                P = rho
            elif m==-1:
                P = z
            elif m==0:
                P = 0
            elif m==1:
                P = -z
            elif m==2:
                P = -rho
                
        elif l==2:
            
            if m==-3:
                P = rho**2
            elif m==-2:
                P = 2*rho*z
            elif m==-1:
                P = 0.25*(4*z**2 - rho**2)
            elif m==0:
                P = 0
            elif m==1:
                P = -0.25*(4*z**2 - rho**2)
            elif m==2:
                P = -2*rho*z
            elif m==3:
                P = -rho**2
            else:
                print("error: m must match condition -l-1 < m < l+1")
                return None
            
        elif l==3:
            
            if m==-4:
                P = rho**3
            elif m==-3:
                P = 3*rho**2*z
            elif m==-2:
                P = 0.5*rho*(6*z**2 - rho**2)
            elif m==-1:
                P = 0.25*z*(4*z**2 - 3*rho**2)
            elif m==0:
                P = 0
            elif m==1:
                P = -0.25*z*(4*z**2 - 3*rho**2)
            elif m==2:
                P = -0.5*rho*(6*z**2 - rho**2)
            elif m==3:
                P = -3*rho**2*z
            elif m==4:
                P = -rho**3
            else:
                print("error: m must match condition -l-1 < m < l+1")
                return None
            
        elif l==4:
            
            if m==-5:
                P = rho**4
            elif m==-4:
                P = 4*rho**3*z
            elif m==-3:
                P = 0.75*(8*rho**2*z**2 - rho**4)
            elif m==-2:
                P = 2*(2*rho*z**3 - rho**3*z)
            elif m==-1:
                P = 0.125*(8*z**4 - 12*rho**2*z**2 + rho**4)
            elif m==0:
                P = 0
            elif m==1:
                P = -0.125*(8*z**4 - 12*rho**2*z**2 + rho**4)
            elif m==2:
                P = -2*(2*rho*z**3 - rho**3*z)
            elif m==3:
                P = -0.75*(8*rho**2*z**2 - rho**4)
            elif m==4:
                P = -4*rho**3*z
            elif m==5:
                P = -rho**4
            else:
                print("error: m must match condition -l-1 < m < l+1")
                return None
        
        elif l==5:
            
            if m==-6:
                P = rho**5
            elif m==-5:
                P = 5*rho**4*z
            elif m==-4:
                P = rho**3*(10*z**2 - rho**2)
            elif m==-3:
                P = 1.25*(8*rho**2*z**3 - 3*rho**4*z)
            elif m==-2:
                P = (5/16)*(16*rho*z**4 - 16*rho**3*z**2 + rho**5)
            elif m==-1:
                P = 0.125*(8*z**5 - 20*rho**2*z**3 + 5*rho**4*z)
            elif m==0:
                P = 0
            elif m==1:
                P = -0.125*(8*z**5 - 20*rho**2*z**3 + 5*rho**4*z)
            elif m==2:
                P = -(5/16)*(16*rho*z**4 - 16*rho**3*z**2 + rho**5)
            elif m==3:
                P = -1.25*(8*rho**2*z**3 - 3*rho**4*z)
            elif m==4:
                P = -rho**3*(10*z**2 - rho**2)
            elif m==5:
                P = -5*rho**4*z
            elif m==6:
                P = -rho**5
            else:
                print("error: m must match condition -l-1 < m < l+1")
                return None
            
        elif l==6:
            
            if m==-7:
                P = rho**6
            elif m==-6:
                P = 6*rho**5*z
            elif m==-5:
                P = 5*rho**4*(12*z**2 - rho**2)
            elif m==-4:
                P = 2*rho**3*z*(10*z**2 - 3*rho**2)
            elif m==-3:
                P = (3/16)*rho**2*(80*z**4 - 60*rho**2*z**2 + 3*rho**4)
            elif m==-2:
                P = 0.125*rho*z*(48*z**4 - 80*rho**2*z**2 + 15*rho**4)
            elif m==-1:
                P = (1/64)*(64*z**6 - 240*rho**2*z**2 + 120*rho**4*z**2 - 5*rho**6)
            elif m==0:
                P = 0
            elif m==1:
                P = -(1/64)*(64*z**6 - 240*rho**2*z**2 + 120*rho**4*z**2 - 5*rho**6)
            elif m==2:
                P = -0.125*rho*z*(48*z**4 - 80*rho**2*z**2 + 15*rho**4)
            elif m==3:
                P = -(3/16)*rho**2*(80*z**4 - 60*rho**2*z**2 + 3*rho**4)
            elif m==4:
                P = -2*rho**3*z*(10*z**2 - 3*rho**2)
            elif m==5:
                P = -5*rho**4*(12*z**2 - rho**2)
            elif m==6:
                P = -6*rho**5*z
            elif m==7:
                P = -rho**6
            else:
                print("error: m must match condition -l-1 < m < l+1")
                return None
            
        elif l==7:
            
            if m==-8: 
                P = pow(rho, 7)
            elif m==-7: 
                P = 7*pow(rho, 6)*z
            elif m==-6: 
                P = (3./2) * (14*pow(rho, 5)*pow(z, 2) - pow(rho, 7))
            elif m==-5: 
                P = (35./4) * (4*pow(rho, 4)*pow(z, 3) - pow(rho, 6)*z)
            elif m==-4: 
                P = (7./8) * (40*pow(rho, 3)*pow(z, 4) - 24*pow(rho, 5)*pow(z, 2) + pow(rho, 7))
            elif m==-3: 
                P = (21./16) * (16*pow(rho, 2)*pow(z, 5) - 20*pow(rho, 4)*pow(z, 3) + 3*pow(rho, 6)*z)
            elif m==-2: 
                P = (7./32) * (32*rho*pow(z, 6) - 80*pow(rho, 3)*pow(z, 4) + 30*pow(rho, 5)*pow(z, 2) - pow(rho, 7))
            elif m==-1: 
                P = (1./64) * (64*pow(z, 7) - 336*pow(rho, 2)*pow(z, 5) + 280*pow(rho, 4)*pow(z, 3) - 35*pow(rho, 6)*z)
            elif m==0: 
                P = 0
            elif m==1: 
                P = -(1./64) * (64*pow(z, 7) - 336*pow(rho, 2)*pow(z, 5) + 280*pow(rho, 4)*pow(z, 3) - 35*pow(rho, 6)*z)
            elif m==2: 
                P = -(7./32) * (32*rho*pow(z, 6) - 80*pow(rho, 3)*pow(z, 4) + 30*pow(rho, 5)*pow(z, 2) - pow(rho, 7))
            elif m==3: 
                P = -(21./16) * (16*pow(rho, 2)*pow(z, 5) - 20*pow(rho, 4)*pow(z, 3) + 3*pow(rho, 6)*z)
            elif m==4: 
                P = -(7./8) * (40*pow(rho, 3)*pow(z, 4) - 24*pow(rho, 5)*pow(z, 2) + pow(rho, 7))
            elif m==5: 
                P = -(35./4) * (4*pow(rho, 4)*pow(z, 3) - pow(rho, 6)*z)
            elif m==6: 
                P = -(3./2) * (14*pow(rho, 5)*pow(z, 2) - pow(rho, 7))
            elif m==7: 
                P = -7*pow(rho, 6)*z
            elif m==8: 
                P = -pow(rho, 7)
            else:
                print("error: m must match condition -l-1 < m < l+1")
                return None
            
        else:
            print("l too high")
            return None
        
        return P
    
    def PiZ(rho, z, l, m):
        
        if l==0:
            
            if m==-1:
                P = 0
            elif m==0:
                P = 1
            elif m==1:
                P = 0
            else:
                print("error: m must match condition -l-1 < m < l+1")
                return None
            
        elif l==1:
            
            if m==-2:
                P = 0
            elif m==-1:
                P = rho
            elif m==0:
                P = z
            elif m==1:
                P = rho
            elif m==2:
                P = 0
            else:
                print("error: m must match condition -l-1 < m < l+1")
                return None
                
        elif l==2:
            
            if m==-3:
                P = 0
            elif m==-2:
                P = rho**2
            elif m==-1:
                P = 2*rho*z
            elif m==0:
                P = -0.5*rho**2 + z**2
            elif m==1:
                P = 2*rho*z
            elif m==2:
                P = rho**2
            elif m==3:
                P = 0
            else:
                print("error: m must match condition -l-1 < m < l+1")
                return None
        
        elif l==3:
            
            if m==-4:
                P = 0
            elif m==-3:
                P = rho**3
            elif m==-2:
                P = 3*rho**2*z
            elif m==-1:
                P = rho*(3*z**2 - 0.75*rho**2)
            elif m==0:
                P = 0.5*z*(2*z**2 - 3*rho**2)
            elif m==1:
                P = rho*(3*z**2 - 0.75*rho**2)
            elif m==2:
                P = 3*rho**2*z
            elif m==3:
                P = rho**3
            elif m==4:
                P = 0
            else:
                print("error: m must match condition -l-1 < m < l+1")
                return None
            
        elif l==4:
            
            if m==-5:
                P = 0
            elif m==-4:
                P = rho**4
            elif m==-3:
                P = 4*rho**3*z
            elif m==-2:
                P = 6*rho**2*z**2 - rho**4
            elif m==-1:
                P = 4*rho*z**3 - 3*rho**3*z
            elif m==0:
                P = 0.125*(8*z**4 - 24*rho**2*z**2 + 3*rho**4)
            elif m==1:
                P = 4*rho*z**3 - 3*rho**3*z
            elif m==2:
                P = 6*rho**2*z**2 - rho**4
            elif m==3:
                P = 4*rho**3*z
            elif m==4:
                P = rho**4
            elif m==5:
                P = 0
            else:
                print("error: m must match condition -l-1 < m < l+1")
                return None
            
        elif l==5:
            
            if m==-6:
                P = 0
            elif m==-5:
                P = rho**5
            elif m==-4:
                P = 5*rho**4*z
            elif m==-3:
                P = 1.25*(8*rho**3*z**2 - rho**5)
            elif m==-2:
                P = 5*(2*rho**2*z**3 - rho**4*z)
            elif m==-1:
                P = (5/8)*(8*rho*z**4 - 12*rho**3*z**2 + rho**5)
            elif m==0:
                P = (1/8)*(8*z**5 - 40*rho**2*z**3 + 15*rho**4*z)
            elif m==1:
                P = (5/8)*(8*rho*z**4 - 12*rho**3*z**2 + rho**5)
            elif m==2:
                P = 5*(2*rho**2*z**3 - rho**4*z)
            elif m==3:
                P = 1.25*(8*rho**3*z**2 - rho**5)
            elif m==4:
                P = 5*rho**4*z
            elif m==5:
                P = rho**5
            elif m==6:
                P = 0
            else:
                print("error: m must match condition -l-1 < m < l+1")
                return None
            
        elif l==6:
            
            if m==-7:
                P = 0
            elif m==-6:
                P = rho**6
            elif m==-5:
                P = 6*rho**5*z
            elif m==-4:
                P = 1.5*rho**4*(10*z**2 - rho**2)
            elif m==-3:
                P = 2.5*rho**3*z*(8*z**2 - 3*rho**2)
            elif m==-2:
                P = (15/16)*rho**2*(16*z**4 - 16*rho**2*z**2 + rho**4)
            elif m==-1:
                P = 0.75*rho*z*(8*z**4 - 20*rho**2*z**2 + 5*rho**4)
            elif m==0:
                P = (1/16)*(16*z**6 - 120*rho**2*z**4 + 90*rho**4*z**2 - 5*rho**6)
            elif m==1:
                P = 0.75*rho*z*(8*z**4 - 20*rho**2*z**2 + 5*rho**4)
            elif m==2:
                P = (15/16)*rho**2*(16*z**4 - 16*rho**2*z**2 + rho**4)
            elif m==3:
                P = 2.5*rho**3*z*(8*z**2 - 3*rho**2)
            elif m==4:
                P = 1.5*rho**4*(10*z**2 - rho**2)
            elif m==5:
                P = 6*rho**5*z
            elif m==6:
                P = rho**6
            elif m==7:
                P = 0
            else:
                print("error: m must match condition -l-1 < m < l+1")
                return None
            
        elif l==7:
            
            if m==-8: 
                P = 0.
            elif m==-7: 
                P = pow(rho, 7)
            elif m==-6: 
                P = 7*pow(rho, 6)*z
            elif m==-5: 
                P = (7./4) * (12*pow(rho, 5)*pow(z, 2) - pow(rho, 7))
            elif m==-4: 
                P = (7./2) * (10*pow(rho, 4)*pow(z, 3) - 3*pow(rho, 6)*z)
            elif m==-3: 
                P = (7./16) * (80*pow(rho, 3)*pow(z, 4) - 60*pow(rho, 5)*pow(z, 2) + 3*pow(rho, 7))
            elif m==-2: 
                P = (7./16) * (48*pow(rho, 2)*pow(z, 5) - 80*pow(rho, 4)*pow(z, 3) + 15*pow(rho, 6)*z)
            elif m==-1: 
                P = (7./64) * (64*rho*pow(z, 6) - 240*pow(rho, 3)*pow(z, 4) + 120*pow(rho, 5)*pow(z, 2) - 5*pow(rho, 7))
            elif m==0: 
                P = (1./16) * (16*pow(z, 7) - 168*pow(rho, 2)*pow(z, 5) + 210*pow(rho, 4)*pow(z, 3) - 35*pow(rho, 6)*z)
            elif m==1: 
                P = (7./64) * (64*rho*pow(z, 6) - 240*pow(rho, 3)*pow(z, 4) + 120*pow(rho, 5)*pow(z, 2) - 5*pow(rho, 7))
            elif m==2: 
                P = (7./16) * (48*pow(rho, 2)*pow(z, 5) - 80*pow(rho, 4)*pow(z, 3) + 15*pow(rho, 6)*z)
            elif m==3: 
                P = (7./16) * (80*pow(rho, 3)*pow(z, 4) - 60*pow(rho, 5)*pow(z, 2) + 3*pow(rho, 7))
            elif m==4: 
                P = (7./2) * (10*pow(rho, 4)*pow(z, 3) - 3*pow(rho, 6)*z)
            elif m==5: 
                P = (7./4) * (12*pow(rho, 5)*pow(z, 2) - pow(rho, 7))
            elif m==6: 
                P = 7*pow(rho, 6)*z
            elif m==7: 
                P = pow(rho, 7)
            elif m==8: 
                P = 0.
            else:
                print("error: m must match condition -l-1 < m < l+1")
                return None
            
        else:
            print("l too high")
            return None
        
        return P
    
    def PiI_cart(iB, x, y, z, l, m):
        if iB==0:
            return Harmonic.PiX(x, y, z, l, m)
        elif iB==1:
            return Harmonic.PiY(x, y, z, l, m)
        elif iB==2:
            return Harmonic.PiZ2(x, y, z, l, m)
    
    def PiX(x, y, z, l, m):
            
        if abs(m)>l+1:
            print("error: m must match condition -l-1 < m < l+1")
            return None
            
        if l==0:
            if m==-1:
                P=0
            elif m==0:
                P=0
            elif m==1:
                P=1
        
        elif l==1:
            if m==-2:
                P=y
            elif m==-1:
                P=0
            elif m==0:
                P=-0.5*x
            elif m==1:
                P=z
            elif m==2:
                P=x

        elif l==2:
            if m==-3:
                P=2*x*y
            elif m==-2:
                P=2*y*z
            elif m==-1:
                P= -0.5*x*y
            elif m==0:
                P= -x*z
            elif m==1:
                P= -0.25*(y**2 + 3*x**2 - 4*z**2)
            elif m==2:
                P=2*x*z
            elif m==3:
                P= x**2 - y**2
                
        elif l==3:
            if m==-4:
                P= 3*x**2*y - y**3
            if m==-3:
                P= 6*x*y*z
            elif m==-2:
                P= -0.5*(3*x**2*y + y**3 - 6*y*z**2)
            elif m==-1:
                P= -1.5*x*y*z
            elif m==0:
                P= (3/8)*(y**2*x + x**3 - 4*x*z**2)
            elif m==1:
                P= -0.25*(9*x**2*z + 3*y**2*z - 4*z**3)
            elif m==2:
                P= -x**3 + 3*x*z**2
            elif m==3:
                P= 3*z*(x**2 - y**2)
            elif m==4:
                P= x**3 - 3*x*y**2
                
        else:
            print("l too high")
            return None
        
        return P
            
    def PiY(x, y, z, l, m):
            
        if abs(m)>l+1:
            print("error: m must match condition -l-1 < m < l+1")
            return None
            
        if l==0:
            if m==-1:
                P=1
            elif m==0:
                P=0
            elif m==1:
                P=0
        
        elif l==1:
            if m==-2:
                P=x
            elif m==-1:
                P=z
            elif m==0:
                P=-0.5*y
            elif m==1:
                P=0
            elif m==2:
                P=-y

        elif l==2:
            if m==-3:
                P= x**2 - y**2
            elif m==-2:
                P=2*x*z
            elif m==-1:
                P= -0.25*(x**2 + 3*y**2 - 4*z**2)
            elif m==0:
                P= -y*z
            elif m==1:
                P= -0.5*x*y
            elif m==2:
                P= -2*y*z
            elif m==3:
                P= -2*x*y
                
        elif l==3:
            if m==-4:
                P= x**3 - 3*x*y**2
            if m==-3:
                P= 3*z*(x**2 - y**2)
            elif m==-2:
                P= -0.5*(x**3 + 3*x*y**2 - 6*x*z**2)
            elif m==-1:
                P= -0.25*(3*x**2*z + 9*y**2*z - 4*z**3)
            elif m==0:
                P= (3/8)*(x**2*y + y**3 - 4*y*z**2)
            elif m==1:
                P= -1.5*x*y*z
            elif m==2:
                P= -3*y*z**2 + y**3
            elif m==3:
                P= -6*x*y*z
            elif m==4:
                P= -3*x**2*y + y**3
                
        else:
            print("l too high")
            return None
        
        return P
    
    def PiZ2(x, y, z, l, m):
            
        if abs(m)>l+1:
            print("error: m must match condition -l-1 < m < l+1")
            return None
            
        if l==0:
            if m==-1:
                P=0
            elif m==0:
                P=1
            elif m==1:
                P=0
        
        elif l==1:
            if m==-2:
                P=0
            elif m==-1:
                P=y
            elif m==0:
                P=z
            elif m==1:
                P=x
            elif m==2:
                P=0
                
        elif l==2:
            if m==-3:
                P=0
            elif m==-2:
                P=2*x*y
            elif m==-1:
                P=2*y*z
            elif m==0:
                P= z**2 - 0.5*(x**2+y**2)
            elif m==1:
                P=2*x*z
            elif m==2:
                P= x**2 - y**2
            elif m==3:
                P=0
                
        elif l==3:
            if m==-4:
                P=0
            if m==-3:
                P= 3*x**2*y - y**3
            elif m==-2:
                P=6*x*y*z
            elif m==-1:
                P= 3*y*z**2 - 0.75*(x**2*y+y**3)
            elif m==0:
                P= z**3 - 1.5*z*(x**2+y**2)
            elif m==1:
                P= 3*x*z**2 - 0.75*(y**2*x+x**3)
            elif m==2:
                P= 3*z*(x**2 - y**2)
            elif m==3:
                P= x**3 - 3*x*y**2
            elif m==4:
                P=0
        
        else:
            print("l too high")
            return None
        
        return P
                   
import matplotlib.pyplot as plt
import numpy as np

Rings = {}
with open("p2p_ringspread.txt","r") as f:
    for line in f:
        parse = line.strip().split()
        trial = int(parse[0])
        r = float(parse[1])
        z = float(parse[2])
        Bphi=float(parse[3])
        Brho=float(parse[4])
        Bz  =float(parse[5])

        if trial not in Rings.keys():
            Rings[trial] = [[],[],[],[],[]]
        Rings[trial][0].append(r)
        Rings[trial][1].append(z)
        Rings[trial][2].append(Bphi)
        Rings[trial][3].append(Brho)
        Rings[trial][4].append(Bz)

# Plot Bz vs z
plt.plot(Rings[1963][1],Rings[1963][2],marker='o',linestyle='--',color='blue',label='Standard')
plt.plot(Rings[1964][1],Rings[1964][2],marker='o',linestyle='--',color='red',label='XsYsZc')
plt.plot(Rings[1965][1],Rings[1965][2],marker='o',linestyle='--',color='green',label='XsYcZs')
plt.xticks(fontsize=13)
plt.yticks(fontsize=13)
plt.legend(numpoints=1,loc='best',fontsize=12)
plt.xlabel(r'$z$ (cm)',fontsize=15)
plt.ylabel(r'Spread of $B_\phi$ (pT)',fontsize=15)
plt.ylim([0,1250])
plt.grid(True)
plt.savefig('ringcomparison.pdf',bbox_inches='tight')
plt.show()

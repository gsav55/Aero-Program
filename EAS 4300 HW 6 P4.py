# A rmajet s to propel an aircraft at mach 1 to 6 at hugh altitude and P_amb = 8.5 kPa and T_amb = 220K
# The turbine inlet Temperature is T_max = 2540K
# Construct plots of d(I)/dnb, d(I)/drn, d(TSFC)/dnb, d(TSFC)/drn as a function of M=1 to M=6.
# The derivatives can be approximated using small but discrete changes, such as
#               d(I)/dnb = Delta_I/Delta_nb
# Use baseline nb and rn values of one, and use delta values of 0.01.  Make plots of
# the normalized sensitivity coefficient (e.g. derivitive) as a function of flight mach number.

import math
import matplotlib.pyplot as plt
import numpy as np

# Conditions given in the problem
P_amb = 8.5
T_amb = 220
Gamma = 1.4
R = 287
Cp = 1.005

Ilist = []
TSFClist = []
rnlist = []
nblist = []
dIlist = []
dTSFClist = []
drnlist = []
dnblist = []

def deriv(list):
    if len(list)>2:
        del list[0]
    if len(list)==2:
        return list[0]-list[1]

for M in np.arange(1,6,0.01):
    To_max = 2540
    To4 = To_max
    
    To2 = T_amb*(1+((Gamma-1)/2)*M)
    T_exit = (T_amb/To2)*To_max

    u_exit = M*math.sqrt(Gamma*R*T_exit)
    u_in = M*math.sqrt(Gamma*R*T_amb)

    I = u_exit-u_in
    TSFC = 1/I

    P6 = P_amb
    Po6 = P6*(1+((Gamma-1)/2)*M**2)
    Po4 = Po6
    rn = (Po6/Po4)
    nb = To4-To2
#    np = (2*(u_in/u_out))/(1+(u_in/u_out))
#    nth = ((u_out**2)-(u_in**2))/(2*Cp*(To_max-To_amb))

    Ilist.append(I)
    TSFClist.append(TSFC)
    rnlist.append(rn)
    nblist.append(nb)

    if len(Ilist)>1:
        dIlist.append(deriv(Ilist))
        dTSFClist.append(deriv(TSFClist))
        drnlist.append(deriv(rnlist))
        dnblist.append(deriv(nblist))

        

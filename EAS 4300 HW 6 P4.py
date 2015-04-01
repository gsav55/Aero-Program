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

for M in np.arange(1,6,0.01):
    To_max = 2540

    To2 = T_amb*(1+((Gamma-1)/2)*M)
    T_exit = (T_amb/To2)*To_max

    U_exit = M*math.sqrt(Gamma*R*T_exit)
    U_in = M*math.sqrt(Gamma*R*T_amb)

    np = (2*(u_in/u_out))/(1+(u_in/u_out))
    nth = ((u_out**2)-(u_in**2))/(2*Cp*(To_max-To_amb))

    rho_exit = P_amb/(R*.001*To_max)
    A_ratio = (1/M)*((2/2.4)*(1+(0.4/2)*M**2))**(2.4/0.8)
    m_exit = rho_exit*U_exit*A_ratio
    m_a = m_exit*F

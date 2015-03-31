# -*- coding: utf-8 -*-
# This code is for EAS 4300 Homework 6 Problem 3 and was written by Grayson Savage on 3/30/2015

# Consider an ideal ramjet engine flying at an altitude of 10,000 m (see
# local atmospheric conditions in appendix III of the textbook). The jet fuel used
# has a heat of combustion of 43,000 kJ/kg, and a stoichiometric fuel to air ratio of
# 0.06. The maximum temperature of the ramjet is 2600 K. Assume the specific
# heat ratio up to the combustion chamber is 1.4, and that the specific heat ratio
# through the rest of the engine (including the combustor) is 1.33. The specific heat
# at constant pressure for the gases in the combustor can be calculated using the
# following equation

# Calculate the specific thrust and TSFC as a function of flight Mach number from a range of 1 to
# 6. Note that the fuel to air ratio can not be above the stoichiometric value (i.e. we
# can’t burn more fuel than that at an equivalence ratio of unity). Some operating
# points may require a higher value of f than the stoichiometric value to reach the
# temperature limit. Under this condition, f should be set to the stoichiometric
# value, and the combustor exit temperature T04 should be calculated from the
# energy equation. Otherwise, the maximum temperature is used and f is calculated
# from the energy equation. Using a procedure and an if/then/else approach in EES
# might be ideal for this. Include the following plots:
# 1. Specific thrust (I) vs. Mflight
# 2. TSFC vs. Mflight
# 3. T04 vs. Mflight
# 4. Aexit/Athroat vs. Mflight
# 5. ηth, ηp, and ηo vs. Mflight

import math
import matplotlib.pyplot as plt
import numpy as np

def Velocity(M,Gamma,T,R):
    u = M*math.sqrt(Gamma*R*T)
    return u
def m_flow(density,u,A):
    m_flow = density*u*A
    return m_flow
def TfromTo(To, Gamma, M):
    T = To/(1+((Gamma-1)/2)*M**2)
    return T
def TofromT(T, Gamma, M):
    To = T*(1+((Gamma-1)/2)*M**2)
    return To
def rho(P, R, T):
    density = P/(R*.001*T)
    return density
def FfromCp(To_amb, To_max, h_c, Cp, Fst):
    F = ((To_max/To_amb)-1)/(((h_c)/(Cp*0.001*To_amb))-(To_max/To_amb))
    if F>Fst:
        F = Fst
    return F
def calcCp(Gamma,R):
    Cp = (Gamma/(Gamma-1))*R
    return Cp
def m_spec(m_e, F):
    m_a = m_e/(1+F)
    return m_a
def calcThrust(m_amb, F, u_in, u_out):
    thrust = m_amb*((1+F)*u_out-u_in)
    return thrust
def calcTo_max (T_amb,h_c,F,Cp):
    To_max = (h_c*F - Cp*0.001*T_amb)/(Cp*0.001*F + Cp*0.001)
    return To_max
def propeffic(u_in,u_out):
    np = (2*(u_in/u_out))/(1+(u_in/u_out))
    return np
def thermeffic(u_in, u_out, To_amb, To_max, Cp):
    nth = (.5*((u_out**2)-(u_in**2)))/(Cp*(To_max-To_amb))
    return nth

R = 287
A_exit = 1 # m^2
P_amb = 26.500
T_amb = 223
P_exit = P_amb
h_c = 43000
Fst = 0.06

Mlist = []
Ilist = []
TSFClist = []
To_maxlist = []
nthlist = []
nplist = []
nolist = []
Flist = []
Cplist = []
Tamblist = []
Toamblist = []
Texitlist = []


# Get the range of Mach numbers to iterate the program through
start = float(raw_input('Please enter starting Mach Number: '))
end = float(raw_input('Please enter ending Mach Number: '))
# Calculate needed values based on given values
for M in np.arange(start,end,0.01):
    To_max = 2600
    Gamma = 1.40
    To_amb = TofromT(T_amb, Gamma, M)
#    mach = str(M)
#    print ('Mach Number = ' + mach)
    u_inlet = Velocity(M, Gamma, T_amb, R)
#    inletVel = str(u_inlet)
#    print('Inlet Velocity = ' + inletVel + ' m/s')
    Gamma = 1.33
    Cp = calcCp(Gamma, R)
    F = FfromCp(To_amb, To_max, h_c, Cp, Fst)
    if F == Fst:
        To_max = calcTo_max(To_amb, h_c, F, Cp)
    T_exit = TfromTo(To_max, Gamma,M)
    density_exit = rho(P_exit,R,T_exit)
    u_exit = Velocity(M, Gamma, T_exit, R)
    m_exit = m_flow(density_exit,u_exit,A_exit)
    m_amb = m_spec(m_exit,F)
    thrust = calcThrust(m_amb, F, u_inlet, u_exit)
    I = thrust/m_amb
    TSFC = (F*m_amb)/thrust
    np = propeffic(u_inlet, u_exit)
    nth = thermeffic(u_inlet, u_exit, To_amb, To_max, Cp)
    no = np*nth

    Mlist.append(M)
    Ilist.append(I)
    TSFClist.append(TSFC)
    To_maxlist.append(To_max)
    nthlist.append(nth)
    nplist.append(np)
    nolist.append(no)
    Cplist.append(Cp)
    Flist.append(F)
    Tamblist.append(T_amb)
    Toamblist.append(To_amb)
    Texitlist.append(T_exit)

#Now to plot these values!
plt.figure(1)
plt.figure(1).subplots_adjust(hspace=0.35)
plt.subplot(211)
plt.plot(Mlist, Ilist)
plt.xlabel('Mach Number, M')
plt.ylabel('Specific Thrust, I')
plt.title('Specific Thrust vs Mach Number')

plt.figure(1)
plt.subplot(212)
plt.plot(Mlist, TSFClist)
plt.xlabel('Mach Number, M')
plt.ylabel('TSFC')
plt.title('TSFC vs Mach Number')

plt.figure(2)
plt.figure(2).subplots_adjust(hspace=0.35)
plt.subplot(211)
plt.plot(Mlist, Flist)
plt.xlabel('Mach Number, M')
plt.ylabel('F')
plt.title('F vs Mach Number')

plt.figure(2)
plt.subplot(212)
plt.plot(Mlist, Cplist)
plt.xlabel('Mach Number, M')
plt.ylabel('Cp')
plt.title('Cp vs Mach Number')

plt.figure(3)
plt.figure(3).subplots_adjust(hspace=0.6)
plt.subplot(311)
plt.plot(Mlist, nthlist)
plt.xlabel('Mach Number, M')
plt.ylabel('nth')
plt.title('nth vs Mach Number')

plt.figure(3)
plt.subplot(312)
plt.plot(Mlist, nplist)
plt.xlabel('Mach Number, M')
plt.ylabel('np')
plt.title('np vs Mach Number')

plt.figure(3)
plt.subplot(313)
plt.plot(Mlist, nolist)
plt.xlabel('Mach Number, M')
plt.ylabel('no')
plt.title('no vs Mach Number')

plt.figure(4)
plt.figure(4).subplots_adjust(hspace=0.35)
plt.figure(4).subplots_adjust(wspace=0.30)
plt.subplot(221)
plt.plot(Mlist, To_maxlist)
plt.xlabel('Mach Number, M')
plt.ylabel('To_max')
plt.title('To_max vs Mach Number')

plt.figure(4)
plt.subplot(222)
plt.plot(Mlist, Tamblist)
plt.xlabel('Mach Number, M')
plt.ylabel('T_ambient')
plt.title('T_ambient vs Mach Number')

plt.figure(4)
plt.subplot(223)
plt.plot(Mlist, Toamblist)
plt.xlabel('Mach Number, M')
plt.ylabel('To_ambient')
plt.title('To_ambient vs Mach Number')

plt.figure(4)
plt.subplot(224)
plt.plot(Mlist, Texitlist)
plt.xlabel('Mach Number, M')
plt.ylabel('T_exit')
plt.title('T_exit vs Mach Number')

plt.show()

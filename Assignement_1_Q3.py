# -*- coding: utf-8 -*-
"""
Created on Thu Sep 14 13:03:44 2023

@author: theodoreln
"""

from BEM import *
from Interpolation_of_coeffs import *
from Blade_geometry import *
from Assignement_1_Q1 import *
# Value of lambda_max and pitch_max come from the file of the first question

# Definition of Pmech which is the cap for the Power of the wind turbine
P_mech = 10640000

# We want to make V0 vary until we find V0_rated which correspond to the Power_rated 
""" We want to iterate the computation of P on the wind speed Vo until we reach the right power """

def assignement_1_Q3(lambda_max, CP_max, CT_max, P_mech, R ) :
    # Definition of maximum, minimum and step for the V0
    wind_min=0
    wind_max=25
    wind_step=0.01
    # Initialisation of wind_list and the number of element
    i=int((wind_max-wind_min)/wind_step)
    wind_list=[0]*i
    # Creation of the table for the power P and the omega
    P_table = np.zeros([3,i])
    # Initialisation of the chord for the iteration
    wind = wind_min
    for k in range(i):
        # We store the value of the wind speed
        wind_list[k] = round(wind,2)
        # If V0 < 4m/s we just take omega = 0 because the wind mill is not running
        if wind < 4 :
            P_table[0,k] = 0
            P_table[1,k] = 0
            P_table[2,k] = 0
        else :
            # Computation of omega and P and we check if it stays below P_mech
            omega = lambda_max * wind / R
            P = 0.5*rho*wind*wind*wind*CP_max*pi*R*R
            T = 0.5*rho*wind*wind*CT_max*pi*R*R
            if P < P_mech :
                P_table[0,k] = P
                P_table[1,k] = omega
                P_table[2,k] = T
                omega_max = round(omega,2)
                V0_rated = round(wind,2)
            else :
                P_table[0,k] = P_mech
                P_table[1,k] = omega_max
                P_table[2,k] = 0
        wind = wind + wind_step
    return (V0_rated, omega_max, P_table, wind_list)

""" Results """
V0_rated = 11.43
omega_max = 1.03


""" Computation of the result """
V0_rated, omega_max, P_table, wind_list = assignement_1_Q3(lambda_max, CP_max, CT_max, P_mech, 89.17 )


print('V0_rated = ',V0_rated)
print('omega_max = ',omega_max)
    
plt.plot(wind_list, P_table[0,:], label = 'Power (MW)')
plt.xlabel('Wind speed (m/s)')
plt.legend()
plt.grid()
plt.show()


plt.plot(wind_list, P_table[1,:], label = 'Omega')
plt.axvline(x = V0_rated, color = 'r', label = 'V0_rated = 11.43 m/s')
plt.xlabel('Wind speed (m/s)')
plt.ylabel('Rotational speed (rad/s)')
plt.legend()
plt.grid()
plt.savefig("plot_omega_wind.pdf", format="pdf", bbox_inches="tight")
plt.show()


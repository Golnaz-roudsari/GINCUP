#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 18 10:56:31 2024

@author: roudsari
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 14:30:24 2024

@author: roudsari
"""

import math
import numpy as np


# Params
def get_params(a298=2.8,          #The energy of lateral bonds between the water molecules within a layer
               b=1.5,     # BFHH is related to how rapidly the interaction energy between the surface and the consecutive layers decays     
               rp=550e-9,
               contangice = 72,
               contang = 28,
               ):
    avo = 6.022 * 1e23
    pi = np.pi
    amu = 1.66e-27
    k = 1.381e-23
    Rg = 8.3143
    planck = 6.62607015e-34
    tcrit = 647.096 
    delhvap = 40.66e3 / avo         #Enthalpy of vaporization of water
    phi = contang * pi / 180.0
    phii = contangice * pi / 180.0
    adot = a298 * k * 298
    m = math.cos(phi)
    mi = math.cos(phii)
    elj = 440.3*k
    slj = 0.31e-9
    return avo, pi, amu, k, Rg, planck, tcrit, rp, delhvap, a298, b, contang, contangice, phi, m, phii, mi, adot, elj, slj 


def tparams(t, tc, avo, k, tcrit, amu, a298, b, elj, slj):
    surt = 0.2358 * (1 - t / tcrit) ** 1.256 * (1 - 0.625 * (1 - t / tcrit))        #surface tention = sigma = A(1-t/tc)^n

    roo = 1.8643535 - 0.0725821489 * t + 2.5194368e-3 * t ** 2      #rho = sum_1^10A_i(T/Tc-1)î #density of water               
    roo += -4.9000203e-5 * t ** 3 + 5.860253e-7 * t ** 4 - 4.5055151e-9 * t ** 5
    roo += 2.2616353e-11 * t ** 6 - 7.3484974e-14 * t ** 7 + 1.4862784e-16 * t ** 8
    roo += -1.6984748e-19 * t ** 9 + 8.3699379e-23 * t ** 10
    roo *= 1000

    a0 = 0.9167
    a1 = -1.75e-4
    a2 = -5e-7
    rooi = (a0 + a1 * tc + a2 * tc ** 2) * 1000         #density of ice
    
    vi = amu*18.015/rooi
#    vi = amu * 18.015 / ((0.9060 - 0.14e-3 * tc) * 1000)            # Espinosa et al. 2008 parametrization for molecular volume in ice
    vw = amu * 18.015 / roo
    dw = vw / 10.5e-20          # water monolayer thickness, 10.5 is the water molecule cross section

    tcrit = 647.096  # Critical temperature for water in Kelvin (adjust if different substance)
    pcrit = 22.064e6
    tt = 1 - t / tcrit

    ew1 = (tcrit / t) * (-7.85951783 * tt + 1.84408259 * tt**1.5 - 11.7866497 * tt**3 +
                     22.6807411 * tt**3.5 - 15.9618719 * tt**4 + 1.80122502 * tt**7.5)
    ew1 = pcrit * math.exp(ew1)


    eice = math.exp(9.550426 - 5723.265 / t + 3.53068 * math.log(t) - 0.00728332 * t) #Saturation vapor pressure of ice
    
#    stsl = (29.986 + 0.25559 * tc - 0.0010465 * tc ** 2 - 4.6503e-6 * tc ** 3 + 2.9065e-7 * tc ** 4) * 1e-3     #ice-water interfacial tension, Espinosa et al. 2018
    stsl = 0.1364-0.15*t*1e-3 -surt                  #hex Young with Plummer & Hale T-dep.
    delmu = (0.00035032 - 0.0046013 * tc - 2.3187e-5 * tc**2 + 6.9536e-8 * tc**3) / avo * 4.184e3 #Chemical potential difference between ice and water, Espinosa et al. 2018

#    rcfr2 = 2 * vi * stsl / (k * t * math.log(ew1 / eice))      #critical size for freezing hexagonal water ice
    rcfr2 = 2*vi*stsl/(k*t*np.log(ew1/eice))

    vfr = 4 / 3 * np.pi * rcfr2 ** 3


#    aa = a298 * k * 298.15  #temperature independent (A')
    aa = a298*k*298.15
    a = aa/k/t            # A(T)
    di = 3.63e-10 #Thickness of one monolayer of ice
#    apu = vw**2/dw**3*(1/vi-1/vw)
    aice = a * (dw / di)**b + np.pi * elj * slj**(b+3) / (3*k * t * di**b) * (vw**(-b / 3) - vi**(-b / 3))
#    aice = a*(dw/di)**b



    return surt, roo, rooi, vw, dw, eice, vi, ew1, stsl, rcfr2, vfr, a, aice

def clustads(t, r, a, surt, dw, vw, k, tcrit, phi, m, phii, mi, b, rp):
    # Constants
    pi = np.pi

    # Calculations
    d = np.sqrt(rp**2 + r**2 - 2.0 * rp * r * m) #distance
    x = r + d - rp  # Max. thickness of ice cluster if ice contact angle >= water contact angle
    cosp = -(r - rp * m) / d #
    cosf = (rp - r * m) / d #Laaksonen et. el. 2020
    v = pi / 3.0 * (r ** 3 * (2 - 3 * cosp + cosp ** 3) - rp ** 3 * (2 - 3 * cosf + cosf ** 3))
    ala = 2 * pi * rp ** 2 * (1 - cosf) #Area of liquid-air interface ala

    fii = math.acos((rp - r * m) / d) #

    if phii < phi:  # Max. thickness of ice cluster if ice contact angle < water contact angle
        u = (rp - r * m) / d
        p = u ** 2 - mi ** 2
        q = 2 * rp * mi * (1 - u * u)
        s = rp ** 2 * (u * u - 1)
        rice = (-q + (q ** 2 - 4 * p * s) ** 0.5) / (2 * p)
        dice = (rp * rp + rice * rice - 2 * rp * rice * mi) ** 0.5
        x = rice + dice - rp

    step = phi / 500.0 #  Initialize step size for gamma, based on contact angle
    gamma = phi - step / 2.0 # set initialize value of gamma, which is used in the loop
    alfaold = fii # angle phi in classical nuc book page 140
    hsum = 0 # Initialize hsum which will accumulate the sum in the loop

    for j17 in range(1, 501):  # Equilibrium saturation ratio for clusterwise adsorption
        h = -rp + r * math.cos(gamma) + math.sqrt(r ** 2 * (math.cos(gamma)) ** 2 - r ** 2 + d ** 2)     # Calculate the height h using the current gamma value

        gamma2= gamma-step/2 # Calculate gamma2, which is gamma offset by half a step

        h2 = -rp + r * math.cos(gamma2) + math.sqrt(r ** 2 * (math.cos(gamma2)) ** 2 - r ** 2 + d ** 2) # Calculate the height h2 using gamma2

        alfa2_i= (rp + h2 - r * math.cos(gamma2)) / d # Calculate the intermediate alfa2_i
#        print(alfa2_i)

#        if alfa2_i > 1: # the Python version is double-percision 
#           alfa2_i = 1  #it gets to 1
            
        alfa2 = np.arccos((alfa2_i))  # Calculate the angle alfa2 from the intermediate alfa2_i

        cosa2 = (rp + h2 - r * math.cos(gamma2)) / d # Calculate cosa2 from alfa2_i

        hsum += h ** b * 2 * pi * rp ** 2 * (cosa2 - math.cos(alfaold))     # Update hsum by adding the current term

        alfaold = alfa2 # Update alfaold for the next iteration

        gamma = gamma - step # Decrease gamma by the step size for the next iteration

    

    theta1 = hsum / (2 * pi * rp ** 2 * (1 - cosf))  # Calculate the thickness of ice 
    thetab = theta1 / dw ** b # N^B,      theta1= sigma^B  Laaksonen_2020_acp

    s1 = np.exp(-a / thetab + 2 * surt * vw / (k * t * r)) # Calculate s1, the equilibrium saturation ratio, using thetab
    theta = thetab ** (1 / b) #theta=N number of monolayers
    

    if s1 > 1:
       s1 = 1
    return x, v, ala, theta, s1

def homrate(t, tc, rcfr2, theta, stsl, adot, delhvap, k, b, avo, amu, sice, aice, a, ew1, eice, vi):

    

    aesp = math.exp(91.656 + 0.11729 * tc - 0.00081401 * tc ** 2) #kinetic prefactor
    delmu = (0.00035032 - 0.0046013 * tc - 2.3187e-5 * tc ** 2 + 6.9536e-8 * tc ** 3) / avo * 4.184e3
    vii = amu * 18.015 / ((0.9060 - 0.14e-3 * tc) * 1000)
    stsle = (29.986 + 0.25559 * tc - 0.0010465 * tc**2 - 4.6503e-6 * tc**3 + 2.9065e-7 * tc**4) * 1e-3
    delg = 16 * np.pi * vii ** 2 * stsle ** 3 / (3 * delmu ** 2) # Delta G in pure water


        
    sigw = 10.5e-20
    sigi = vii/3.63e-10
    delmuad = delmu+ k*t*theta**(-b)*(aice*(sigw/sigi)**(b)-a)
    delmuadr = k*t*np.log(ew1/eice) + k*t*theta**(-b)*(aice*(sigw/sigi)**(b)-a)
    rcfr = rcfr2*(k*t*np.log(ew1/eice)/delmuadr)
    nrw = aesp*np.exp(-delg/(k*t))
    expcom = -16*np.pi*stsle**3*vii**2/(k*t*3)

    jw = np.exp(expcom/delmu**2)


    ja = np.exp(expcom/delmuad**2)

    jratio = ja/jw
    nra = nrw*jratio
    return np.float64(nra), np.float64(nrw), np.float64(rcfr), np.float64(jratio)


def icevapor(t, tc, a, ew1, aice, rind, amu, rp, m, phi, b, k):
    surt = 0.1364 - 0.15 * t * 1e-3
    vi = amu * 18.015 / ((0.9060 - 0.14e-3 * tc) * 1000)
    dw = 3.63e-10
    if rind <= rp:
        r = 0.1e-10
        s = 0
        for i1 in range(1, 15001):
            r += 1e-10
            d = (rp * rp + r * r - 2 * rp * r * m) ** 0.5
            cosp = -(r - rp * m) / d
            cosf = (rp - r * m) / d
            step = phi / 100
            gamma = phi - step / 2
            alfaold = math.acos((rp - r * m) / d)
            hsum = 0
            for j17 in range(1, 101):
                h = -rp + r * math.cos(gamma) + math.sqrt(r ** 2 * (math.cos(gamma)) ** 2 - r ** 2 + d ** 2)
                gamma2 = gamma - step / 2
                h2 = -rp + r * math.cos(gamma2) + math.sqrt(r ** 2 * (math.cos(gamma2)) ** 2 - r ** 2 + d ** 2)
                alfa2 = math.acos((rp + h2 - r * math.cos(gamma2)) / d)
                cosa2 = (rp + h2 - r * math.cos(gamma2)) / d
                hsum += h ** b * 2 * np.pi * rp ** 2 * (cosa2 - math.cos(alfaold))
                alfaold = alfa2
                gamma = gamma - step
            theta1 = hsum / 2 / np.pi / rp ** 2 / (1 - cosf)
            theta = theta1 / dw ** b
            s1 = math.exp(-aice / theta + 2 * surt * vi / (k * t * r))
            eice = math.exp(9.550426 - 5723.265 / t + 3.53068 * math.log(t) - 0.00728332 * t)
            sw = s1 * (eice / ew1)
            if sw > 1:
                sicevmax = ew1 / eice
                return sicevmax
            if s1 < s:
                sicevmax = s
                rmax = r
                thetamax = theta
                return sicevmax
            else:
                s = s1
    else:
        r = rp + 1e-10
        s = 0
        for i2 in range(1, 80001):
            r += 0.1e-10
            theta = (r - rp) / dw
            s1 = math.exp(-aice / theta ** b + 2 * surt * vi / (k * t * r))
            if s1 < s:
                sicevmax = s
                rmax= r
                thetamax = theta
                return sicevmax
            else:
                s = s1
        
        return sicevmax         

def filmads(t, tc, r, a, surt, dw, vw, avo, amu, k, Rg, planck, tcrit, a298, b, rp ):
    pi = math.pi


    # Calculations
    theta = (r - rp) / dw
    s1 = math.exp(-a / theta ** b + 2 * surt * vw / (k * t * r))
    x = r - rp
    v = 4.0 / 3.0 * pi * (r ** 3 - rp ** 3)
    ala = 4.0 * pi * rp ** 2

    # Check condition
    if s1 > 1.0:
        s1 = 1.0

    return x, v, ala, theta, s1

def get_constants():
    avo = 6.022 * 1e23
    pi = np.pi
    amu = 1.66e-27
    k = 1.381e-23
    Rg = 8.3143
    planck = 6.62607015e-34
    tcrit = 647.096 
    rp = 50e-9
    delhvap = 40.66e3 / avo       
    return avo, pi, amu, k, Rg, planck, tcrit, rp, delhvap










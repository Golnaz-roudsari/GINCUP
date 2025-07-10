#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  1 17:01:25 2024

@author: roudsari
"""

import numpy as np
from HomFilmAlgorithm import get_params, tparams, icevapor, homrate, filmads


def HomFilmSim(a_i, b_i, rp_i, output_file): 

    avo, pi, amu, k, Rg, planck, tcrit, rp, delhvap, a298, b, contang, contangice, phi, m, phii, mi, adot, elj, slj = get_params(a298 = a_i,
                                                                                                                       b = b_i,
                                                                                                                      rp = rp_i )
    ti = 273.15 - 1

    with open(output_file, 'w') as file:
        for j1 in range(1, 74):
       
            tci = ti - 273.15
            surti, roo, rooi, vw, dwi, eicei, vicei, ew1i, stsli, rcfr2i, vfri, a, aice = tparams(ti, tci, avo, k, tcrit, amu, a298, b, elj, slj)
            
            n = 75e-9/dwi
            
            ri = rp + rcfr2i*0.05 #for filmwise adsoption
            

            for i1 in range(1, 30001):
                ri += 0.1e-10
                xi, vi, alai, thetai, s1i = filmads(ti, tci, ri, a, surti, dwi, vw, avo, amu, k, Rg, planck, tcrit, a298, b, rp )


                sicei = s1i * (ew1i / eicei)

                if sicei <= 1:
                    continue
                
                nrhoi, nrw, rcfri, jratioi = homrate(ti, tci, rcfr2i, thetai, stsli, adot, delhvap, k, b, avo, amu, sicei, aice, a, ew1i, eicei,vi)
                
                
                dcrit = (rp ** 2 + rcfri ** 2 - 2 * rp * rcfri * mi) ** 0.5
                xcrit = rcfri + dcrit - rp 
                nprobhe = 1 - np.exp(-nrhoi * 10 * alai) 
                
                if s1i >= 1:
                    break

                vaddi = (4.0 / 3.0) * np.pi * ((rp + 2 * dwi)**3 - rp**3)

                if (xi / dwi >= 2 and xi >= 2 * rcfri + 3 * dwi and nrhoi * (vi - vaddi) >= 0.001):
                    sicevmaxi = icevapor(ti, tci, a, ew1i, aice, ri, amu, rp, m, phi, b, k)
                    if sicevmaxi>sicei:
                        break
                    else:
                        break
            
            file.write(f"{ti}, {sicei*100},{nrhoi}, {nrw}, {rcfri}, {jratioi}\n")
            ti -= 1

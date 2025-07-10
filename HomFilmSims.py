#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  1 21:11:28 2024

@author: roudsari
"""

import numpy as np
from HomFilm import HomFilmSim

#path1 = './HomFilmResults/RH_results/RH_a_b_results/RH_original/'
#path2 = './HomFilmResults/RH_results/RH_a_b_results/RH_NewEw/'
path2 = './HomFilmResults/NewResultsSilica/'

a = 3.8
b = 1.2
rp = 200
dw = 2

file_path =  path2 + 'homfilm' +'_a' + str(a) + '_b' + str(b) + '_rp' + str(rp) + '_ML'+ str(dw)+'.txt'  
#file_path =  path2 + 'homfilm' +'water_saturation'+'.txt'              
HomFilmSim(a, b, rp*1e-9, file_path)

            
import pandas as pd
import matplotlib.pyplot as plt

def plot_columns(file_path):

    df = pd.read_csv(file_path, delimiter=',', header=None)
    

    df.columns = ['col1', 'col2', 'col3', 'col4', 'col5', 'col6']
    

    df.dropna(subset=['col1', 'col2'], inplace=True)
    
    # Plot the second column vs the first column
    plt.figure(figsize=(10, 6))
    plt.plot(df['col1'], df['col2'], '.-', label='Sice')
    plt.xlabel('Column 1')
    plt.ylabel('Column 2')
    plt.title('Column 2 vs Column 1 Plot')
    plt.legend()
    plt.grid(True)
    plt.show()

# Example usage
plot_columns(file_path)

            
            
            
            


# INPUT-HomFilm
This repository provides the Python implementation for simulating homogeneous filmwise ice nucleation in the upper troposphere. 
This work supports the findings in our publication:
Homogeneous ice nucleation in adsorbed water films: A theoretical approach
Authors: Ari Laaksonen, Golnaz Roudsari, Ana A. Piedehierro, and André Welti
Published in: Atmospheric Physics and Chemistry (ACP), 2025, https://doi.org/10.5194/egusphere-2024-4095, 2025

## 👥 Authors:
Ari Laaksonen - [@ajlaaksonen] (https://github.com/ajlaaksonen)


The INPUT code calculates nucleation rates and saturation ratios by simulating molecular interactions in supercooled water films. It captures detailed physical parameters like surface tension, density, and vapor pressures over temperature ranges.

This code is designed to reproduce the results shown in our published paper. 
- Computes critical cluster size and nucleation probabilities.
- Implements adsorption, homogeneous nucleation rate models, and film growth.
- Includes plotting routines for visualization.

## 📦 Requirements
  Python 3.x  
  Libraries:
    - numpy
    - matplotlib
    - pandas



## ▶️ How to Run:
python HomFilmSims.py
Modify parameters a, b, and rp in HomFilmSims.py to suit your material or system.

## 🧾 Citation
If you use this code, please cite:
Laaksonen, A., Roudsari, G., Piedehierro, A. A., and Welti, A.: Homogeneous ice nucleation in adsorbed water films: A theoretical approach, EGUsphere [preprint], https://doi.org/10.5194/egusphere-2024-4095, 2025.



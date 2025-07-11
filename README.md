# GINCUP
This repository provides the Python implementation for simulating homogeneous filmwise ice nucleation in the upper troposphere. 
This work supports the findings in our publication:
Homogeneous ice nucleation in adsorbed water films: A theoretical approach
Authors: Ari Laaksonen, Golnaz Roudsari, Ana A. Piedehierro, and André Welti
Published in: Atmospheric Physics and Chemistry (ACP), 2025, https://doi.org/10.5194/egusphere-2024-4095, 2025

## 👥 Authors:
Ari Laaksonen - [@ajlaaksonen] (https://github.com/ajlaaksonen)
Golnaz Roudsari  - [@Golnaz-roudsari]


This repository provides code to simulate **homogeneous ice nucleation in multilayer water films** adsorbed on insoluble atmospheric particles. It is based on the **Frenkel–Halsey–Hill (FHH) adsorption theory** combined with **classical nucleation theory (CNT)**.

The model explains how the freezing temperature, critical cluster size, and nucleation rates vary with film thickness and substrate interaction strength. It supports experimental validation using silica particles and is designed to capture deposition ice nucleation behavior in the upper troposphere. We will develop this model for heterogenous and droplet wise nucleation


## 📊 Features
- Computes nucleation rate `J` and critical nucleus radius as a function of temperature
- Models melting point depression in nanoscale water films
- Simulates relative humidity thresholds (`RH_ice*`) for ice nucleation
- Compares theoretical predictions with SPIN experiment data on silica
- Includes plotting functionality for results visualization

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



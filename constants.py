#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  5 10:33:25 2024

@author: yogehs
"""

#Value of thermal energy at room temperature (25deg C) in pN nm
kbT = 4.114  # pn nm
k_bT_nm = 4.114

mode_roll = 'reflect'

#limits for plotting 
L_p_min_lt, L_p_max_lt = 0, 95
y_min, y_max = 1, 1000
x_renorm_limits = [0.030644186650635584, 4996.008916478943]
y_renorm_limit = [0.06738070130927563, 2.9694807353946584]

x_un_limits = [0.23042308789504504, 78.19536741138113]
y_un_limits = [0.7854736945816706, 148.95038821847717]
y_min_norm, y_max_norm = 0.00008, 0.8
x_bend_lt = [0.23042308789504504, 78.19536741138113]
y_bend_lt = [10, 6000]

fit_l_p_header = ['file_label', 'config_label', 'fil_label', 'avg_arc_L',
                  'L_fitted', 'L_P_fitted', 'zeta_fitted', "avg_L_p_WLC", 'tau_fitted', 'r_c']
#dict of fps to convert the frames to time in seconds for each configuration 
long_vids_fps = {"S_1": 1}
#dict of hex color values for each configuration
color_saha_configs_hex = {'S': '#ff0000'}
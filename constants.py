#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  5 10:33:25 2024

@author: yogehs
"""

L_p_min_lt, L_p_max_lt = 0, 85

kbT = 4.114  # pn nm
k_bT_nm = 4.114
mode_roll = 'reflect'
y_min, y_max = 1, 1000

x_renorm_limits = [0.030644186650635584, 4996.008916478943]
y_renorm_limit = [0.06738070130927563, 2.9694807353946584]

x_un_limits = [0.23042308789504504, 78.19536741138113]
y_un_limits = [0.7854736945816706, 148.95038821847717]
y_min_norm, y_max_norm = 0.00008, 0.8

fit_l_p_header = ['file_label', 'config_label', 'fil_label', 'avg_arc_L',
                  'L_fitted', 'L_P_fitted', 'zeta_fitted', "avg_L_p_WLC", 'tau_fitted', 'r_c']


long_vids_fps = {"C_1_long": 1, "C_2_long": 1, "C_3_long": 1, "C_4_long": 1, "C_5_long": 3.33, "C_6_long": 0.5, 'C_7_long': 1,

                 "S_1": 1, "S_2": 1, "S_3": 1, "S_4": 1, "S_5": 1, "S_6": 1,
                 "S_7": 1, "S_8": 2, "S_9": 2,


                 "S_1_long": 1, "S_2_long": 1, "S_3_long": 1, "S_4_long": 3.33, "S_5_long": 2, "S_6_long": 2,
                 'ring_noATP_1': 2, 'ring_noATP_2': 2, 'ring_noATP_3': 2,

                 'ring_ATP_1': 1, 'ring_ATP_2': 1,
                 'ring_ATP_3': 1, 'ring_ATP_4': 1,


                 'mixed_1_long': 0.5, 'mixed_2_long': 1.3, 'HMRN_H_closed_1': 1,
                 'HMRN_H_closed_2': 1,
                 'HMRN_H_closed_3': 1,
                 'HMRN_H_closed_4': 1,
                 'HMRN_H_open_1': 1,
                 'HMRN_H_open_2': 1,

                 'braided_1': 1, 'braided_2': 1,
                 'braided_3': 1,

                 'monomer_1': 1, 'monomer_2': 2,
                 'monomer_3': 3.33, 'monomer_4': 2, 'monomer_5': 1,


                 'vol_diff_S_1': 2, 'vol_diff_S_2': 1, 'vol_diff_S_3': 1, 'vol_diff_S_4': 1, 'vol_diff_S_5': 3.33,
                 'hRm_open_1': 1, 'hRm_open_2': 1,

                 'hMRm_closed_1': 1, 'hMRm_closed_2': 1, 'hMRm_closed_3': 1, 'hMRm_closed_4': 1,

                 'sbcC_ATP_1': 1, 'sbcC_ATP_2': 1, 'sbcC_ATP_3': 1, 'sbcC_ATP_4': 1, 'sbcC_ATP_5': 1,
                 'hMR_closed_1': 1, 'hMR_closed_2': 1, 'hMR_closed_3': 1,
                 'hMR_closed_4': 1, 'hMR_closed_5': 1, 'hMR_closed_6': 1,

                 'hMR_open_1': 1, 'hMR_open_2': 1, 'hMR_open_3': 1,
                 'hMR_open_4': 1,

                 'hMRm_closed_1': 1, 'hMRm_closed_2': 1, 'hMRm_closed_3': 1,
                 'hMRm_open_1': 2, 'hMRm_open_2': 2, 'hMRm_open_3': 1,

                 'hR_open_1': 1, 'hR_open_2': 0.5,

                 'hRm_closed_1': 1, 'hRm_closed_2': 1, 'hRm_closed_3': 1,

                 'hRm_open_1': 1, 'hRm_open_2': 1,
                 }
color_saha_configs_hex = {'S': '#ff0000',
                          'C': '#FF00FF',

                          'closed': '#008000',
                          'open': '#00bfff',
                          'b': '#9400d3',


                          'm': '#000000',
                          'noATP': '#4169e1',
                          'ATP': '#0000ff',

                          'v': '#980000',
                          'x': '#084a91',
                          'y': '#006428ff',
                          's': '#ff4500',

                          "hMR_closed": '#556B2F', 'hMR_open': '#000080',
                          'hMRm_closed': '#32CD32', "hMRm_open": '#00CED1',

                          'hR_open': '#7d5f04', "hRm_open": '#FFAA33',
                          "hRm_closed": '#C87137',
                          'parallel': '#084a91', "perpendicular": '#006428ff'

                          }